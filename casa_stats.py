import sys
import os
import json
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Ellipse
from argparse import ArgumentParser

from casatasks import imstat, imhead
import casatools as ct

def prune_mask(mask, start_x, start_y):
    '''
    Take mask and only retain the region containing
    the starting x and y coordinates
    '''
    sys.setrecursionlimit(mask.size)

    def find_nearest_nonzero(matrix, x, y):
        idx = np.where(matrix > 0)

        x_off = idx[0]-x
        y_off = idx[1]-y

        min_offset = np.argmin(x_off**2+y_off**2)

        return idx[0][min_offset], idx[1][min_offset]

    def floodfill(matrix, x, y):
        # Flood fill algorithm with 8-connectivity
        if matrix[x,y] == 1:
            matrix[x,y] = 2
            if x > 0:
                floodfill(matrix,x-1,y)
            if x < matrix.shape[0] - 1:
                floodfill(matrix,x+1,y)
            if y > 0:
                floodfill(matrix,x,y-1)
            if y < matrix.shape[1] - 1:
                floodfill(matrix,x,y+1)
            if x > 0 and y > 0:
                floodfill(matrix,x-1,y-1)
            if x > 0 and y < matrix.shape[1] - 1:
                floodfill(matrix,x-1,y+1)
            if y > 0 and x < matrix.shape[0] - 1:
                floodfill(matrix,x+1,y-1)
            if x < matrix.shape[0] - 1 and y < matrix.shape[1] - 1:
                floodfill(matrix,x+1,y+1)

    try:
        mask[start_x, start_y]
    except IndexError:
        print('Given coordinates are outside the image, falling back to default')
        start_x = int(mask.shape[0]/2)
        start_y = int(mask.shape[1]/2)

    # If starting position not in mask, find nearest pixel that is
    if mask[start_x, start_y] == 0:
        start_x, start_y = find_nearest_nonzero(mask, start_x, start_y)

    floodfill(mask, start_x, start_y)
    mask[mask == 1] = 0
    mask[mask == 2] = 1

    return mask

def include_ellipses(gaussians, wcs):
    '''
    If Gaussians are specified, define matplotlib ellipses

    Keyword arguments:
    guassians -- List of Gaussians to turn to matplotlib ellipses
    wcs -- Coordinate system of the image
    '''
    if 'Gaus_id' in gaussians.colnames:
        ids = gaussians['Gaus_id']
    else:
        ids = gaussians['Source_id']
        gaussians = [gaussians]

    # Get ids
    increment = max(wcs.increment(format='n', type='direction')['numeric'])

    ell = []
    for source in gaussians:
        xy = wcs.topixel([np.deg2rad(source['RA']), np.deg2rad(source['DEC'])])
        ell.append(Ellipse(xy = xy['numeric'][:2],
                           width = source['Min']/np.rad2deg(increment),
                           height = source['Maj']/np.rad2deg(increment),
                           angle = source['PA']))

    return ell, ids

def transform_axis(axis_len, pix_size):
    '''
    Get correct axes (in arcsec)
    '''
    interval = 10
    nticks = np.inf
    while nticks > 10:
        image_size = axis_len*pix_size
        nearest_round = int((image_size/2-1)/interval)*interval

        labels = np.arange(-nearest_round, nearest_round+interval, interval, dtype=int)
        ticks = (labels+image_size/2)/pix_size

        nticks = len(ticks)
        interval += 10

    return ticks, labels

class Image:

    def __init__(self, input_image, threshold, start_coordinates=None, source_id=0):
        self.image_file = input_image
        self.threshold = threshold
        self.source_id = source_id

        # Open image and get needed data
        ia = ct.image()

        ia.open(input_image)
        self.image = np.squeeze(ia.getregion())
        self.wcs = ia.coordsys([0,1])

        if start_coordinates:
            # Get the maximum position
            max_ra = np.deg2rad(float(start_coordinates[0]))
            max_dec = np.deg2rad(float(start_coordinates[1]))

            pixel_coord = self.wcs.topixel([max_ra, max_dec])
            self.start_x = int(pixel_coord['numeric'][0])
            self.start_y = int(pixel_coord['numeric'][1])
        else:
            self.start_x = -1
            self.start_y = -1
        ia.close()

        # Initialize central coordinates
        self.x_mean = 0
        self.y_mean = 0

        # Get beam parameters
        self.bmaj = imhead(input_image, mode='get', hdkey='BMAJ')['value']
        self.bmin = imhead(input_image, mode='get', hdkey='BMIN')['value']
        self.bpa = imhead(input_image, mode='get', hdkey='BPA')['value']

    def plot_image(self, mask_2d, gaussians=None, source=None):

        # Get pixel size in arcsec
        increment = max(self.wcs.increment(format='n', type='direction')['numeric'])
        pix_size = np.rad2deg(increment)*3600

        # Plot image
        fig = plt.figure()
        fig.add_subplot(111)

        cmap_norm = TwoSlopeNorm(vcenter=self.threshold,
                                 vmax=np.max(self.image[mask_2d.astype(bool)]))
        plt.imshow(self.image.T, interpolation='none', origin='lower',
                   norm=cmap_norm, cmap='jet')
        plt.contour(mask_2d.T, colors='k', levels=[0,1], linewidths=1)
        plt.scatter(self.x_mean, self.y_mean, color='w', marker='x', s=5, zorder=2)

        # Add beam
        pix_bmaj = self.bmaj/pix_size #arcsec to deg
        pix_bmin = self.bmin/pix_size #arcsec to deg
        beam = Ellipse(xy=(pix_bmaj/2,pix_bmaj/2),
                       width=pix_bmin,
                       height=pix_bmaj,
                       angle=self.bpa)

        plt.gca().add_artist(beam)
        beam.set_facecolor('white')
        beam.set_edgecolor('black')

        if gaussians:
            ell, ids = include_ellipses(gaussians, self.wcs)
            for i, e in enumerate(ell):
                plt.gca().add_artist(e)
                e.set_facecolor('none')
                e.set_edgecolor('m')
                e.set_lw(1)

                # Get ellipse vertices to annotate the id at the correct location
                path = e.get_path()
                vertices = path.vertices.copy()
                vertices = e.get_patch_transform().transform(vertices)
                xy = (vertices[0][0]-0.5, vertices[0][1]+0.5)
                plt.annotate(text=ids[i],
                             xy=xy,
                             color='m')

        if source:
            ell, source_id = include_ellipses(source, self.wcs)
            for i, e in enumerate(ell):
                plt.gca().add_artist(e)
                e.set_facecolor('none')
                e.set_edgecolor('b')
                e.set_lw(1)

                # Get ellipse vertices to annotate the id at the correct location
                path = e.get_path()
                vertices = path.vertices.copy()
                vertices = e.get_patch_transform().transform(vertices)
                xy = (vertices[0][0]-0.5, vertices[0][1]+0.5)
                plt.annotate(text=source_id,
                             xy=xy,
                             color='b')

        xticks, xlabels = transform_axis(self.image.T.shape[1], pix_size)
        yticks, ylabels = transform_axis(self.image.T.shape[0], pix_size)

        plt.xticks(xticks,xlabels)
        plt.yticks(yticks,ylabels)

        plt.xlabel('dRA (arcsec)')
        plt.ylabel('dDEC (arcsec)')
        plt.savefig(self.image_file.split('.')[0]+'.png')
        plt.close()

    def plot_alpha_image(self, alpha_image, mask_2d, gaussians=None, source=None):

        # Get pixel size in arcsec
        increment = max(self.wcs.increment(format='n', type='direction')['numeric'])
        pix_size = np.rad2deg(increment)*3600

        # Open alpha image
        ia = ct.image()

        ia.open(alpha_image)
        image = np.squeeze(ia.getregion())
        ia.close()

        fig = plt.figure()
        fig.add_subplot(111)

        plt.imshow(image.T, interpolation='none', origin='lower',
                   vmin=-3, vmax=3, cmap='coolwarm')
        plt.contour(mask_2d.T, colors='k', levels=[0,1], linewidths=1)
        plt.scatter(self.x_mean, self.y_mean, color='w', marker='x', s=5, zorder=2)

        # Add beam
        increment = max(self.wcs.increment(format='n', type='direction')['numeric'])

        pix_bmaj = self.bmaj/np.rad2deg(increment)/3600 #arcsec to deg
        pix_bmin = self.bmin/np.rad2deg(increment)/3600 #arcsec to deg
        beam = Ellipse(xy=(pix_bmaj/2,pix_bmaj/2),
                       width=pix_bmin,
                       height=pix_bmaj,
                       angle=self.bpa)

        plt.gca().add_artist(beam)
        beam.set_facecolor('black')
        beam.set_edgecolor('white')

        if gaussians:
            ell, ids = include_ellipses(gaussians, self.wcs)
            for i, e in enumerate(ell):
                plt.gca().add_artist(e)
                e.set_facecolor('none')
                e.set_edgecolor('m')
                e.set_lw(1)

                # Get ellipse vertices to annotate the id at the correct location
                path = e.get_path()
                vertices = path.vertices.copy()
                vertices = e.get_patch_transform().transform(vertices)
                xy = (vertices[0][0]-0.5, vertices[0][1]+0.5)
                plt.annotate(text=ids[i],
                             xy=xy,
                             color='m')

        if source:
            ell, source_id = include_ellipses(source, self.wcs)
            for i, e in enumerate(ell):
                plt.gca().add_artist(e)
                e.set_facecolor('none')
                e.set_edgecolor('g')
                e.set_lw(1)

                # Get ellipse vertices to annotate the id at the correct location
                path = e.get_path()
                vertices = path.vertices.copy()
                vertices = e.get_patch_transform().transform(vertices)
                xy = (vertices[0][0]-0.5, vertices[0][1]+0.5)
                plt.annotate(text=source_id,
                             xy=xy,
                             color='b')

        xticks, xlabels = transform_axis(self.image.T.shape[1], pix_size)
        yticks, ylabels = transform_axis(self.image.T.shape[0], pix_size)

        plt.xticks(xticks,xlabels)
        plt.yticks(yticks,ylabels)

        plt.xlabel('dRA (arcsec)')
        plt.ylabel('dDEC (arcsec)')
        plt.savefig(alpha_image.split('.')[0]+'.png')
        plt.close()

    def get_stats(self, plot=False, source=None, gaussians=None, alpha_image=None):
        '''
        Get the stats of the source in an image

        Keyword arguments:
        input_image (string) -- Filename of the input image
        threshold (float) -- Threshold for stats (in general in units of Jy/beam)
        coord (tuple) -- Starting coordinates (RA, DEC) to identify correct island
        source_id (string) -- Source id to associate source with
        plot (bool) -- Output plot with the source, threshold, and optionally Gaussians
        gaussians (table object) -- List of Gaussians associated with the source
        alpha_image (string) -- If specified, measure spectral index from this file
        '''
        im = ct.imager()
        ia = ct.image()

        # Create mask using the threshold
        im.mask(self.image_file, mask='stats_mask', threshold=self.threshold)

        # Open mask and retain only pixels connected to the source
        ia.open('stats_mask')
        mask = ia.getregion()

        mask_2d = mask[:,:,0,0]
        mask_2d = prune_mask(mask_2d, self.start_x, self.start_y)

        mask[:,:,0,0] = mask_2d
        ia.putregion(mask)
        ia.close()

        # Run imstat to get integrated flux
        image_stats = imstat(self.image_file, mask='stats_mask')

        # Open image and get center of mass
        valid_pix = np.where(mask_2d)
        self.x_mean = np.average(valid_pix[0], weights=self.image[valid_pix])
        self.y_mean = np.average(valid_pix[1], weights=self.image[valid_pix])

        center_of_mass = self.wcs.toworld([self.x_mean,self.y_mean], 'n')
        ra_deg = np.rad2deg(center_of_mass['numeric'][0])
        dec_deg = np.rad2deg(center_of_mass['numeric'][1])

        is_inmask = bool(mask_2d[int(self.x_mean),int(self.y_mean)])
        is_maxpos = ((self.start_x - image_stats['maxpos'][0])**2
                     + (self.start_y - image_stats['maxpos'][1])**2)**0.5 < 2.

        os.system('rm -r stats_mask')

        if alpha_image:
            # Get spectral index value from image
            ia.open(alpha_image)
            alpha = np.squeeze(ia.getregion())
            mean_alpha = np.average(alpha[valid_pix], weights=self.image[valid_pix])
            ia.close()

        if plot:
            self.plot_image(mask_2d, gaussians, source)
            if alpha_image:
                self.plot_alpha_image(alpha_image, mask_2d, gaussians, source)

        source_attr = {'Source_id':         self.source_id,
                       'RA_mean':           ra_deg,
                       'DEC_mean':          dec_deg,
                       'Cutout_Total_flux': image_stats['flux'][0],
                       'Isinmask':          is_inmask,
                       'Ismaxpos':          is_maxpos}
        if gaussians:
            source_attr['N_Gaus'] = len(gaussians)
        if alpha_image:
            source_attr['Cutout_Spectral_index'] = mean_alpha

        return source_attr

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    input_image = args.in_image
    threshold = args.threshold
    source_id = args.id
    coord = args.coord
    plot = args.plot

    image = Image(input_image, threshold, coord, source_id)
    source_attr = image.get_stats(plot)

    print(f'Measured properties for source with source id {source_id} in {input_image}: \n')
    print(f"Total flux above threshold ({float(threshold):.2e}): {source_attr['Cutout_Total_flux']}")
    print(f"Mean RA:                                             {source_attr['RA_mean']:.6f}")
    print(f"Mean DEC:                                            {source_attr['DEC_mean']:.6f}\n")
    if source_attr['Isinmask'] == False:
        print('Mean coordinates of the source fall outside the threshold\n')
    if source_attr['Ismaxpos'] == False:
        print('The position of the maximum pixel of the source does not correspond to the input coordinates\n')


def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("in_image",
                        help="Input image.")
    parser.add_argument("threshold",
                        help="Perform statistics on everything above threshold (Jy/beam).")
    parser.add_argument("--coord", nargs='+', default=None,
                        help="""Source coordinates (RA, DEC) of the source,
                                input as RA DEC (default=center of the image)""")
    parser.add_argument("-i", "--id", default=0,
                        help="Optional source id to keep track of the stats (default=0)")
    parser.add_argument("--plot", action='store_true',
                        help="Plot the source along with threshold")

    return parser

if __name__ == '__main__':
	main()