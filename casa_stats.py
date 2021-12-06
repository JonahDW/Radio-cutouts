import os
import json
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
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
            if x < len(matrix[y]) - 1:
                floodfill(matrix,x+1,y)
            if y > 0:
                floodfill(matrix,x,y-1)
            if y < len(matrix) - 1:
                floodfill(matrix,x,y+1)
            if x > 0 and y > 0:
                floodfill(matrix,x-1,y-1)
            if x > 0 and y < len(matrix) - 1:
                floodfill(matrix,x-1,y+1)
            if x < len(matrix) -1 and y > 0:
                floodfill(matrix,x+1,y-1)
            if x < len(matrix) - 1 and y < len(matrix) - 1:
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
        id_col = 'Gaus_id'
    else:
        id_col = 'Source_id'

    # Get ids
    ids = gaussians[id_col]
    increment = max(wcs.increment(format='n', type='direction')['numeric'])

    ell = []
    for source in gaussians:
        xy = wcs.topixel([np.deg2rad(source['RA']), np.deg2rad(source['DEC'])])
        ell.append(Ellipse(xy = xy['numeric'][:2],
                           width = source['Min']/np.rad2deg(increment),
                           height = source['Maj']/np.rad2deg(increment),
                           angle = -source['PA']))

    return ell, ids

def get_stats(input_image, threshold, coord, source_id=0, plot=False, gaussians=None):
    '''
    Get the stats of the source in an image

    Keyword arguments:
    input_image (string) -- Filename of the input image
    threshold (float) -- Threshold for stats (in general in units of Jy/beam)
    coord (tuple) -- Starting coordinates (RA, DEC) to identify correct island
    source_id (string) -- Source id to associate source with
    plot (bool) -- Output plot with the source, threshold, and optionally Gaussians
    gaussians (table object) -- List of Gaussians associated with the source
    '''
    im = ct.imager()
    ia = ct.image()

    # Create mask using the threshold
    im.mask(input_image, mask='stats_mask', threshold=threshold)

    # Get the maximum position
    max_ra = np.deg2rad(float(coord[0]))
    max_dec = np.deg2rad(float(coord[1]))

    ia.open(input_image)
    wcs = ia.coordsys([0,1])
    pixel_coord = wcs.topixel([max_ra, max_dec])
    start_x = int(pixel_coord['numeric'][0])
    start_y = int(pixel_coord['numeric'][1])

    ia.close()

    # Open mask and retain only pixels connected to the source
    ia.open('stats_mask')
    mask = ia.getregion()

    mask_2d = mask[:,:,0,0]
    mask_2d = prune_mask(mask_2d, start_x, start_y)

    mask[:,:,0,0] = mask_2d
    ia.putregion(mask)
    ia.close()

    # Run imstat to get integrated flux
    image_stats = imstat(input_image, mask='stats_mask')

    # Open image and get center of mass
    ia.open(input_image)
    image = ia.getregion()[:,:,0,0]

    valid_pix = np.where(mask_2d)
    x_cent = np.average(valid_pix[0], weights=image[valid_pix])
    y_cent = np.average(valid_pix[1], weights=image[valid_pix])

    center_of_mass = ia.toworld([x_cent,y_cent], 'n')
    ra_deg = np.rad2deg(center_of_mass['numeric'][0])
    dec_deg = np.rad2deg(center_of_mass['numeric'][1])

    is_inmask = bool(mask_2d[int(x_cent),int(y_cent)])
    is_maxpos = ((start_x - image_stats['maxpos'][0])**2
                 + (start_y - image_stats['maxpos'][1])**2)**0.5 < 2.

    os.system('rm -r stats_mask')

    if plot:
        fig = plt.figure()
        fig.add_subplot(111)

        plt.imshow(image.T, interpolation='none', origin='lower',
                   norm=TwoSlopeNorm(vcenter=float(threshold)),
                   cmap='jet')
        plt.contour(mask_2d.T, colors='k', levels=[0,1], linewidths=1)
        plt.scatter(x_cent, y_cent, color='w', marker='x', s=5, zorder=2)

        if gaussians:
            ell, ids = include_ellipses(gaussians, wcs)
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

        plt.xticks([])
        plt.yticks([])
        #plt.tight_layout()
        plt.savefig(input_image.split('.')[0]+'.png')
        plt.close()

    source_attr = {'Source_id':source_id,
                   'RA_mean': ra_deg,
                   'DEC_mean':dec_deg,
                   'Cutout_Total_flux':image_stats['flux'][0],
                   'Isinmask':is_inmask,
                   'Ismaxpos':is_maxpos}
    return source_attr

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    input_image = args.in_image
    threshold = args.threshold
    source_id = args.id
    coord = args.coord
    plot = args.plot

    source_attr = get_stats(input_image, threshold, coord, source_id, plot)

    print(f'Measured properties for source with source id {source_id} in {input_image}: \n')
    print(f"Total flux above threshold ({threshold:.2f}):    {source_attr['Cutout_Total_flux']}")
    print(f"Mean RA:                                         {source_attr['RA_mean']:.6f}")
    print(f"Mean DEC:                                        {source_attr['DEC_mean']:.6f}\n")
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