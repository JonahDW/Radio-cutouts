import os
import json
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from argparse import ArgumentParser

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

    # If starting position not in mask, find nearest pixel that is
    if mask[start_x, start_y] == 0:
        start_x, start_y = find_nearest_nonzero(mask, start_x, start_y)

    floodfill(mask, start_x, start_y)
    mask[mask == 1] = 0
    mask[mask == 2] = 1

    return mask

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    input_image = args.in_image
    threshold = args.threshold
    source_id = args.id
    coord = args.coord
    plot = args.plot

    # Create mask using the threshold
    im.mask(input_image, mask='stats_mask', threshold=threshold)

    # Get the maximum position
    max_ra = np.deg2rad(float(coord[0]))
    max_dec = np.deg2rad(float(coord[1]))

    ia.open(input_image)
    pixel_coord = ia.topixel([max_ra, max_dec])
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

    is_inmask = mask_2d[int(x_cent),int(y_cent)]
    is_maxpos = ((start_x - image_stats['maxpos'][0])**2
                 + (start_y - image_stats['maxpos'][1])**2)**0.5 < 2.

    if plot:
        fig = plt.figure()
        fig.add_subplot(111)

        plt.imshow(image.T, interpolation='none', 
                   norm=TwoSlopeNorm(vcenter=float(threshold)),
                   cmap='jet')
        plt.contour(mask_2d.T, colors='k', levels=[0,1], linewidths=1)
        plt.scatter(x_cent, y_cent, color='w', marker='x', s=5, zorder=2)

        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        plt.savefig(input_image.split('.')[0]+'.png')

    if os.path.isfile('source_stats.csv'):
        with open('source_stats.csv','a') as f:
            f.write(f"{source_id},{ra_deg},{dec_deg},{image_stats['flux'][0]},{is_inmask},{is_maxpos}")
            f.write('\n')
    else:
        with open('source_stats.csv','w') as f:
            f.write('Source_id,RA_Mean,DEC_Mean,Isl_Int_flux,Isinmask,Ismaxpos')
            f.write('\n')
            f.write(f"{source_id},{ra_deg},{dec_deg},{image_stats['flux'][0]},{is_inmask},{is_maxpos}")
            f.write('\n')

    os.system('rm -r stats_mask')

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