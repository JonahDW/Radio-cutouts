import os
import sys
import aplpy

import matplotlib.axes as maxes
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import PercentileInterval, ImageNormalize, AsinhStretch

from argparse import ArgumentParser

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    image_file = args.image_file
    percentile = args.percentile
    cmap = args.cmap
    dpi = args.dpi

    im_id = image_file.split('.')[0]

    hdu = fits.open(image_file)[0]
    wcs = WCS(hdu.header)

    norm = ImageNormalize(hdu.data*1e3, interval=PercentileInterval(percentile), stretch=AsinhStretch())

    ax = plt.subplot(projection=wcs)
    im = ax.imshow(hdu.data*1e3, origin='lower', cmap=cmap, norm=norm)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label('Flux (mJy/beam)')

    #ax.grid(color='white', ls='dotted')
    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('DEC (J2000)')

    plt.savefig(im_id+'.png', dpi=dpi, bbox_inches='tight')

def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("image_file",
                        help="""Input image in fits format.""")
    parser.add_argument("--percentile", default=99.9,
                        help="Percentile interval to use for normalization (default=99.9).")
    parser.add_argument("--cmap", default='gist_heat',
                        help="What colormap to use (default=gist_heat)")
    parser.add_argument("--dpi", default=300,
                        help="Dpi of the image (default=300).")
    return parser

if __name__ == '__main__':
    main()