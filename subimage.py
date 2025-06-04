import os
import sys
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

from argparse import ArgumentParser

class CutoutRegion():

    def __init__(self, pixel=False, x=None, y=None, xsize=None, ysize=None, region=None):

        self.pixel = pixel

        if region:
            self.region = region

            self.center = self.region.center
            self.ysize = self.region.width
            self.xsize = self.region.height
        else:
            self.xsize = xsize
            self.ysize = ysize

            if self.pixel:
                self.center = (x, y)
            else:
                self.center = SkyCoord(x, y, unit=u.deg)

    def do_astropy_subimage(self, input_image, output_image):
        '''
        Make an image cutout

        Keyword arguments:
        input_image (string) -- Filename of the image
        output_image (string) -- Destination filename of cutout
        x, y (float) -- Central coordinates (pixel or sky in degrees)
        xsize, ysize (float) -- Size of the cutout (pixel or sky in degrees)
        pixel (bool) -- Whether input coordinates are in pixel coordinates
        '''
        from astropy.io import fits
        from astropy.wcs import WCS
        from astropy.nddata.utils import Cutout2D

        if self.pixel:
            size = u.Quantity((self.xsize,self.ysize), u.pix)
        else:
            size = u.Quantity((self.xsize,self.ysize), u.deg)

        # Open image
        hdu = fits.open(input_image)[0]
        wcs = WCS(hdu.header)
        img = np.squeeze(hdu.data)

        # make cutout
        cutout = Cutout2D(img, self.center, size, wcs=wcs.celestial)

        # Update the FITS header with the cutout WCS, preserve axes
        new_wcs = cutout.wcs.celestial
        new_header = cutout.wcs.to_header()
        new_header["WCSAXES"] = hdu.header['NAXIS']
        new_hdu = fits.PrimaryHDU(data=cutout.data[np.newaxis,np.newaxis,:,:], 
                                  header=new_header)

        # Write the cutout to a new FITS file
        cutout_filename = output_image + '.fits'
        new_hdu.writeto(cutout_filename, overwrite=True)

    def do_casa_subimage(self, input_image, output_image, fits):
        '''
        Make an image cutout

        Keyword arguments:
        input_image (string) -- Filename of the image
        output_image (string) -- Destination filename of cutout
        x, y (float) -- Central coordinates (pixel or sky in degrees)
        xsize, ysize (float) -- Size of the cutout (pixel or sky in degrees)
        pixel (bool) -- Whether input coordinates are in pixel coordinates
        fits (bool) -- Whether to additionally out a fits file
        '''
        from casatasks import imsubimage, exportfits

        if self.region:
            reg = self.region.serialize(format='crtf')
        else:
            if self.pixel:
                reg = (f"centerbox[[{int(self.x)}pix,{int(self.y)}pix],"
                       f"[{int(self.xsize)}pix,{int(self.ysize)}pix]]")
            else:
                reg = (f"centerbox[[{self.x}deg,{self.y}deg],"
                       f"[{self.xsize}deg,{self.ysize}deg]]")

        imsubimage(imagename=input_image,
                   outfile=output_image+'.image',
                   region=reg,
                   overwrite=True)

        if fits:
            new_input = output_image+'.image'
            if not output_image.endswith('.fits'):
                output_image = output_image + '.fits'
            exportfits(new_input, output_image, dropdeg=True, overwrite=True)

def regions_from_file(region_file):
    from regions import Regions

    regions = Regions.read(region_file)
    return regions

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    input_image = args.in_image
    out_images = args.out_image
    region_file = args.regions
    casa = args.casa

    pixel = args.pixel
    x = args.x
    y = args.y
    xsize = args.xsize
    ysize = args.ysize

    fits = args.fits

    if region_file is not None:
        out_image_base = region_file.rsplit('.',1)[0]

        regions = regions_from_file(region_file)
        cutout_regions = [CutoutRegion(pixel, region=region) for region in regions]
    else:
        out_image_base = input_image.rsplit('.',1)[0]

        cutout_regions = []
        for i in range(len(x)):
            cutout_region = CutoutRegion(pixel, x=x, y=y, 
                                         xsize=xsize, ysize=ysize)
            cutout_regions.append(cutout_region)

    if out_images is None:
        out_images = [out_image_base+'_reg'+str(i+1) for i in range(len(cutout_regions))]

    if casa:
        for i, cr in enumerate(cutout_regions):
            cr.do_casa_subimage(input_image, out_images[i], fits)
    else:
        for i, cr in enumerate(cutout_regions):
            cr.do_astropy_subimage(input_image, out_images[i])

def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("in_image",
                        help="Input image, either CASA or fits format.")
    parser.add_argument("--out_image", nargs='*', default=None,
                        help="Name(s) of the output cutout image(s).")
    parser.add_argument("--regions", default=None, type=str,
                        help="""Region file to read from. Ignores all other 
                              position and size options""")
    parser.add_argument("--casa", action='store_true',
                        help="Whether to use casa functionality (otherwise astropy)")
    parser.add_argument("--pixel", action='store_true',
                        help="""If this option is selected, input should be 
                                given in pixel coordinates (default=sky coordinates).""")
    parser.add_argument("-x", nargs='+', type=float,
                        help="RA/X central coordinate(s).")
    parser.add_argument("-y", nargs='+', type=float,
                        help="DEC/Y central coordinate(s).")
    parser.add_argument("--xsize", nargs='+', type=float,
                        help="Size(s) of image(s) in RA/X.")
    parser.add_argument("--ysize", nargs='+', type=float,
                        help="Size(s) of image(s) in DEC/Y.")
    parser.add_argument("--fits", action='store_true',
                        help="Store the output image(s) as a fits file.")
    return parser

if __name__ == '__main__':
	main()