import sys
import os
from argparse import ArgumentParser

from casatasks import imsubimage, exportfits

def do_subimage(input_image, output_image, x, y, xsize, ysize, pixel, fits):
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
    if pixel:
        reg = f'box[[{int(x-xsize/2)}pix,{int(y-ysize/2)}pix],[{int(x+xsize/2)}pix,{int(y+ysize/2)}pix]]'
    else:
        reg = f'box[[{x-xsize/2}deg,{y-ysize/2}deg],[{x+xsize/2}deg,{y+ysize/2}deg]]'

    imsubimage(imagename=input_image,
               outfile=output_image+'.image',
               region=reg,
               overwrite=True)

    if fits:
        new_input = output_image+'.image'
        if not output_image.endswith('.fits'):
            output_image = output_image + '.fits'
        exportfits(new_input, output_image, dropdeg=True, overwrite=True)

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    input_image = args.in_image

    output_image = args.out_image
    output_image = output_image.rsplit('.',1)[0]

    pixel = args.pixel
    x = args.x
    y = args.y

    xsize = args.xsize
    ysize = args.ysize
    fits = args.fits

    do_subimage(input_image, output_image, 
                x, y, xsize, ysize, 
                pixel, fits)

def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("in_image",
                        help="Input image, either CASA or fits format.")
    parser.add_argument("out_image",
                        help="Name of the output cutout image.")
    parser.add_argument("--pixel", action='store_true',
                        help="""If this option is selected, input should be 
                                given in pixel coordinates (default=sky coordinates).""")
    parser.add_argument("-x",  default=100, type=float,
                        help="RA/X central coordinate (default=100).")
    parser.add_argument("-y", default=100, type=float,
                        help="DEC/Y central coordinate (default=100).")
    parser.add_argument("--xsize", default=100, type=float,
                        help="Size of image in RA/X (default=100).")
    parser.add_argument("--ysize", default=100, type=float,
                        help="Size of image in DEC/Y (default=100).")
    parser.add_argument("--fits", action='store_true',
                        help="Store the output image as a fits file.")
    return parser

if __name__ == '__main__':
	main()