import os
from argparse import ArgumentParser

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    input_image = args.in_image

    output_image = args.out_image
    output_image = output_image.split('.')[0]

    pixel = args.pixel
    x = args.x
    y = args.y

    xsize = args.xsize
    ysize = args.ysize
    fits = args.fits

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
        os.system('rm *.last')

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