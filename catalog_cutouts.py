import os
import sys
import numpy as np

from astropy.table import Table, join
from astropy.io import fits

from argparse import ArgumentParser
from pathlib import Path

from casa_subimage import do_subimage
from casa_stats import Image

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    image_file = args.image
    catalog_file = args.catalog
    alpha_file = args.alpha_image
    gaul_file = args.gaul
    stats = args.stats
    threshold = args.threshold

    catalog = Table.read(catalog_file)
    name = catalog.meta['OBJECT'].replace("'","")

    # Create cutout folder
    out_folder = os.path.join(os.path.dirname(catalog_file),name+'_cutouts')
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    # Check if quality flag column is present
    if 'Quality_flag' in catalog.columns:
        selected_sources = catalog[catalog['Quality_flag'] == 0]
    else:
        print('Quality flag column not found, selecting all sources in the catalog')
        selected_sources = catalog

    if gaul_file:
        gaul = Table.read(gaul_file)
        print(f'Reading in gaussian list from {gaul_file}')
    else:
        gaul = None

    if not stats:
        print(f'Making cutouts for {len(selected_sources)} sources')
        for source in selected_sources:

            # Get correct filenames for sources
            source_filename = source['Source_name'].replace(' ','_')
            cutout_file = os.path.join(out_folder,source_filename)

            do_subimage(image_file, cutout_file, source['RA'], source['DEC'], 
                        xsize=5*source['Maj'], ysize=5*source['Maj'], pixel=False, fits=True)

            os.system(f'rm -r {cutout_file}.image')

            if alpha_file:
                cutout_file = os.path.join(out_folder,source_filename+'_alpha')

                do_subimage(alpha_file, cutout_file, source['RA'], source['DEC'],
                            xsize=5*source['Maj'], ysize=5*source['Maj'], pixel=False, fits=True)

                os.system(f'rm -r {cutout_file}.image')


    if stats:
        print(f'Making cutouts and measuring fluxes for {len(selected_sources)} sources')\
        # Create table with for holding the stats

        all_stats = []
        for source in selected_sources:

            # Get correct filenames for sources
            source_filename = source['Source_name'].replace(' ','_')
            cutout_file = os.path.join(out_folder,source_filename)

            do_subimage(image_file, cutout_file, x=source['RA'], y=source['DEC'], 
                        xsize=5*source['Maj'], ysize=5*source['Maj'], pixel=False, fits=True)

            if alpha_file:
                alpha_cutout_file = os.path.join(out_folder,source_filename+'_alpha')

                do_subimage(alpha_file, alpha_cutout_file, source['RA'], source['DEC'],
                            xsize=5*source['Maj'], ysize=5*source['Maj'], pixel=False, fits=True)

                alpha_cutout_file = alpha_cutout_file+'.image'
            else:
                alpha_cutout_file = None

            rms = source['Isl_rms']
            max_coord = (source['RA_max'], source['DEC_max'])
            gaussians = None
            if gaul:
                gaussians = gaul[gaul['Source_id'] == source['Source_id']]

            image = Image(cutout_file+'.image', threshold*rms, start_coordinates=max_coord,
                          source_id=int(source['Source_id']))
            source_stats = image.get_stats(plot=True, source=source, gaussians=gaussians,
                                           alpha_image=alpha_cutout_file)

            cutout_flag = ''
            if source_stats['Isinmask'] == False:
                cutout_flag +='M'
            if source_stats['Ismaxpos'] == False:
                cutout_flag += 'C'
            if abs(source['Isl_Total_flux']/source_stats['Cutout_Total_flux'] - 1) > 0.2:
                cutout_flag += 'F'

            source_stats['Source_name'] = source['Source_name']
            # Change dictionary to include new flag
            source_stats['Cutout_flag'] = cutout_flag
            del source_stats['Isinmask']
            del source_stats['Ismaxpos']
            all_stats.append(source_stats)

            os.system(f'rm -r {cutout_file}.image')
            if alpha_file:
                os.system(f'rm -r {alpha_cutout_file}')

        source_table = Table(rows=all_stats)

        outfile = os.path.join(out_folder,'source_cutout_stats.csv')
        # Preserve classification if it exists
        if os.path.exists(outfile):
            old_table = Table.read(outfile)
            if 'Cutout_class' in old_table.colnames and len(old_table) == len(source_table):
                print("Preserving classification from old table")
                source_table['Cutout_class'] = old_table['Cutout_class']

        # Write to file
        source_table.write(outfile, overwrite=True)

def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("image", type=str,
                        help="Input image.")
    parser.add_argument("catalog", type=str,
                        help="Input source catalog")
    parser.add_argument("-a", "--alpha_image", default=None,
                        help="Spectral index image")
    parser.add_argument("-g", "--gaul", default=None,
                        help="Input gaussian list catalog, only used in plotting sources")
    parser.add_argument("-s", "--stats", action='store_true',
                        help="""Measure the statistics for all 
                                the sources and write to the catalog""")
    parser.add_argument("-t", "--threshold", default=3.0, type=float,
                        help="Threshold for measuring statistics to use (in sigma, default=3.0)")

    return parser

if __name__ == '__main__':
    main()