import os
import sys
import numpy as np

from astropy.table import Table, join
from astropy.io import fits

from argparse import ArgumentParser
from pathlib import Path

from casa_subimage import do_subimage
from casa_stats import get_stats

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    image_file = args.image
    catalog_file = args.catalog
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
            cutout_file = os.path.join(out_folder,source['Source_name'].replace(' ','_'))

            do_subimage(image_file, cutout_file, source['RA'], source['DEC'], 
                        xsize=5*source['Maj'], ysize=5*source['Maj'], pixel=False, fits=True)

            os.system(f'rm -r {cutout_file}.image')

        os.system('rm *.last')
        os.system('rm casa-*.log')

    if stats:
        print(f'Making cutouts and measuring fluxes for {len(selected_sources)} sources')\
        # Create table with for holding the stats
        source_table = Table(names=['Source_id','RA_mean','DEC_mean',
                                    'Cutout_Total_flux','Cutout_flag'],
                             dtype=[float,float,float,
                                    float, 'S3'])

        for source in selected_sources:
            cutout_file = os.path.join(out_folder,source['Source_name'].replace(' ','_'))

            do_subimage(image_file, cutout_file, x=source['RA'], y=source['DEC'], 
                        xsize=5*source['Maj'], ysize=5*source['Maj'], pixel=False, fits=True)

            rms = source['Isl_rms']
            max_coord = (source['RA_max'], source['DEC_max'])
            gaussians = None
            if gaul:
                gaussians = gaul[gaul['Source_id'] == source['Source_id']]

            source_stats = get_stats(cutout_file+'.image', threshold*rms, max_coord, 
                                     source_id=source['Source_id'], plot=True,
                                     gaussians=gaussians)

            cutout_flag = ''
            if source_stats['Isinmask'] == False:
                cutout_flag +='M'
            if source_stats['Ismaxpos'] == False:
                cutout_flag += 'C'
            if abs(source['Isl_Total_flux']/source_stats['Cutout_Total_flux'] - 1) > 0.2:
                cutout_flag += 'F'

            source_table.add_row([source_stats['Source_id'], source_stats['RA_mean'], 
                                  source_stats['DEC_mean'], source_stats['Cutout_Total_flux'],
                                  cutout_flag])

            os.system(f'rm -r {cutout_file}.image')

        # Match catalogs and write to file
        catalog = join(catalog, source_table, join_type='left', keys='Source_id')
        catalog = catalog.filled(np.nan)
        catalog.write(catalog_file, overwrite=True)

def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("image",
                        help="Input image.")
    parser.add_argument("catalog",
                        help="Input source catalog")
    parser.add_argument("-g", "--gaul", default=None,
                        help="Input gaussian list catalog, only used in plotting sources")
    parser.add_argument("-s", "--stats", action='store_true',
                        help="""Measure the statistics for all 
                                the sources and write to the catalog""")
    parser.add_argument("-t", "--threshold", default=3.0,
                        help="Threshold for measuring statistics to use (in sigma, default=3.0)")

    return parser

if __name__ == '__main__':
    main()