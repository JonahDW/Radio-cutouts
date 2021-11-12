import os
import sys
import numpy as np

from astropy.table import Table, join
from astropy.io import fits

from argparse import ArgumentParser
from pathlib import Path

def do_subimage(source, image_file, cutout_file):
    '''
    Call the subimage script
    '''
    xsize = 5*source['Maj']
    ysize = 5*source['Maj']

    casa_subimage_script = Path(__file__).parent / 'casa_subimage.py'
    os.system(f"casa --nologfile --nologger -c {casa_subimage_script} {image_file} {cutout_file}"
              +f" -x {source['RA']} -y {source['DEC']} --xsize {xsize} --ysize {ysize} --fits")

def do_stats(source, threshold, cutout_file):
    '''
    Call the stats script
    '''
    rms = source['Isl_rms']

    casa_stats_script = Path(__file__).parent / 'casa_stats.py'
    os.system(f"casa --nologfile --nologger -c {casa_stats_script} {cutout_file}.image {threshold*rms}"
              +f" -i {source['Source_id']} --coord {source['RA_max']} {source['DEC_max']} --plot")
    os.system(f'rm -r {cutout_file}.image')

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    image_file = args.image
    catalog_file = args.catalog
    stats = args.stats
    threshold = args.threshold

    catalog = Table.read(catalog_file)
    name = catalog.meta['OBJECT'].replace("'","")

    out_folder = name+'_cutouts'
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    if os.path.exists('source_stats.csv'):
        os.system('rm source_stats.csv')

    if 'Quality_flag' in catalog.columns:
        selected_sources = catalog[catalog['Quality_flag'] == 0]
    else:
        print('Quality flag column not found, selecting all sources in the catalog')
        selected_sources = catalog

    if not stats:
        print(f'Making cutouts for {len(selected_sources)} sources')
        for source in selected_sources:
            cutout_file = os.path.join(out_folder,source['Source_name'].replace(' ','_'))
            do_subimage(source, image_file, cutout_file)

        os.system('rm *.last')

    if stats:
        print(f'Making cutouts and measuring fluxes for {len(selected_sources)} sources')
        for source in selected_sources:
            cutout_file = os.path.join(out_folder,source['Source_name'].replace(' ','_'))
            do_subimage(source, image_file, cutout_file)
            do_stats(source, cutout_file)

        os.system('rm *.last')

        # Make sure these columns are unused
        if 'Cutout_flag' in catalog.columns:
            del catalog['Cutout_flag']
            del catalog['RA_Mean']
            del catalog['DEC_Mean']
            del catalog['Isl_Int_flux']

        # Read in created CSV file
        source_stats = Table.read('source_stats.csv')
        catalog = join(catalog, source_stats, join_type='left', keys='Source_id')

        catalog[np.where(catalog['is_inmask'])]['Cutout_flag'] = 1
        catalog[np.where(catalog['is_maxpos'])]['Cutout_flag'] = 2
        catalog[abs(catalog['Isl_Total_flux']/catalog['Isl_Int_flux'] - 1) > 0.2]['Cutout_flag'] = 3
        catalog.write(catalog_file, overwrite=True)


def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("image",
                        help="""Input image.""")
    parser.add_argument("catalog",
                        help="""Input catalog""")
    parser.add_argument("-s", "--stats", action='store_true',
                        help="""Measure the statistics for all 
                                the sources and write to the catalog""")
    parser.add_argument("-t", "--threshold", default=3.0,
                        help="Threshold for measuring statistics to use (in sigma, default=3.0)")

    return parser

if __name__ == '__main__':
    main()