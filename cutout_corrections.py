import os
import sys
import numpy as np

from astropy.table import Table, join

from argparse import ArgumentParser
from pathlib import Path


def flux_correction(catalog, flux_col, reset_flux):

    I_flag = catalog['Cutout_class'] == 'I'
    P_flag = catalog['Cutout_class'] == 'P'

    replace_col = flux_col+'_measured'
    if reset_flux:
        if replace_col in catalog.colnames:
            catalog[flux_col] = catalog[replace_col].copy()
        else:
            flux_col_index = catalog.colnames.index(flux_col)
            orig_flux_col = catalog[flux_col].copy()
            catalog.rename_column(flux_col, replace_col)

            catalog.add_column(orig_flux_col,
                               name=flux_col,
                               index=flux_col_index)

    catalog[flux_col][I_flag] = catalog['Isl_Total_flux'][I_flag]
    catalog[flux_col][P_flag] = catalog['Peak_flux'][P_flag]

    return catalog

def alpha_correction(catalog, alpha_col, reset_alpha):

    I_flag = catalog['Cutout_class'] == 'I'

    replace_col = alpha_col+'_measured'
    if reset_alpha:
        if replace_col in catalog.colnames:
            catalog[alpha_col] = catalog[replace_col].copy()
        else:
            alpha_col_index = catalog.colnames.index(alpha_col)
            orig_alpha_col = catalog[alpha_col].copy()
            catalog.rename_column(alpha_col, replace_col)

            catalog.add_column(orig_alpha_col,
                               name=alpha_col,
                               index=alpha_col_index)

    catalog[alpha_col][I_flag] = catalog['Cutout_Spectral_index'][I_flag]

    return catalog

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    catalog_file = args.catalog
    flux_col = args.flux_col
    reset_flux = args.reset_flux
    alpha_col = args.alpha_col
    reset_alpha = args.reset_alpha

    catalog = Table.read(catalog_file)

    if flux_col:
        catalog = flux_correction(catalog, flux_col, reset_flux)
    if alpha_col:
        catalog = alpha_correction(catalog, alpha_col, reset_alpha)

    # Fill N_Gaus column while we're at it
    catalog['N_Gaus'][catalog['S_Code'] == 'S'] = 1

    catalog.write(catalog_file, overwrite=True)


def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("catalog",
                        help="Input catalog to correct spectral index and flux")
    parser.add_argument('-f', '--flux_col', default=None,
                        help="""Flux column for to apply flux correction on
                                (default = none, don't apply flux corrections)""")
    parser.add_argument('-a', '--alpha_col', default=None,
                        help="""Spectral index column for to apply flux correction on
                                (default = none, don't apply flux corrections)""")
    parser.add_argument('--reset_flux', action='store_true',
                        help="""Reset flux column to their original values before
                                applying corrections.""")
    parser.add_argument('--reset_alpha', action='store_true',
                        help="""Reset alpha column to their original values before
                                applying corrections.""")

    return parser


if __name__ == '__main__':
	main()