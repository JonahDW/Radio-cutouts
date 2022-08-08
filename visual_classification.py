import os
import sys

import numpy as np
from astropy.table import Table, join

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, CheckButtons

from argparse import ArgumentParser

class Classification():

    def __init__(self):
        self.checkbuttons = None

        self.options = ''
        self.results = []

    def add_new_check(self,checkbuttons):
        self.checkbuttons = checkbuttons
        self.options = ''

    def go_next(self,event):
        status = self.checkbuttons.get_status()

        if status[0]:
            self.options += 'I'
        if status[1]:
            self.options += 'G'
        if status[2]:
            self.options += 'P'
        if status[3]:
            self.options += 'A'

        self.results.append(self.options)
        plt.close()

    def get_results(self):
        return self.results

def crop_plot(image):
    ''' Crop image such that whitespace is gone'''
    center_x = image.shape[0]//2
    center_y = image.shape[1]//2

    top = np.where((image[:,center_y,0] == 0) & (image[:,center_y,1] == 0) & (image[:,center_y,2] == 0))[0][0]
    bottom = np.where((image[:,center_y,0] == 0) & (image[:,center_y,1] == 0) & (image[:,center_y,2] == 0))[0][-1]

    left = np.where((image[center_x,:,0] == 0) & (image[center_x,:,1] == 0) & (image[center_x,:,2] == 0))[0][0]
    right = np.where((image[center_x,:,0] == 0) & (image[center_x,:,1] == 0) & (image[center_x,:,2] == 0))[0][-1]

    image = image[top:bottom+1,left-1:right+1,:]
    return image

def make_cutout_button(plotfile, classification):
    '''
    Make a collage of cutouts of a given source

    source: Name of the source
    file_id: File where the image of the source is stored
    ps1_table, wise_table: Input tables with all relevant data
    j: Position of the source in the input table

    Plots a collage of the source in all relevant filters
    '''
    fig, ax = plt.subplots(figsize=(10,5))

    image = plt.imread(plotfile, format='png')
    image = crop_plot(image)

    ax.imshow(image)
    ax.set_axis_off()

    plt.subplots_adjust(right=0.5)

    activated = [False, False, False, False]
    labels = ['I: Complex source, use island flux\n'
              '   and gravitational mean',
              'G: Source well described by Gaussians,\n'
              '   use as is',
              'P: Single valid Gaussian,\n'
              '   use peak flux instead',
              'A: Artifact']

    # Add check buttons
    ax_check = plt.axes([0.5, 0.3, 0.4, 0.5])
    classify_button = CheckButtons(ax_check, labels, activated)
    classification.add_new_check(classify_button)

    # Add next button
    ax_next = plt.axes([0.5, 0.15, 0.1, 0.1])
    next_button = Button(ax_next,'Next')
    next_button.on_clicked(classification.go_next)

    plt.show()

def main():

    parser = new_argument_parser()
    args = parser.parse_args()

    cutout_folder = args.cutout_dir
    result_name = args.result_name
    append_results = args.append_results

    source_stats_file = os.path.join(cutout_folder, 'source_cutout_stats.csv')
    source_stats = Table.read(source_stats_file)

    # If append_results is a filename, read in and write results to source cutout stats
    if isinstance(append_results, str):
        results_file = os.path.join(cutout_folder, append_results)
        if os.path.exists(results_file):
            results = Table.read(results_file)
            if 'Cutout_class' in source_stats.colnames:
                source_stats.remove_column('Cutout_class')
            source_stats = join(source_stats, results, keys='Source_name')
            source_stats.write(source_stats_file, overwrite=True)
        else:
            print(f'Specified results file {append_results} not found')
        sys.exit()

    # Do visual classification of sources
    classification = Classification()
    source_names = []
    for source in source_stats:
        source_names.append(source['Source_name'])
        source_file = os.path.join(cutout_folder, source['Source_name'].replace(' ','_')+'.png')
        make_cutout_button(source_file, classification)

    # Put results in table and write to result file
    results = Table()
    flags = classification.get_results()
    results['Source_name'] = source_names
    results['Cutout_class'] = flags
    if not result_name.endswith('.csv'):
        result_name = result_name+'.csv'

    result_file = os.path.join(cutout_folder, result_name)
    results.write(result_file, overwrite=True)

    # Put results in the source stats file
    if append_result:
        if 'Cutout_class' in source_stats.colnames:
            source_stats.remove_column('Cutout_class')
        source_stats = join(source_stats, results, keys='Source_name')
        source_stats.write(source_stats_file, overwrite=True)

def new_argument_parser():

    parser = ArgumentParser()

    parser.add_argument("cutout_dir", type=str,
                        help="Directory containing source cutouts")
    parser.add_argument("result_name", type=str,
                        help="""Unique name to write results to i.e. jonah_classification.
                                Results will be written to a csv file with that name, if it 
                                already exists it will be overwritten.""")
    parser.add_argument("--append_results", nargs="?", const=True,
                        help="""Append results to the source cutout stats csv file. If another
                                file than the result_name file is specified (which exists in cutout_dir) 
                                this file is appended instead and visual inspection is skipped.""")

    return parser

if __name__ == '__main__':
    main()