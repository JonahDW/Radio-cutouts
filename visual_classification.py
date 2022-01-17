import os
import sys

import numpy as np

from astropy.table import Table, Column

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, CheckButtons

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
    cutout_folder = sys.argv[1]

    source_stats = Table.read(os.path.join(cutout_folder, 'source_cutout_stats.csv'))
    classification = Classification()

    for source in source_stats:
        source_file = os.path.join(cutout_folder, source['Source_name'].replace(' ','_')+'.png')
        make_cutout_button(source_file, classification)

    flags = classification.get_results()
    source_stats.add_column(flags, name='Cutout_class')

    source_stats.write(os.path.join(cutout_folder, 'source_cutout_stats.csv'), overwrite=True)


if __name__ == '__main__':
    main()