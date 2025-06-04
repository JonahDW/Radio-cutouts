# Radio-cutouts

This is a module for getting poststamps and statistics from CASA images, either manually or in an automated fashion. It makes use of CASA, matplotlib and some astropy routines. 

## subimage.py

Simple script to make cutouts. Make one or multiple subimages from an image, using either CASA or Astropy. Coordinates and size can be given in either pixel or sky coordinates, or a (DS9) region file can be specified.

```
usage: subimage.py [-h] [--out_image [OUT_IMAGE [OUT_IMAGE ...]]]
                   [--regions REGIONS] [--casa] [--pixel] [-x X [X ...]]
                   [-y Y [Y ...]] [--xsize XSIZE [XSIZE ...]]
                   [--ysize YSIZE [YSIZE ...]] [--fits]
                   in_image

positional arguments:
  in_image              Input image, either CASA or fits format.

optional arguments:
  -h, --help            show this help message and exit
  --out_image [OUT_IMAGE [OUT_IMAGE ...]]
                        Name(s) of the output cutout image(s).
  --regions REGIONS     Region file to read from. Ignores all other position
                        and size options
  --casa                Whether to use casa functionality (otherwise astropy)
  --pixel               If this option is selected, input should be given in
                        pixel coordinates (default=sky coordinates).
  -x X [X ...]          RA/X central coordinate(s).
  -y Y [Y ...]          DEC/Y central coordinate(s).
  --xsize XSIZE [XSIZE ...]
                        Size(s) of image(s) in RA/X.
  --ysize YSIZE [YSIZE ...]
                        Size(s) of image(s) in DEC/Y.
  --fits                Store the output image(s) as a fits file.
```

## catalog_cutouts.py

Script that wraps around the `subimage.py` and `stats.py` scripts in the module to create cutouts and measure statistics of sources from a catalog. If a 'Quality_flag' column, which is created and can be manipulated with the [Image-processing]{https://github.com/JonahDW/Image-processing} module, is present in the catalog, this is used to determine which sources are selected, sources with a quality flag of 0 will be measured. 

```
usage: catalog_cutouts.py [-h] [-a ALPHA_IMAGE] [-g GAUL] [-s] [-t THRESHOLD]
                          image catalog

positional arguments:
  image                 Input image.
  catalog               Input source catalog

optional arguments:
  -h, --help            show this help message and exit
  -a ALPHA_IMAGE, --alpha_image ALPHA_IMAGE
                        Spectral index image
  -g GAUL, --gaul GAUL  Input gaussian list catalog, only used in plotting
                        sources
  -s, --stats           Measure the statistics for all the sources and write
                        to the catalog
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for measuring statistics to use (in sigma,
                        default=3.0)
```

## visual_classification.py

Script that can be ran after `catalog_cutouts.py` in order to visually inspect the cutouts and assign flags indicating the quality of the sources. Input is simply the directory containing cutouts files as well as the `source_cutout_stats.csv` file which should have been created by `catalog_cutouts.py`.

```
usage: visual_classification.py [-h] [--append_results [APPEND_RESULTS]]
                                cutout_dir result_name

positional arguments:
  cutout_dir            Directory containing source cutouts
  result_name           Unique name to write results to i.e.
                        jonah_classification. Results will be written to a csv
                        file with that name, if it already exists it will be
                        overwritten.

optional arguments:
  -h, --help            show this help message and exit
  --append_results [APPEND_RESULTS]
                        Append results to the source cutout stats csv file. If
                        another file than the result_name file is specified
                        (which exists in cutout_dir) this file is appended
                        instead and visual inspection is skipped.
```

## casa_stats.py

Measure the statistics of a CASA (cutout) image and store them in a csv file. Statistics are measured at given RA and DEC above given threshold, and thus is only done for single islands of emission. Output stored in the `source_stats.csv` consists at the moment of: the intensity weighted mean RA and DEC, the integrated flux of the island, and columns indicating whether the mean RA/DEC is located within the bounds set by the threshold, and whether the maximum RA/DEC match with the input RA/DEC. A source is assigned a source ID if that is specified.

```
usage: casa_stats.py [-h] [--coord COORD [COORD ...]] [-i ID] [--plot]
                     in_image threshold

positional arguments:
  in_image              Input image.
  threshold             Perform statistics on everything above threshold
                        (Jy/beam).

optional arguments:
  -h, --help            show this help message and exit
  --coord COORD [COORD ...]
                        Source coordinates (RA, DEC) of the source, input as
                        RA DEC (default=center of the image)
  -i ID, --id ID        Optional source id to keep track of the stats
                        (default=0)
  --plot                Plot the source along with threshold
```
