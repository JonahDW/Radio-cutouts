# CASA-Poststamp

This is a module for getting poststamps and statistics from CASA images, either manually or in an automated fashion. It makes use of CASA, matplotlib and some astropy routines. 

## catalog_cutouts.py

Script that wraps around the `casa_subimage.py` and `casa_stats.py` scripts in the module to create cutouts and measure statistics of sources from a catalog. If a 'Quality_flag' column is present in the catalog this is used to determine which sources are selected, sources with a quality flag of 0 will be measured. 

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

Script that can be ran after `catalog_cutouts.py` in order to visually inspect the cutouts and assign flags indicating the quality of the source and how well it is represented by the Gaussian catalog. Input is simply the directory containing cutouts files as well as the `source_cutout_stats.csv` file which should have been created by `catalog_cutouts.py`.

```
usage: python visual_classification.py cutout_dir
```

## casa_subimage.py

Make a subimage in CASA. Pass coordinates and size, in either pixel or sky coordinates. Standard output is another CASA image, but this can be converted to fits.

```
usage: casa_subimage.py [-h] [--pixel] [-x X] [-y Y] [--xsize XSIZE]
                        [--ysize YSIZE] [--fits]
                        in_image out_image

positional arguments:
  in_image       Input image, either CASA or fits format.
  out_image      Name of the output cutout image.

optional arguments:
  -h, --help     show this help message and exit
  --pixel        If this option is selected, input should be given in pixel
                 coordinates (default=sky coordinates).
  -x X           RA/X central coordinate (default=100).
  -y Y           DEC/Y central coordinate (default=100).
  --xsize XSIZE  Size of image in RA/X (default=100).
  --ysize YSIZE  Size of image in DEC/Y (default=100).
  --fits         Store the output image as a fits file.
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
