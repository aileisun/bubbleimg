# Bubbleimg

Bubbleimg is a package for handeling astronomy SDSS and Subaru HSC Survey Data. 

## What it does

The primary goal is to image emission line regions of type 2 AGN and it is designed to work on large samples of galaxies. It is also useful for downloading data (catalog, spectrum, and cutout images), measure spectral line strengths, analyze image, and do certain types of data visualization. 


## How to use it

See documentation:

    * [Start with one object](docs/obsobj.rst)

    * [How to download images](docs/imgdownload.rst)

    * [Apply it to the entire sample](docs/batch.rst)


## Downloading

Bubbleimg is still in beta version and `python setup.py install` may not work properly. The recommended method is to clone or download the source code and add the path to your .bashrc.. 

~~~~
export "PYTHONPATH=$PYTHONPATH:/where/you/put/the/code/" 
~~~~

Currently bubbleimg works under [lsst's anaconda environment](https://pipelines.lsst.io/install/conda.html) with Python 2.7.  

