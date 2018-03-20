# Bubbleimg

Bubbleimg is a package for handeling Sloan Digital Sky Survey (SDSS) and Subaru Hyper Suprime-Cam (HSC) Survey Data. 

## What it does

The primary goal is to image emission line regions of type 2 AGN and it is designed to work on large samples of galaxies. It is also useful for downloading data (catalog, spectrum, and cutout images), measuring spectral line strengths, analysing image, and for certain types of data visualization. 


## How to use it

See documentation:

* [Start with one object](docs/obsobj.rst)

* [How to download images](docs/imgdownload.rst)

* [Apply it to the entire sample](docs/batch.rst)


## Downloading

The recommended way of running bubbleimg is to clone or download the source code and add the path to your `~/.bashrc`. 

~~~~
export "PYTHONPATH=$PYTHONPATH:/where/you/put/the/code/" 
~~~~

There are a number of required packages, as listed in `setup.py`, including a few that needs to be installed manually. 

Bubbleimg runs under Python 3

