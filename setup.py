from setuptools import setup

"""
For full functionality, please download and install the following software: 

modelBC03
  https://github.com/aileisun/modelBC03
within ~/.bashrc please add:
  export "PYTHONPATH=$PYTHONPATH:~/where/is/your/modelBC03/"

readAtlasImages
  http://www.sdss.org/dr12/algorithms/read_psf/
  readAtlasImages-v5_4_11.tar.gz
within ~/.bashrc please add:
  DIRREADATLAS='~/where/is/your/readAtlas/readAtlasImages-v5_4_11/'
  export DIRREADATLAS
  alias dir_readatlas="$DIRREADATLAS"

humvi
  https://github.com/drphilmarshall/HumVI
within ~/.bashrc please add:
  HUMVI='~/where/is/your/HumVI-master/compose.py'
  export HUMVI
  alias humvi="$HUMVI"

"""

setup(name='bubbleimg',
      version='1.0',
      description='Produce emission-line images from broadband images',
      url='http://github.com/aileisun/bubbleimg',
      author='Ai-Lei Sun',
      author_email='aileisundr@gmail.com',
      license='MIT',
      packages=['bubbleimg'],
      install_requires=[
          'numpy',
          'matplotlib',
          'astropy',
          'scipy',
          'astroquery',
          'scikit-image',
          'shapely',
          'PyAstronomy',
          'paramiko',
          'requests',
      ],
      zip_safe=False)