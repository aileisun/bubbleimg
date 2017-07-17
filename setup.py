from setuptools import setup

"""
  To use the hsc part one needs to be in lsst environment
  # https://pipelines.lsst.io/install/conda.html
  source activate lsst
  source eups-setups.sh
  setup lsst_distrib
"""

setup(name='bubbleimg',
      version='0.1',
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
          'os',
          'scipy',
          'sys',
          'astroquery',
          'pywt',
          'copy',
          'collections',
          'scikit-image',
          'shapely',
          'PyAstronomy',
          'paramiko',
          'requests',

      ],
      zip_safe=False)