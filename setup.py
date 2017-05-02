from setuptools import setup


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
          'pickle',
          'skimage',
          'shapely',
          'PyAstronomy'
      ],
      zip_safe=False)