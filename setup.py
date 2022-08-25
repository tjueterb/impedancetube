# coding=UTF-8
#------------------------------------------------------------------------------
# Copyright (c) Acoular Development Team.
#------------------------------------------------------------------------------

from setuptools import setup, find_packages
from os.path import join, abspath, dirname
import os

bf_version = "0.2.5"
bf_author = "Tobias Jüterbock"
bf_email = "a.jueterbock@tu-berlin.de"


# Get the long description from the relevant file
here = abspath(dirname(__file__))
with open(join(here, 'README.md')) as f:
    long_description = f.read()


install_requires = list([
    'acoular',
    'numpy',
    'setuptools',
    'acoular',
    'traits>=6.0',
])

setup_requires = list([
    'acoular',
    'numpy',
    'setuptools',
    'acoular',
    'traits>=6.0',
])

setup(name="impedancetube",
      version=bf_version,
      description="Library for Impedance/Transmission tube evaluation with the transfer function method",
      long_description_content_type='text/markdown',
      long_description=long_description,
      url='https://github.com/tjueterb/impedancetube',
      license="BSD",
      author=bf_author,
      author_email=bf_email,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Physics',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
      ],
      keywords='acoustic impedance transmission tube microphone measurement',
      packages=['impedancetube','tests'],
    #   package_dir={'': 'src'},

      install_requires=install_requires,

      setup_requires=setup_requires,

      include_package_data=True,
      package_data={
		    'impedancetube': ['tests/*.*']},
      #to solve numba compiler
      zip_safe=False
      )
