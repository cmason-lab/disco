#!/usr/bin/env python

# from distutils.core import setup
from setuptools import setup

setup(
      name='disco',
      version='0.3',
      description='Python tools for analyzing single cell alternative splicing',
      author='Priyanka Vijay',
      author_email='prv2004@med.cornell.edu',
      url='https://pbtech-vc.med.cornell.edu/git/mason-lab/disco.git',
      packages=['disco'],
      entry_points={'console_scripts': ['disco = disco.run:main']},
      requires=['numpy', 'pandas', 'scipy', 'matplotlib', 'seaborn', 'statsmodels']
     )
#      scripts=['scripts/Disco'],
#      package_dir={'': 'lib'},
#scripts=['scripts/Disco'],