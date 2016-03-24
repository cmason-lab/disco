from distutils.core import setup

setup(
      name='disco',
      version='0.1.0',
      description='Python tools for analyzing single cell alternative splicing',
      author='Priyanka Vijay',
      author_email='prv2004',
      url='',
      packages=[
                'disco',
                ],
      package_dir={'disco': 'lib/disco'},
      scripts=['scripts/disco'],
      requires=['numpy', 'pandas', 'scipy', 'matplotlib', 'seaborn']
     )