from distutils.core import setup

setup(
      name='disco',
      version='0.0.1',
      description='Python tools for analyzing single cell alternative splicing',
      author='Priyanka Vijay',
      author_email='prv2004',
      url='',
      package_dir={'': 'lib'},
      entry_points={
                    'console_scripts': [
                                        'disco = disco.disco:main',
                                        'Disco = disco.disco:main',
                                        ]
                   },
      packages=[
                'disco',
                ],
      requires=['numpy', 'pandas', 'scipy', 'matplotlib', 'seaborn']
     )
