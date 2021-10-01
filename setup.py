from setuptools import setup

setup(
    name = 'rapyuta',
    version = '2.0.2',
    author = 'D. HU',
    author_email = 'dangning.hu@outlook.com',
    description = 'libraRy of Astronomical PYthon UTilities for Astrophysics nerds',
    license = 'BSD',
    keywords = 'astronomy astrophysics astrometry imaging spectroscopy spitzer akari jwst',
    url = 'https://github.com/kxxdhdn/RAPYUTA',
    project_urls={
        'IDL': 'https://github.com/kxxdhdn/RAPYUTA/tree/main/idl',
        'SwING': 'https://github.com/kxxdhdn/RAPYUTA/tree/main/swing',
        'Tests': 'https://github.com/kxxdhdn/RAPYUTA/tree/main/tests',
    },

    python_requires='>=3.6',
    install_requires = [
        'numpy', 'scipy', 'matplotlib', 
        'astropy', 'reproject>=0.7.1', 'h5py', 'tqdm', 
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
    ],
    
    ## Plugins
    entry_points={
        # Installation test with command line
        'console_scripts': [
            'rapyutest = rapyuta:iTest',
        ],
    },

    ## Packages
    packages = ['rapyuta'],

    ## Package data
    package_data = {
        # include files in rapyuta/lib
        'rapyuta': ['lib/*.txt','lib/data/*.h5'],
    },
)
