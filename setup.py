from setuptools import setup, find_packages

setup()

'''
exec( open('rapyuta/version.py').read() )
__version__ = version

setup(
    name = 'rapyuta',
    author = 'D. HU',
    author_email = 'dangning.hu@outlook.com',
    description = 'Regular AstroPY UTility Assemblage',
    long_description = 'file: README.org',
    license = 'BSD 3-Clause License',
    keywords = 'astronomy, astrophysics, astrometry, imaging, photometry, spectroscopy, space, spitzer, akari, jwst',
    url = 'https://github.com/kxxdhdn/RAPYUTA',
    project_urls={
        'Reports': 'https://github.com/kxxdhdn/RAPYUTA/issues',
        'Source': 'https://github.com/kxxdhdn/RAPYUTA/rapyuta',
    },

    ## Versions
    version=__version__,
    # setup_requires=['setuptools_scm'],
    # use_scm_version=True,
    ## in case setup.py is not in Git root
    # use_scm_version={
    #     "root": "..", # Git root
    #     "relative_to": __file__, # setup.py
    # },

    python_requires='>=3.8',
    install_requires = [
        'astropy>=5.1',
        'colorama',
        'h5py',
        'matplotlib', 
        'numpy',
        'photutils',
        'reproject>=0.7.1',
        'scipy',
        'specutils',
        'tqdm',
    ],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Fortran',
        'Programming Language :: IDL',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    
    ## Plugins
    entry_points={
        # Installation test with command line
        'console_scripts': [
            'rapyboost = rapyuta:boost_rapyuta',
        ],
    },

    ## Packages
    ## Set packages to find: to automatically find all sub-packages
    packages = find_packages(),
    exclude = [
        'rapyuta._dev',
    ],
    ## Set packages manually
    # packages = [
    #     'rapyuta',
    #     'rapyuta.impro',
    # ],

    ## Package data
    package_data = {
        # include files in rapyuta/lib
        'rapyuta': [
            'lib/*.h5',
            'lib/*.txt',
            'lib/filt/*.h5',
            'lib/filt/*.txt',
            'lib/logo/*.txt',],
    },
)
'''