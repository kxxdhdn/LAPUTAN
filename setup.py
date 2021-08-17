from setuptools import setup

setup(
    name = 'laputan',
    version = '1.2.3',
    author = 'D. HU',
    author_email = 'dangning.hu@cea.fr',
    description = 'Library of Astronomical Python Utility Tool for Astrophysics Nerds',
    license = 'BSD',
    keywords = 'astronomy astrophysics',
    url = 'https://github.com/kxxdhdn/laputan',
    project_urls={
        'IDL': 'https://github.com/kxxdhdn/laputan/tree/master/idl',
        'SwING': 'https://github.com/kxxdhdn/laputan/tree/master/swing',
        'Tests': 'https://github.com/kxxdhdn/laputan/tree/master/tests',
    },

    python_requires='>=3.6',
    install_requires = [
        'numpy', 'scipy', 'matplotlib', 
        'astropy', 'reproject', 'h5py', 'tqdm', 
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
            'laputest = laputan:iTest',
        ],
    },

    ## Packages
    packages = ['laputan'],

    ## Package data
    package_data = {
        # include files in laputan/lib
        'laputan': ['lib/*.txt','lib/data/*.h5'],
    },
)
