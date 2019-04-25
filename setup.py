import setuptools

setuptools.setup(
    name='dnabc',
    version="0.0.5",
    description='Demultiplex pooled DNA sequencing data',
    author='Kyle Bittinger',
    author_email='kylebittinger@gmail.com',
    url='https://github.com/PennChopMicrobiomeProgram',
    packages=['dnabclib'],
    entry_points={
        'console_scripts': [
            'dnabc.py=dnabclib.main:main',
            'split_samplelanes.py=dnabclib.split_samplelanes:main',
        ],
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    license='GPLv2+',
)
