import setuptools

setuptools.setup(
    name='dnabc',
    version="0.0.4",
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
)
