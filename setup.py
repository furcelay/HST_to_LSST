import os
from setuptools import setup, find_packages


setup(
    name='hsc_to_lsst',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'astropy',
        'photutils',
        'reproject',
        'lenstronomy',
        'scipy',
        'scikit-learn',
        'tqdm'
    ],
    # entry_points={
    #     'console_scripts': [
    #         'hsc_to_lsst = hsc_to_lsst.__main__:main'
    #     ]
    # },
    package_data={
        'hsc_to_lsst': ['data/*']
    },
    author='Felipe Urcelay, LSST Strong Lensing Science Collaboration',
)