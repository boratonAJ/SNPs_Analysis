#!/usr/bin/env python

from setuptools import setup

# Utility function to read the README file.  
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

with open("README.md", 'r') as f:
    long_description = f.read()


setup(
    name='SNPS_Analysis',
    version='0.0.1',
    description='complexo_pipeline is a pipeline system for bioinformatics workflows\
     with support for running pipeline stages on a distributed compute cluster.',
    author='Ajayi Olabode',
    author_email='boraton2010@gmail.com',
    packages=['src'],
    entry_points={
        'console_scripts': ['SNPs_Analysis = src.main:main']
    },
    url='https://github.com/boratonAJ/SNPs_Analysis',
    license='LICENSE.txt',
    keywords = "example documentation tutorial",
    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
