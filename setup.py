#!/usr/bin/env python

from setuptools import setup

# Utility function to read the README file.  
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

with open("README.md", 'r') as f:
    long_description = f.read()


setup(
    name='SNPs_Analysis',
    version='0.0.3',
    description='SNPs_Analysis is a pipeline analysis bioinformatics workflows\
     that is specifically designed for functional genomic scientists or researchers in the field of biomedical or bioinformatics. This analysis pipeline is supported with high performance computing systems (HPCs). That is, the HPC is for running pipeline stages on a distributed compute nodes.',
    author='Ajayi Olabode',
    author_email='boraton2010@gmail.com',
    packages=['src'],
    entry_points={
        'console_scripts': ['SNPs_Analysis = src.main:main']
    },
    url='https://github.com/boratonAJ/SNPs_Analysis',
    license='LICENSE.txt',
    keywords = "example documentation tutorial",
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
     ],


    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
    long_description = """\
Universal character encoding detector
-------------------------------------

Detects
 - ASCII, UTF-8, UTF-16 (2 variants), UTF-32 (4 variants)
 - Big5, GB2312, EUC-TW, HZ-GB-2312, ISO-2022-CN (Traditional and Simplified Chinese)
 - EUC-JP, SHIFT_JIS, ISO-2022-JP (Japanese)
 - EUC-KR, ISO-2022-KR (Korean)
 - KOI8-R, MacCyrillic, IBM855, IBM866, ISO-8859-5, windows-1251 (Cyrillic)
 - ISO-8859-2, windows-1250 (Hungarian)
 - ISO-8859-5, windows-1251 (Bulgarian)
 - windows-1252 (English)
 - ISO-8859-7, windows-1253 (Greek)
 - ISO-8859-8, windows-1255 (Visual and Logical Hebrew)
 - TIS-620 (Thai)

This version requires Python 3 or later; a Python 2 version is available separately.
"""





)
