#!/usr/bin/env python3

#########################################################################
# Author: Hechuan Yang
# Created Time: 2017-06-21 16:14:38
# File Name: setup.py
# Description: 
#########################################################################

"""setup.py: setuptools control."""
 
import re
from setuptools import setup
 
 
version = re.search(
    "^__version__\s*=\s*'(.*)'",
    open('psite/psite.py').read(),
    re.M
    ).group(1)
 
with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")
 
setup(
    name = "PSiTE",
    packages = ["psite"],
    install_requires=['numpy','pyfaidx','pyyaml'],
    entry_points = {
        "console_scripts": ['psite = psite.psite:main']
        },
    version = version,
    description = "PSiTE is a Python program for jointly simulating Single Nucleotide Variants (SNVs) and Copy Number Variants (CNVs) for a sample of tumor cells.",
    long_description = long_descr,
    author = "Hechuan Yang",
    author_email = "yanghch@outlook.com",
    url = "https://github.com/hchyang/PSiTE.git",
#https://docs.python.org/3/library/subprocess.html
#The subprocess.run() function was added in Python 3.5
    python_requires='>=3.5'
    )


