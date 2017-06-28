#!/usr/bin/env python

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
    open('csite/csite.py').read(),
    re.M
    ).group(1)
 
with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")
 
setup(
    name = "CSiTE",
    packages = ["csite"],
    install_requires=['numpy'],
    entry_points = {
        "console_scripts": ['csite = csite.csite:main']
        },
    version = version,
    description = "CSiTE is a Python program for jointly simulating Single Nucleotide Variants (SNVs) and Copy Number Variants (CNVs) for a sample of tumor cells.",
    long_description = long_descr,
    author = "Hechuan Yang",
    author_email = "yanghechuan@gmail.com",
    url = "https://github.com/hchyang/CSiTE.git",
    python_requires='>=3'
    )


