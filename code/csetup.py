#!/usr/bin/env python

#instructions: python csetup.py build
 
from distutils.core import setup
from distutils.extension import Extension
 
setup(name="h_combinations",
    ext_modules=[
        Extension("h_combinations", ["combinations.cpp"],
        libraries = ["boost_python"])
    ])

