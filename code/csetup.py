#!/usr/bin/env python

#instructions: python csetup.py build
#to see whether it works:

"""
import string
import c_combinations
c_combinations.combinations(5, len(string.letters), 
                            'test.out', [l for l in string.letters] )
"""

#these 2.5M combinations should be written to "test.out" in a couple of seconds
#compared to the itertools.combinations there is mainly a speed advantage if the
#strings to be written are short
 
from distutils.core import setup
from distutils.extension import Extension
 
setup(name="c_combinations",
    ext_modules=[
        Extension("c_combinations", ["combinations.cpp"],
        libraries = ["boost_python"])
    ])

