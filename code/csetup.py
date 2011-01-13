#!/usr/bin/env python

"""
instructions: python csetup.py build
sudo apt-get install python-dev
sudo apt-get install libcgal-dev
sudo apt-get install libboost-python-dev 
sudo apt-get install 
"""

#these 2.5M combinations should be written to "test.out" in a couple of seconds
#compared to the itertools.combinations there is mainly a speed advantage if the
#strings to be written are short
 
from distutils.core import setup
from distutils.extension import Extension
 
setup(name="srmcollider",
    url = "http://www.srmcollider.org", 
    version = "0.7",
    author = "Hannes Roest",
    ext_modules=[
        Extension("c_combinations", ["combinations.cpp"],
            libraries = ["boost_python"]),
        Extension("c_getnonuis", ["getNonUis.cpp"],
            libraries = ["boost_python"]),
        Extension("r_rangetree", ["rangetree.cpp"],
            libraries = ["boost_python", "CGAL"]),
    ]) 


