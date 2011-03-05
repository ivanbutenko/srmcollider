#!/usr/bin/env python

"""
Dependencies in Python (the following packages must be available):
    MySQLdb
    sqlite

instructions: python csetup.py build
sudo apt-get install python-dev
sudo apt-get install libboost-dev
sudo apt-get install libcgal-dev
sudo apt-get install libboost-python-dev 

#or all at the same time
sudo apt-get install -y python-dev libcgal-dev libboost-python-dev


On Ubuntu, it might fail with 
 /usr/bin/ld: cannot find -lboost_python
then try to do this
sudo ln -s /usr/lib/libboost_python-mt.so /usr/lib/libboost_python.so


gcc -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I/usr/local/include/python2.7 -c combinations.cpp -o build/temp.linux-x86_64-2.7/combinations.o
g++ -pthread -shared build/temp.linux-x86_64-2.7/combinations.o -lboost_python -o build/lib.linux-x86_64-2.7/c_combinations.so


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
    requires=["MySQLdb", "sqlite"],
    ext_modules=[
        Extension("c_combinations", ["combinations.cpp"],
            libraries = ["boost_python"]),
        Extension("c_getnonuis", ["getNonUis.cpp"],
            libraries = ["boost_python"]),
        Extension("c_rangetree", ["rangetree.cpp"],
            libraries = ["boost_python", "CGAL"]),
        Extension("c_integrated", ["integratedrun.cpp"],
            libraries = ["boost_python", "CGAL"]),
        Extension("libsrmcolliderLib.a", ["srmcolliderLib.cpp"],
            libraries = ["boost_python" ]),
    ],
     ) 


