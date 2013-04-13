from distutils.core import setup
from distutils.extension import Extension

setup(name="p4",
    ext_modules=[
        Extension("p4stm", ['p4stm.cpp'],
            # Adjust the following to be able to find pyublas, numpy, and boost.  Maybe need library_dirs as well?
            include_dirs = ["/Users/peter/lib/python/PyUblas-2011.1-py2.6-macosx-10.6-x86_64.egg/include",
              "/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/numpy/core/include",
              "/usr/local/include"],
                  libraries = ["boost_python-mt"]) 
    ])
