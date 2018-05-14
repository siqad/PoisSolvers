#!/usr/bin/env python

'''
setup.py file for SWIG phys connector
'''

from distutils.core import setup, Extension
import os
os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"
phys_connector_module = Extension('_phys_connector',
                                sources=['phys_connector_wrap.cxx', 'phys_connector.cc'],
                                )

setup (
        name    = 'phys_connector',
        version = '0.01',
        author  = 'Samuel Ng',
        description = '''Python wrapper for physics connector''',
        ext_modules = [phys_connector_module],
        py_modules = [],
    )
