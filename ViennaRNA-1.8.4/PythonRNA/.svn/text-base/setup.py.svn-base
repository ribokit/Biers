#!/usr/bin/env python

"""
setup.py file for SWIG RNA
"""

from distutils.core import setup, Extension


RNA_module = Extension('_RNA',
                           sources=['RNA_wrap.c'],
			   library_dirs=[ '../'],
			   libraries=['RNA'],
                           )

setup (name = 'RNA',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig RNA from docs""",
       ext_modules = [RNA_module],
       py_modules = ["RNA"],
       )
