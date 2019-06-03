#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='pywigxjpf',
      version='1.9',
      description='Wigner Symbols',
      license = "GPLv3",
      author='C. Forssen and H. T. Johansson',
      url='http://fy.chalmers.se/subatom/wigxjpf',
      packages=['pywigxjpf'],
      package_data={'pywigxjpf': ['../lib/libwigxjpf_shared.*']},
      setup_requires=["cffi>=1.0.0"],
      cffi_modules=["pywigxjpf/pywigxjpf_ffi_builder.py:ffibuilder"],
      install_requires=["cffi>=1.0.0"],
     )
