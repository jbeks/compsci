#!/usr/bin/env python3
from distutils.core import setup, Extension

# set module parameters
system_interpolation_module = Extension(
    "system_interpolation",
    sources=["system_interpolation.cpp"],
    include_dirs=["."],
    language="c++"
)

# create module
setup(
    name="system_interpolation",
    version="0.1.0",
    description="Interpolation functions for nbody simulation",
    ext_modules=[system_interpolation_module]
)

