#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
SPDX-License-Identifier: BSD-3-Clause
--------------------------------------------------------------------------
$Maintainer: Hannes Roest $
$Authors: Hannes Roest $
--------------------------------------------------------------------------
"""

"""
This program can be used to only compile the pyopenms.pyx to a cpp file using
Cython. This is an intermediate step in the pyOpenMS build process and it
should only be used for debugging.

It may be useful if something went wrong during the Cython compilation. One can
then try to fix the pyopenms.pyx file, run this script and see whether Cython
parses it.

Author: Hannes Roest
"""

# Py2/3 fix
try:
    import cPickle
except ImportError:
    import _pickle as cPickle

persisted_data_path = "include_dir.bin"
try:
    autowrap_include_dirs = cPickle.load(open(persisted_data_path, "rb"))
except IOError:
    print("The file include_dir.bin does not yet exist, please run setup.py first to create it.")

from Cython.Compiler.Main import compile, CompilationOptions
import Cython
print("Will try to compile with Cython version", Cython.__version__)

# Prepare options
print ("include:", autowrap_include_dirs)
options = dict(include_path=autowrap_include_dirs,
               #output_dir=".",
               #gdb_debug=True,
               cplus=True)

# Do Cython compile
import sys
out  = sys.argv[1]
options = CompilationOptions(**options)
compile(out, options=options)

