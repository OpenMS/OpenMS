#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
--------------------------------------------------------------------------
                  OpenMS -- Open-Source Mass Spectrometry
--------------------------------------------------------------------------
Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
ETH Zurich, and Freie Universitaet Berlin 2002-2020.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
For a full list of authors, refer to the file AUTHORS.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

