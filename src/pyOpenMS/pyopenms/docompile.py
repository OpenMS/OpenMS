# Python script to compile Cython to C++
#
# this script compiles the *.pyx infile (pyopenms/pyopenms.pyx) to .cpp
#
# It can be run in the build process of pyopenms when manual changes to the pyx
# files need to be made (not recommended!).
#
# == taken from autowrap/Main.py run(), lower part
from __future__ import print_function

infile  = "pyopenms/pyopenms.pyx"

import cPickle
persisted_data_path = "include_dir.bin"
autowrap_include_dirs = cPickle.load(open(persisted_data_path, "rb"))
# autowrap_include_dirs = ['/usr/local/lib/python2.7/dist-packages/autowrap/data_files/boost', '/usr/local/lib/python2.7/dist-packages/autowrap/data_files', '/path/to/pyOpenMS/pxds', 'from Map cimport Map as _Map']


from Cython.Compiler.Main import compile as cy_compile, CompilationOptions
from Cython.Compiler.Options import directive_defaults
directive_defaults["boundscheck"] = False
directive_defaults["wraparound"] = False

options = dict(include_path=autowrap_include_dirs,
               compiler_directives=directive_defaults,
               #output_dir=".",
               #gdb_debug=True,
               cplus=True)

print("Compiling with Cython the file", infile)
print("Using include_path", autowrap_include_dirs)
options = CompilationOptions(**options)
cy_compile(infile, options=options)
print("Success!")

