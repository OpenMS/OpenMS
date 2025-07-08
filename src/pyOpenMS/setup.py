
#!/usr/bin/env python

# THIS IS ADAPTED FROM PYARROW
# https://github.com/apache/arrow/blob/main/python/setup.py


# All the package metadata and build+run dependencies goes to
# pyproject.toml. This setup.py is only for building/registering
# the C++ extensions.

# Idea is:
#  - if run outside of the general OpenMS CMake script,
#  the setup.py just calls cmake --build pyopenms, collects the
#  extension modules and declares them as output
#  - if run inside the general OpenMS CMake script, the setup.py
#  reads the generated env.py and uses that to add the
#  modules to be expected from the build to ext_modules



# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

import contextlib
import os
import os.path
from os.path import join as pjoin
import re
import shlex
import sys
import warnings
try:
    import env
except ImportError:
    print("Info: env.py not found. Which means running outside of a CMake process.")

if sys.version_info >= (3, 10):
    import sysconfig
else:
    # Get correct EXT_SUFFIX on Windows (https://bugs.python.org/issue39825)
    from distutils import sysconfig

from setuptools import setup, Extension, Distribution

from Cython.Distutils import build_ext as _build_ext
import Cython

# Check if we're running 64-bit Python
is_64_bit = sys.maxsize > 2**32


if Cython.__version__ < '3':
    raise Exception(
        'Please update your Cython version. Supported Cython >= 3')

setup_dir = os.path.abspath(os.path.dirname(__file__))

ext_suffix = sysconfig.get_config_var('EXT_SUFFIX')


@contextlib.contextmanager
def changed_dir(dirname):
    oldcwd = os.getcwd()
    os.chdir(dirname)
    try:
        yield
    finally:
        os.chdir(oldcwd)


def strtobool(val):
    """Convert a string representation of truth to true (1) or false (0).

<<<<<<< HEAD
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    # Copied from distutils
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
=======
if iswin:
    # /EHs is important. It sets _CPPUNWIND which causes boost to
    # set BOOST_NO_EXCEPTION in <boost/config/compiler/visualc.hpp>
    # such that  boost::throw_excption() is declared but not implemented.
    # The linker does not like that very much ...
    extra_compile_args = ["/EHs", "/bigobj"]
    extra_compile_args.append("/std:c++17")
    extra_link_args.append("/std:c++17")

elif sys.platform.startswith("linux"):
    extra_link_args = ["-Wl,-s"]
    if OMP:
        libraries.append("libomp")
        libraries.append("pthread")
elif sys.platform == "darwin":
    library_dirs.insert(0,j(OPEN_MS_BUILD_DIR,"pyOpenMS","pyopenms"))
    if OMP:
        libraries.append("omp")
    # we need to manually link to the Qt Frameworks
    python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    pyopenms_path = f"@loader_path/../lib/python{python_version}/site-packages/pyopenms"
    extra_compile_args = ["-Qunused-arguments", "-fopenmp"]
    extra_link_args = ["-Wl,-rpath,@loader_path/../lib", f"-Wl,-rpath,{pyopenms_path}", "-lomp"]
if IS_DEBUG:
    extra_compile_args.append("-g2")
if OMP and OPENMP_CXX_FLAGS:
    extra_compile_args.extend(OPENMP_CXX_FLAGS.split(";"))

if not iswin:
    extra_link_args.append("-std=c++17")
    extra_compile_args.append("-std=c++17")
    if isosx: # MacOS
        extra_compile_args.append("-stdlib=libc++")
        extra_link_args.append("-stdlib=libc++") # MacOS libstdc++ does not include c++11+ lib support.
        extra_link_args.append("-mmacosx-version-min=10.9") # due to libc++
        extra_compile_args.append("-Wno-deprecated")
        extra_compile_args.append("-Wno-nullability-completeness")
        if (osx_ver >= "10.14.0" and SYSROOT_OSX_PATH): # since macOS Mojave
            extra_link_args.append("-isysroot" + SYSROOT_OSX_PATH)
            extra_compile_args.append("-isysroot" + SYSROOT_OSX_PATH)
>>>>>>> develop
    else:
        raise ValueError("invalid truth value %r" % (val,))


class build_ext(_build_ext):
    _found_names = ()

    CYTHON_MODULE_NAMES = [ 
        '_all_modules',
        '_pyopenms_1',
        '_pyopenms_2',
        '_pyopenms_3',
        '_pyopenms_4',
        '_pyopenms_5',
        '_pyopenms_6',
        '_pyopenms_7',
        '_pyopenms_8',
        '_version'
    ]

    def build_extensions(self):
        import numpy
        numpy_incl = numpy.get_include()

        self.extensions = [ext for ext in self.extensions
                           if ext.name != '__dummy__']

        for ext in self.extensions:
            if (hasattr(ext, 'include_dirs') and
                    numpy_incl not in ext.include_dirs):
                ext.include_dirs.append(numpy_incl)
        _build_ext.build_extensions(self)

    def run(self):
        self._run_cmake()
        _build_ext.run(self)

    # adapted from cmake_build_ext in dynd-python
    # github.com/libdynd/dynd-python

    description = "Build the Cpp-extensions for OpenMS"
    user_options = ([('cmake-generator=', None, 'CMake generator'),
                     ('extra-cmake-args=', None, 'extra arguments for CMake'),
                     ('num_modules=', None, 'number of modules to split cython cpps into'),
                     ('num_threads=', None, 'number of threads to use for building, i.e. how many modules to build in parallel'),
                     ('build-type=', None, 'build type (debug or release), default release')] +
                    _build_ext.user_options)

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.cmake_generator = os.environ.get('OPENMS_CMAKE_GENERATOR')
        if not self.cmake_generator and sys.platform == 'win32':
            self.cmake_generator = 'Visual Studio 15 2017 Win64'
        self.extra_cmake_args = os.environ.get('OPENMS_CMAKE_OPTIONS', '')
        self.build_type = os.environ.get('OPENMS_BUILD_TYPE',
                                         'release').lower()

        self.cmake_cxxflags = os.environ.get('OPENMS_CXXFLAGS', '')

        if sys.platform == 'win32':
            # Cannot do debug builds in Windows unless Python itself is a debug
            # build
            if not hasattr(sys, 'gettotalrefcount'):
                self.build_type = 'release'


    def _run_cmake(self):
        # check if build_type is correctly passed / set
        if self.build_type.lower() not in ('release', 'debug',
                                           'relwithdebinfo'):
            raise ValueError("--build-type (or OPENMS_BUILD_TYPE) needs to "
                             "be 'release', 'debug' or 'relwithdebinfo'")

        # The directory containing this setup.py
        source = os.path.dirname(os.path.abspath(__file__)) + "/../../"

        # The staging directory for the module being built
        build_cmd = self.get_finalized_command('build')
        saved_cwd = os.getcwd()
        build_temp = pjoin(saved_cwd, build_cmd.build_temp)
        build_lib = pjoin(saved_cwd, build_cmd.build_lib)

        if not os.path.isdir(build_temp):
            self.mkpath(build_temp)

        if self.inplace:
            # a bit hacky
            build_lib = saved_cwd

        install_prefix = pjoin(build_lib, "OPENMS")

        # Change to the build directory
        with changed_dir(build_temp):
            # Detect if we built elsewhere
            if os.path.isfile('CMakeCache.txt'):
                cachefile = open('CMakeCache.txt', 'r')
                cachedir = re.search('CMAKE_CACHEFILE_DIR:INTERNAL=(.*)',
                                     cachefile.read()).group(1)
                cachefile.close()
                if (cachedir != build_temp):
                    build_base = pjoin(saved_cwd, build_cmd.build_base)
                    print(f"-- Skipping build. Temp build {build_temp} does "
                          f"not match cached dir {cachedir}")
                    print("---- For a clean build you might want to delete "
                          f"{build_base}.")
                    return

            cmake_options = [
                f'-DCMAKE_INSTALL_PREFIX={install_prefix}',
                f'-DPYTHON_EXECUTABLE={sys.executable}',
                f'-DPython3_EXECUTABLE={sys.executable}',
                f'-DOPENMS_CXXFLAGS={self.cmake_cxxflags}',
            ]

            def append_cmake_bool(value, varname):
                cmake_options.append('-D{0}={1}'.format(
                    varname, 'on' if value else 'off'))

            def append_cmake_component(flag, varname):
                # only pass this to cmake if the user pass the --with-component
                # flag to setup.py build_ext
                if flag is not None:
                    flag_name = (
                        "--with-"
                        + varname.removeprefix("OPENMS_").lower().replace("_", "-"))
                    warnings.warn(
                        MSG_DEPR_SETUP_BUILD_FLAGS.format(flag_name),
                        UserWarning, stacklevel=2
                    )
                    append_cmake_bool(flag, varname)

            if self.cmake_generator:
                cmake_options += ['-G', self.cmake_generator]

            #append_cmake_bool(self.bundle_cython_cpp,
            #                  'OPENMS_BUNDLE_CYTHON_CPP')
            #append_cmake_bool(self.generate_coverage,
            #                  'OPENMS_GENERATE_COVERAGE')

            cmake_options.append(
                f'-DCMAKE_BUILD_TYPE={self.build_type.lower()}')

            extra_cmake_args = shlex.split(self.extra_cmake_args)

            build_tool_args = []
            if sys.platform == 'win32':
                if not is_64_bit:
                    raise RuntimeError('Not supported on 32-bit Windows')
            else:
                build_tool_args.append('--')
                if os.environ.get('OPENMS_BUILD_VERBOSE', '0') == '1':
                    cmake_options.append('-DCMAKE_VERBOSE_MAKEFILE=ON')
                parallel = os.environ.get('OPENMS_PARALLEL')
                if parallel:
                    build_tool_args.append(f'-j{parallel}')

            # Generate the build files
            print("-- Running cmake for OPENMS")
            self.spawn(['cmake'] + extra_cmake_args + cmake_options + [source])

            print("-- Finished cmake for OPENMS")

            print("-- Running cmake --build for OPENMS")
            print(f"Current working directory: {os.getcwd()}")
            self.spawn(['cmake', '--build', source, '--config', self.build_type] +
                       build_tool_args)
            print("-- Finished cmake --build for OPENMS")

            print("-- Running cmake --build --target install for OPENMS")
            self.spawn(['cmake', '--build', source, '--config', self.build_type] +
                       ['--target', 'install'] + build_tool_args)
            print("-- Finished cmake --build --target install for OPENMS")

            self._found_names = []
            for name in self.CYTHON_MODULE_NAMES:
                built_path = pjoin(install_prefix, name + ext_suffix)
                if os.path.exists(built_path):
                    self._found_names.append(name)

    def _get_build_dir(self):
        # Get the package directory from build_py
        build_py = self.get_finalized_command('build_py')
        return build_py.get_package_dir('OPENMS')

    def _get_cmake_ext_path(self, name):
        # This is the name of the arrow C-extension
        filename = name + ext_suffix
        return pjoin(self._get_build_dir(), filename)

    def get_ext_generated_cpp_source(self, name):
        if sys.platform == 'win32':
            head, tail = os.path.split(name)
            return pjoin(head, tail + ".cpp")
        else:
            return pjoin(name + ".cpp")

    def get_ext_built_api_header(self, name):
        if sys.platform == 'win32':
            head, tail = os.path.split(name)
            return pjoin(head, tail + "_api.h")
        else:
            return pjoin(name + "_api.h")

    def get_names(self):
        return self._found_names

    def get_outputs(self):
        # Just the C extensions
        # regular_exts = _build_ext.get_outputs(self)
        return [self._get_cmake_ext_path(name)
                for name in self.get_names()]


class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True


setup(
    distclass=BinaryDistribution,
    # Dummy extension to trigger build_ext
    ext_modules=[Extension('__dummy__', sources=[])],
    cmdclass={
        'build_ext': build_ext
    },
)
