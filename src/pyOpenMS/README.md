General
-------

pyopenms is a Python library for the analysis of mass spectrometry data.
It is mainly based on Cython wrappers around the OpenMS C++ library. To see which classes and functions are
currently wrapped, please check the pxd files under "./pxds" or consult our
[API documentation](https://pyopenms.readthedocs.io/en/latest/apidocs/index.html).
Additionally, it provides some convenience functions for plotting or converting from/to dataframes or numpy arrays.

Wrapping new classes
--------------------

See [README_WRAPPING_NEW_CLASSES](./README_WRAPPING_NEW_CLASSES)

Build instructions
------------------

0. (optional) Create a virtual python environment:
    
    ```bash
    python -m venv /path/to/myenv
    
    # ... and activate it, e.g.
    # Linux:
    source <venv>/bin/activate
    # Windows:
    c:\path\to\myenv\Scripts\activate.bat
    ```
    
1. Get Python 3.7+ and the following Python libraries:
   
   pip install -r requirements_bld.txt

2. If running from an OpenMS build tree (recommended), just reconfigure with

   ```bash
   cmake -DPYOPENMS=ON .
   ```
   
   If it does not find the python that you just installed or complains about libraries not located although
   you just installed them, help it find the correct python executable by adding `-DPython_EXECUTABLE="/path/to/python(.exe)"`
   If your computer has a lot of RAM (16GB+) you can add `-DPY_NUM_THREADS=2`
   (or up to the number of split modules, which are by default PY_NUM_MODULES=8).

   Building with an existing, installed OpenMS library is possible but not well-tested. All you need is the current pyOpenMS
   directory (where this README is located). Configure CMake in a build dir of your choice with general CMake
   and pyopenms-related options only. It will try to find OpenMS based on its CMake config files that should be installed
   with newer versions of OpenMS (around 2.8+) and hopefully parse all necessary options and find OpenMS' transitive dependencies.
   The rest is the same.

4. Build CMake target "pyopenms" build-system agnostic with

   ```bash
   cmake --build . --target pyopenms
   ```

5. Run tests with

   ```bash
   ctest -R pyopenms
   ```

   "-R" to restrict to pyopenms* tests. If running out of the OpenMS build tree, this should not be necessary.

6. Install locally (and in-place for live edits [option -e]) into current Python with

   ```bash
   pip install -e pyopenms --no-cache-dir --no-binary=pyopenms
   ```
   
   `--no-binary` is used because the binaries are/were built with CMake.
