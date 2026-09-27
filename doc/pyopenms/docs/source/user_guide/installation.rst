Installation
============


Command Line
------------

To :index:`install` pyOpenMS from the command line using the binary wheels, you
can type

.. code-block:: bash

  pip install numpy
  pip install pyopenms


We have binary packages for OSX, Linux and Windows (64 bit only) available from
`PyPI <https://pypi.org/project/pyopenms>`_. Make sure to download
the 64bit Python release for Windows. Currently we only support
Python 3.11, 3.12, 3.13 and 3.14.

You can install Python first from `here <https://www.python.org/downloads/>`_,
again make sure to download the 64bit release. You can then open a shell and
type the two commands above (on Windows you may potentially have to use
``C:\Python37\Scripts\pip.exe`` in case ``pip`` is not in your system path).

Nightly/ CI wheels
------------------

If you want the newest features you can also install nightly builds of pyOpenMS with the following shell command:

.. code-block:: bash

  pip install --index-url https://pypi.openms.de/simple/ pyopenms

Type Stubs
----------

Since version 3.6, the pyOpenMS wheels contain type stubs (``.pyi`` files with
the parameter types, return types and docstrings of the classes and functions)
and the ``py.typed`` marker (:pep:`561`). Editors use the stubs for code
completion and to show method signatures, and type checkers such as mypy or
pyright use them to find wrong types before the code runs; no further setup is
needed. For example, for this script ``check_types.py``:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    exp = oms.MSExperiment()
    oms.MzMLFile().load(42, exp)
    n_spectra: str = exp.getNrSpectra()

mypy, installed in the same environment as pyOpenMS (``pip install mypy``),
reports both errors:

.. code-block:: bash

  mypy check_types.py

.. code-block:: output

    check_types.py:4: error: Argument 1 to "load" of "MzMLFile" has incompatible type "int"; expected "str"  [arg-type]
    check_types.py:5: error: Incompatible types in assignment (expression has type "int", variable has type "str")  [assignment]
    Found 2 errors in 1 file (checked 1 source file)

Source (advanced users)
-----------------------

To install pyOpenMS from :index:`source`, you will first have to compile OpenMS
successfully on your platform of choice and then follow the `building from
source <../community/build_from_source.html>`_ instructions. Note that this may be
non-trivial and *is not recommended* for most users.

Wrap Classes (advanced users)
-----------------------------

In order to wrap new classes in pyOpenMS, read the following `guide
<../community/wrapping_workflows_new_classes.html>`_.
