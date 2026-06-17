Wrapping Workflows and New Classes
**********************************

How pyOpenMS Wraps C++ Classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pyOpenMS uses `nanobind <https://github.com/wjakob/nanobind>`_ to expose OpenMS
C++ classes to Python. The binding code is hand-maintained in C++ source files
under ``src/pyOpenMS/bindings/``.

The wrapping process works as follows:

- **Step 1:** The developer writes nanobind binding code in the appropriate
  ``bind_<domain>.cpp`` file under ``src/pyOpenMS/bindings/``. Each file
  corresponds to a domain of the OpenMS library (e.g., ``bind_kernel.cpp``
  for kernel classes, ``bind_format.cpp`` for file format classes).
- **Step 2:** CMake compiles each binding file into a separate Python extension
  module (e.g., ``_pyopenms_kernel.so``). All modules share type information
  via the ``NB_DOMAIN "pyopenms"`` mechanism.
- **Step 3:** The top-level ``pyopenms/__init__.py`` imports and re-exports all
  classes, so users simply write ``import pyopenms``.

In addition to the C++ bindings, pure Python convenience methods can be added
via the **addon system** in ``pyopenms/addons/``. These are injected into the
wrapped classes at import time.


Binding File Organization
^^^^^^^^^^^^^^^^^^^^^^^^^

Each binding file covers a domain of the OpenMS C++ library:

.. list-table::
   :header-rows: 1

   * - C++ Header Path
     - Binding File
   * - ``KERNEL/``
     - ``bind_kernel.cpp``
   * - ``KERNEL/MSSpectrum.h``
     - ``bind_spectrum.cpp``
   * - ``KERNEL/MSChromatogram.h``
     - ``bind_chromatogram.cpp``
   * - ``KERNEL/MSExperiment.h``
     - ``bind_experiment.cpp``
   * - ``FORMAT/``
     - ``bind_format.cpp``
   * - ``ANALYSIS/``
     - ``bind_analysis.cpp``
   * - ``CHEMISTRY/``
     - ``bind_chemistry.cpp``
   * - ``METADATA/``
     - ``bind_metadata.cpp``
   * - ``PROCESSING/``
     - ``bind_processing.cpp``
   * - ``FEATUREFINDER/``
     - ``bind_featurefinder.cpp``
   * - ``DATASTRUCTURES/``, ``MATH/``, ``CONCEPT/``
     - ``bind_datastructures.cpp``
   * - ``ML/``
     - ``bind_ml.cpp``
   * - Everything else
     - ``bind_misc.cpp``

Each class in a binding file is preceded by a section comment for navigation:

.. code-block:: cpp

   // --- ClassName ---


Maintaining Existing Wrappers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the C++ API of a wrapped class changes (e.g., a method is renamed, removed,
or has its signature changed), the corresponding binding code in
``src/pyOpenMS/bindings/bind_<domain>.cpp`` must be updated to match.


How to Wrap New Methods in Existing Classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To expose a new C++ method to Python, find the class in the appropriate
``bind_<domain>.cpp`` file and add a ``.def()`` call.

For a simple method:

.. code-block:: cpp

   .def("methodName", &OpenMS::ClassName::methodName,
        "Short description of the method")

For methods that need argument adaptation (e.g., filling output parameters by
reference), use a lambda wrapper:

.. code-block:: cpp

   .def("getValues", [](const OpenMS::ClassName& self) {
       std::vector<double> result;
       self.getValues(result);
       return result;
   }, "Returns the values as a list")

Use named arguments for clarity:

.. code-block:: cpp

   .def("setValue", &OpenMS::ClassName::setValue,
        "value"_a, "Sets the value")


How to Wrap New Classes
^^^^^^^^^^^^^^^^^^^^^^^

A Simple Example
----------------

To wrap a new OpenMS class:

1. **Choose the binding file** based on the C++ header location (see table above).
2. **Add the include** for the C++ header.
3. **Add the class binding** with constructors, copy support, and methods.

.. code-block:: cpp

   #include <OpenMS/MODULE/MyClass.h>

   // --- MyClass ---
   nb::class_<OpenMS::MyClass>(m, "MyClass", "Short class description")
       .def(nb::init<>())
       .def(nb::init<const OpenMS::MyClass &>())
       .def("__copy__", [](const OpenMS::MyClass& self) {
           return OpenMS::MyClass(self);
       })
       .def("__deepcopy__", [](const OpenMS::MyClass& self, nb::dict) {
           return OpenMS::MyClass(self);
       }, "memo"_a)
       .def("getValue", &OpenMS::MyClass::getValue,
            "Gets value (between 0 and 5)")
       .def("setValue", &OpenMS::MyClass::setValue,
            "v"_a, "Sets value (between 0 and 5)")
       ;

Key points:

- Always provide both a default constructor and a copy constructor.
- Always add ``__copy__`` and ``__deepcopy__`` methods.
- Use ``"param"_a`` for named parameters.
- Add a short docstring to each method. For longer docstrings, use raw strings:

  .. code-block:: cpp

     .def("myMethod", ..., R"doc(
         Longer description of the method.

         :param x: Description of x
         :returns: Description of return value
     )doc")


Inheritance
-----------

**Simple inheritance (no virtual destructor mismatch):**

If neither the class nor its base introduces a virtual destructor mismatch,
declare the base class directly:

.. code-block:: cpp

   nb::class_<OpenMS::Acquisition, OpenMS::MetaInfoInterface>(
       m, "Acquisition", "An acquisition")
       .def(nb::init<>())
       // ...
       ;

**Virtual destructor mismatch (critical):**

If the derived class has a virtual destructor (``virtual ~ClassName()`` or
``~ClassName() override``) but the base class (e.g., ``MetaInfoInterface``)
does **not**, you must **not** declare the base in nanobind. Instead, use
the helper templates from ``binding_utils.h``:

.. code-block:: cpp

   // WRONG - will segfault!
   nb::class_<OpenMS::PeptideHit, OpenMS::MetaInfoInterface>(m, "PeptideHit", ...)

   // CORRECT - use helper template
   auto peptidehit = nb::class_<OpenMS::PeptideHit>(m, "PeptideHit", "A peptide hit")
       .def(nb::init<>())
       // ... other methods ...
       ;
   def_MetaInfoInterface<OpenMS::PeptideHit>(peptidehit);

This pattern is required because nanobind computes incorrect pointer offsets
when a derived class introduces a vtable pointer that the base class lacks.

**How to check:** Look for ``virtual ~ClassName()`` or ``~ClassName() override``
in the C++ header file. If present and the base class has no virtual destructor,
use the helper template pattern.

Helper Templates
----------------

Reusable helper templates in ``binding_utils.h`` bind common OpenMS interfaces:

- ``def_MetaInfoInterface<T>(cls)`` -- binds ``getMetaValue``, ``setMetaValue``, ``isMetaEmpty``, etc.
- ``def_UniqueIdInterface<T>(cls)`` -- binds ``getUniqueId``, ``setUniqueId``, etc.
- ``def_CVTermList<T>(cls)`` -- binds ``addCVTerm``, ``getCVTerms``, ``hasCVTerm``, etc.
- ``def_DefaultParamHandler<T>(cls)`` -- binds ``setParameters``, ``getParameters``, ``getDefaults``, etc.
- ``def_ProgressLogger<T>(cls)`` -- binds ``setLogType``, ``startProgress``, ``endProgress``, etc.
- ``def_DocumentIdentifier<T>(cls)`` -- binds ``setIdentifier``, ``getLoadedFilePath``, etc.

Usage:

.. code-block:: cpp

   auto cls = nb::class_<OpenMS::MyAlgorithm>(m, "MyAlgorithm", "...")
       .def(nb::init<>())
       // ... class-specific methods ...
       ;
   def_DefaultParamHandler<OpenMS::MyAlgorithm>(cls);
   def_ProgressLogger<OpenMS::MyAlgorithm>(cls);


Enums
-----

Basic enum:

.. code-block:: cpp

   nb::enum_<OpenMS::MyEnum>(m, "MyEnum")
       .value("VALUE1", OpenMS::MyEnum::VALUE1)
       .value("VALUE2", OpenMS::MyEnum::VALUE2)
       ;

Enum supporting ``int()`` conversion (add ``nb::is_arithmetic()``):

.. code-block:: cpp

   nb::enum_<OpenMS::DriftTimeUnit>(m, "DriftTimeUnit", nb::is_arithmetic())
       .value("NONE", OpenMS::DriftTimeUnit::NONE)
       .value("MILLISECOND", OpenMS::DriftTimeUnit::MILLISECOND)
       ;

Nested enum (scoped under a class):

.. code-block:: cpp

   nb::enum_<OpenMS::MyClass::InnerEnum>(myclass, "InnerEnum")
       .value("A", OpenMS::MyClass::InnerEnum::A)
       ;


Operator Overloading
--------------------

Common operators can be bound using nanobind's operator support:

.. code-block:: cpp

   .def(nb::self == nb::self)
   .def(nb::self != nb::self)
   .def(nb::self < nb::self)
   .def(nb::self <= nb::self)
   .def(nb::self > nb::self)
   .def(nb::self >= nb::self)


Methods Returning NumPy Arrays
------------------------------

For C++ methods that fill output vectors by reference, wrap them with a lambda
that returns the data as numpy arrays:

.. code-block:: cpp

   .def("getData", [](const OpenMS::MyClass& self) {
       std::vector<double> mz, intensity;
       self.getData(mz, intensity);
       // Return as numpy arrays using nanobind's ndarray support
       // ...
   }, "Returns (mz, intensity) arrays")


Adding Pure Python Methods via Addons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For convenience methods that don't need C++ performance, use the addon system
in ``pyopenms/addons/``. Addon methods are pure Python functions injected into
wrapped classes at import time.

Create a file ``pyopenms/addons/classname.py`` (snake_case filename):

.. code-block:: python

   from __future__ import annotations
   from . import addon

   @addon("ClassName")
   def to_tuple(self) -> tuple:
       """Return as a (mz, intensity) tuple."""
       return (self.getMZ(), self.getIntensity())

The ``@addon("ClassName")`` decorator (with PascalCase class name) attaches the
function as a method on the specified class.

Guidelines for addons:

- Keep addons minimal -- only for non-performance-critical convenience methods.
- Performance-critical methods should be C++ lambdas in the binding files.
- Each addon file should document its methods with standard Python docstrings.


Build and Test
^^^^^^^^^^^^^^

.. code-block:: bash

   # Build pyOpenMS
   cmake --build OpenMS-build --target pyopenms -j$(nproc)

   # Run tests (from /tmp to avoid import shadowing with source pyopenms/)
   cd /tmp && PYTHONPATH=/path/to/OpenMS-build/pyOpenMS \
       python3 -m pytest /path/to/src/pyOpenMS/tests/ -v

   # Run a specific test file
   PYTHONPATH=OpenMS-build/pyOpenMS \
       python3 -m pytest src/pyOpenMS/tests/unittests/test_MyClass.py -v
