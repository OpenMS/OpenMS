Wrapping Workflows and New Classes
**********************************

How pyOpenMS Wraps C++ Classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

General concept of how the wrapping is done (all files are in ``src/pyOpenMS/``): 

- Step 1: The author declares which C++ classes and which functions of these
  classes she/he wants to wrap (expose to Python). This is done by writing the
  class and function declaration in a ``.pxd`` file in the ``pxds/`` folder.
- Step 2: The Python tool "autowrap" (developed for this project) creates the
  wrapping code automatically from the function declaration - see
  https://github.com/OpenMS/autowrap for an explanation of the autowrap
  tool. 
  Since not all code can be wrapped automatically using the pxd file, manual 
  Cython code sometimes need to be written in the ``addons/`` folder.
  Autowrap will create an output file at ``pyopenms/pyopenms.pyx`` 
  which can be interpreted by Cython.
- Step 3: Cython "cythonizes", i.e., translates the ``pyopenms/pyopenms.pyx`` to C++ code at
  ``pyopenms/pyopenms.cpp`` that contains everything to build a Python module.
- Step 4: A compiler, nowadays driven by CMake in OpenMS, compiles the C++ code to a Python
  extension module (pyopenms.so/dylib/dll/pyd) which is then importable in Python with
  ``import pyopenms``

Maintaining existing wrappers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the C++ API is changed, then pyOpenMS will
not build any more.  Thus, find the corresponding file in the ``pyOpenMS/pxds/``
folder and adjust the function declaration accordingly.

How to Wrap New Methods in Existing Classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lets say you have written a new method for an existing OpenMS class and you
would like to expose this method to pyOpenMS. First, identify the correct
``.pxd`` file in the ``src/pyOpenMS/pxds`` folder (for example for
:py:class:`~.Adduct` that would be `Adduct.pxd
<https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/pxds/Adduct.pxd>`_).
Open it and add your new function *with the correct indentation*:

- Place the full function declaration into the file (indented as the other functions)
- Check whether you are using any classes that are not yet imported, if so add a corresponding ``cimport`` statement to the top of the file. E.g., if your method is using :py:class:`~.MSExperiment`, then add ``from MSExperiment cimport *`` to the top (note it's ``cimport``, not ``import``).
- Most of the time it is clearer for a python user if you replace const-references with values in the function signature because the python function will always make a copy in
  that case. You should add `nogil except +` to the end of the signature to indicate that it can be multithreaded and throw exceptions.
  - Ex: ``void setType(Int a);`` becomes ``void setType(Int a) nogil except +`` 
  - Ex: ``const T& getType() const;`` becomes ``T getType() nogil except +``

- Remove any qualifiers (e.g. `const`) from the argument signatures, but leave reference and pointer indicators

  - Ex: ``const T&`` becomes ``T``, preventing an additional copy operation
  - Ex: ``T&`` will stay ``T&`` (indicating ``T`` needs to copied back to Python, make a note in the function docs that this argument may be changed)
  - Ex: ``T*`` will stay ``T*`` (indicating ``T`` needs to copied back to Python, make a note in the function docs that this argument may be changed)
  - One exception is ``OpenMS::String``. You can leave ``const String&`` as-is, since we have special handling for this.
  
- STL constructs are replaced with Cython constructs: ``std::vector<X>`` becomes ``libcpp_vector[ X ]`` etc. 
- Most complex STL constructs can be wrapped even if they are nested, however mixing them with user-defined types does not always work, see `Limitations <#Limitations>`_ below. Nested ``std::vector`` constructs work well even with user-defined (OpenMS-defined) types. However, ``std::map<String, X>`` does not work (since ``String`` is user-defined, however a primitive C++ type such as ``std::map<std::string, X>`` would work).
- Python cannot pass primitive data types by reference (therefore no ``int& res1``)
- Python does not allow default values for parameters in C(++) functions. Therefore you either accept the fact that the "optional" value always needs to be specified in python
  or you wrap it twice, once with and once without the optional parameter.
- Replace ``boost::shared_ptr<X>`` with ``shared_ptr[X]`` and add ``from smart_ptr cimport shared_ptr`` to the top
- Public members are simply added with ``Type member_name``. Private ones should be left out.
- You can inject documentation that will be shown when calling ``help()`` in the function by adding ``wrap-doc:Your documentation`` as a comment after the function:

  - Ex: ``void modifyWidget() nogil except + # wrap-doc:This changes your widget``
  - Warning: For a single-line comment, there should not be a space between wrap-doc and the following comment.
  - Note: The space between the hash and wrap-doc (# wrap-doc) is not necessary, but used for consistency.
  - Note: Please start the comment with a capital letter.
  
See the next section for a `simple example <#a-simple-example>`_ and a more `advanced example <#a-further-example>`_ of a wrapped class with several functions.

How to Wrap New Classes
^^^^^^^^^^^^^^^^^^^^^^^

A Simple Example
----------------

To wrap a new OpenMS class: Create a new ".pxd" file in the folder ``./pxds``. As
a small example, look at the `Adduct.pxd 
<https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/pxds/Adduct.pxd>`_ 
to get you started. Start with the following structure:

.. code-block:: bash

  from xxx cimport *
  cdef extern from "<OpenMS/path/to/header/Classname.h>" namespace "OpenMS":

      cdef cppclass ClassName(DefaultParamHandler):
          # wrap-inherits:
          #    DefaultParamHandler

          # default constructor
          ClassName() nogil except +
          # copy constructor
          ClassName(ClassName &) nogil except +

          Int getValue() nogil except + # wrap-doc:Gets value (between 0 and 5)
          void setValue(Int v) nogil except + # wrap-doc:Sets value (between 0 and 5)


- make sure to use ``ClassName:`` instead of ``ClassName(DefaultParamHandler):`` to
  wrap a class that does not inherit from another class (e.g., DefaultParamHandler)
  and also remove the two comments regarding inheritance below that line.
- always use ``cimport`` and not Python ``import``
- always add default constructor AND copy constructor to the code (note that the C++
  compiler will add a default copy constructor to any class even if not 
  specified in the header)
  This will be wrapped as __copy__ and __deepcopy__ method in python:

  .. code-block:: bash

    def __copy__(self):
     cdef ClassName rv = ClassName.__new__(ClassName)
     rv.inst = shared_ptr[ClassName](new ClassName(deref(self.inst.get())))
     return rv

    def __deepcopy__(self, memo):
       cdef ClassName rv = ClassName.__new__(ClassName)
       rv.inst = shared_ptr[ClassName](new ClassName(deref(self.inst.get())))
       return rv

- Remember to include a copy constructor (even if none was declared in the C++
  header file) since Cython will need it for certain operations. Otherwise you
  might see error messages like ``item2.inst = shared_ptr[_ClassName](new _ClassName(deref(it_terms))) Call with wrong number of arguments``.
- to expose a function to Python, copy the signature to your pxd file, e.g.
  ``DataValue getValue()`` and make sure you ``cimport`` all corresponding classes.
  Replace ``std::vector`` with the corresponding Cython vector, in this case
  ``libcpp_vector`` (see for example `PepXMLFile.pxd
  <https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/pxds/PepXMLFile.pxd>`_)
- you can add documentation that will show up in the interactive Python documentation (using ``help()``) using the ``wrap-doc`` qualifier

A Further Example
-----------------

A slightly more complicated class could look like this, where we demonstrate
how to handle a templated class with template ``T`` and static methods:

.. code-block:: cython

  from xxx cimport *
  from AbstractBaseClass cimport *
  from AbstractBaseClassImpl1 cimport *
  from AbstractBaseClassImpl2 cimport *
  cdef extern from "<OpenMS/path/to/header/Classname.h>" namespace "OpenMS":

      cdef cppclass ClassName[T](DefaultParamHandler):
          # wrap-inherits:
          #    DefaultParamHandler
          # 
          # wrap-instances:
          #   ClassName := ClassName[X]
          #   ClassNameY := ClassName[Y]

          ClassName() nogil except +
          ClassName(ClassName[T] &) nogil except + # wrap-ignore

          void method_name(int param1, double param2) nogil except +
          T method_returns_template_param() nogil except +

          size_t size() nogil except +
          T operator[](int) nogil except + # wrap-upper-limit:size()

          libcpp_vector[T].iterator begin() nogil except +  # wrap-iter-begin:__iter__(T)
          libcpp_vector[T].iterator end()   nogil except +  # wrap-iter-end:__iter__(T)

          void getWidgets(libcpp_vector[String] & keys) nogil except +
          void getWidgets(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getWAsInt

          # C++ signature: void process(AbstractBaseClass * widget)
          void process(AbstractBaseClassImpl1 * widget) nogil except +
          void process(AbstractBaseClassImpl2 * widget) nogil except +

  cdef extern from "<OpenMS/path/to/header/Classname.h>" namespace "OpenMS::Classname<OpenMS::X>":

      void static_method_name(int param1, double param2) nogil except + # wrap-attach:ClassName

  cdef extern from "<OpenMS/path/to/header/Classname.h>" namespace "OpenMS::Classname<OpenMS::Y>":

      void static_method_name(int param1, double param2) nogil except + # wrap-attach:ClassNameY


Here the copy constructor will not be wrapped but the Cython parser will import
it from C++ so that is is present (using ``wrap-ignore``). The ``operator[]``
will return an object of type ``X`` or ``Y`` depending on the template
argument ``T`` and contain a guard that the number may not be exceed ``size()``.

The wrapping of iterators allows for iteration over the objects inside the
``Classname`` container using the appropriate Python function (here
``__iter__`` with the indicated return type ``T``).

The ``wrap-as`` keyword allows the Python function to assume a different
name. 

Note that pointers to abstract base classes can be passed as arguments but the
classes have to be known at compile time, e.g. the function ``process``
takes a pointer to ``AbstractBaseClass`` which has two known
implementations ``AbstractBaseClassImpl1`` and
``AbstractBaseClassImpl2``. Then, the function needs to declared and
overloaded with both implementations as arguments as shown above.

An Example with Handwritten Addon Code
--------------------------------------

A more complex examples requires some hand-written wrapper code
(``pxds/Classname.pxd``), for example for singletons that implement a ``getInstance()``
method that returns a pointer to the singleton resource. Note that in this case
it is quite important to not let autowrap take over the pointer and possibly
delete it when the lifetime of the Python object ends. This is done through
``wrap-manual-memory`` and failing to doing so could lead to segmentation
faults in the program.

.. code-block:: cython

  from xxx cimport *
  cdef extern from "<OpenMS/path/to/header/Classname.h>" namespace "OpenMS":

      cdef cppclass ModificationsDB "OpenMS::ModificationsDB":
          # wrap-manual-memory
          # wrap-hash:
          #   getFullId().c_str()

          ClassName(ClassName[T] &) nogil except + # wrap-ignore

          void method_name(int param1, double param2) nogil except +

          int process(libcpp_vector[Peak1D].iterator, libcpp_vector[Peak1D].iterator) nogil except + # wrap-ignore

  cdef extern from "<OpenMS/path/to/header/Classname.h>" namespace "OpenMS::Classname":

      const ClassName* getInstance() nogil except + # wrap-ignore


Here the ``wrap-manual-memory`` keyword indicates that memory management
will be handled manually and autowrap can assume that a member called
``inst`` will be provided which implements a ``gets()`` method to
obtain a pointer to an object of C++ type ``Classname``.

We then have to provide such an object (``addons/Classname.pyx``):

.. code-block:: cython

    # This will go into the header (no empty lines below is *required*)
    # NOTE: _Classname is the C++ class while Classname is the Python class
    from Classname cimport Classname as _Classname
    cdef class ClassnameWrapper:
        # A small utility class holding a ptr and implementing get()
        cdef const _Classname* wrapped
        cdef setptr(self, const _Classname* wrapped): self.wrapped = wrapped
        cdef const _Classname* get(self) except *: return self.wrapped

        # This will go into the class (after the first empty line)
        # NOTE: we use 4 spaces indent
        # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
        cdef ClassnameWrapper inst

        def __init__(self):
          self.inst = ClassnameWrapper()
          # the following require some knowledge of the internals of autowrap:
          # we call the getInstance method to obtain raw ptr
          self.inst.setptr(_getInstance_Classname())

        def __dealloc__(self):
          # Careful here, the wrapped ptr is a single instance and we should not
          # reset it (which is why we used 'wrap-manual-dealloc')
          pass

        def process(self, Container c):
          # An example function here (processing Container c):
          return self.inst.get().process(c.inst.get().begin(), c.inst.get().end())


Note how the manual wrapping of the process functions allows us to
access the ``inst`` pointer of the argument as well as of the object
itself, allowing us to call C++ functions on both pointers. This makes it easy
to generate the required iterators and process the container efficiently.


.. _Limitations Section:

Considerations and Limitations
------------------------------

Further considerations and limitations:

- Inheritance: there are some limitations, see for example ``Precursor.pxd``

- Reference: arguments by reference may be copied under some circumstances. For
  example, if they are in an array then not the original argument is handed
  back, so comparisons might fail. Also, simple Python types like int, float
  etc cannot be passed by reference.

- operator+=: see for example ``AASequence.iadd`` in ``AASequence.pxd``

- operator==, !=, <=, <, >=, > are wrapped automatically

- Iterators: some limitations apply, see MSExperiment.pxd for an example

- copy-constructor becomes __copy__/__deepcopy__ in Python

- shared pointers: is handled automatically, check DataAccessHelper using ``shared_ptr[Spectrum]``. Use ``from smart_ptr cimport shared_ptr`` as import statement

These hints can be given to autowrap classes (also check the autowrap documentation):

- ``wrap-ignore`` is a hint for autowrap to not wrap the class (but the declaration might still be important for Cython to know about) 
- ``wrap-instances:`` for templated classes (see MSSpectrum.pxd)
- ``wrap-hash:`` hash function to use for ``__hash__`` (see Residue.pxd)
- ``wrap-manual-memory:`` hint that memory management will be done manually

These hints can be given to autowrap functions (also check the autowrap documentation):

- ``wrap-ignore`` is a hint for autowrap to not wrap the function (but the declaration might still be important for Cython to know about) 
- ``wrap-as:`` see for example AASequence
- ``wrap-iter-begin:``, ``wrap-iter-end:`` (see ConsensusMap.pxd)
- ``wrap-attach:`` enums, static methods (see for example VersionInfo.pxd)
- ``wrap-upper-limit:size()`` (see MSSpectrum.pxd)


Wrapping Code Yourself in ./addons
----------------------------------

Not all code can be wrapped automatically (yet). Place a file with the same (!)
name in the addons folder (e.g. ``myClass.pxd`` in ``pxds/`` and ``myClass.pyx`` in ``addons/``)
and leave two lines empty on the top (this is important). Start with 4 spaces
of indent and write your additional wrapper functions, adding a wrap-ignore
comment to the pxd file. See the example above, some additional examples, look into the ``src/pyOpenMS/addons/`` folder:

- `IDRipper.pyx <https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/addons/IDRipper.pyx>`_

  - for an example of both input and output of a complex STL construct (``map< String, pair<vector<>, vector<> >`` )

- `MSQuantifications.pyx <https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/addons/MSQuantifications.pyx>`_

  - for a ``vector< vector< pair <String,double > > >`` as input in registerExperiment
  - for a ``map< String, Ratio>`` in getRatios to get returned

- `QcMLFile.pyx <https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/addons/QcMLFile.pyx>`_
  - for a ``map< String, map< String,String> >`` as input

- `SequestInfile.pyx <https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/addons/SequestInfile.pyx>`_

  - for a ``map< String, vector<String> >`` to get returned

- `Attachment.pyx <https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/addons/Attachment.pyx>`_

  - for a ``vector< vector<String> >`` to get returned

- `ChromatogramExtractorAlgorithm.pxd <https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/pxds/ChromatogramExtractorAlgorithm.pxd>`_

  - for an example of an abstract base class (``ISpectrumAccess``) in the function
    ``extractChromatograms`` - this is solved by copy-pasting the function
    multiple times for each possible implementation of the abstract base class.

Make sure that you *always* declare your objects (all C++ and all Cython
objects need to be declared) using ``cdef`` Type name. Otherwise you get ``Cannot
convert ... to Python object`` errors.
