from Types cimport *
from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from ParamValue cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":

    cdef cppclass ParamEntry:

        libcpp_string name
        libcpp_string description
        ParamValue value

        # The next two properties are tricky, see https://github.com/OpenMS/autowrap/issues/173
        # Cython's autoconversion will be used here -> it accepts str and bytes and returns bytes
        # Therefore we have a mismatch/narrowing https://github.com/python/mypy/issues/3004
        # For now we use the narrow type bytes -> Although we could accept str the typing will suggest not to.
        # but the output will be correctly typed.
        
        libcpp_set[libcpp_string] tags
        libcpp_vector[libcpp_string] valid_strings
        double  max_float
        double  min_float
        Int  max_int
        Int  min_int

        ParamEntry() except + nogil  # TODO
        ParamEntry(ParamEntry &) except + nogil 
        ParamEntry(libcpp_string n, ParamValue v, libcpp_string d, libcpp_vector[libcpp_string] t) except + nogil 
        ParamEntry(libcpp_string n, ParamValue v, libcpp_string d) except + nogil 

        # TODO: wrap bool isValid(libcpp_string & message) maybe as libcpp_pair[bool, libcpp_str] isValid() if possible
        bool operator==(ParamEntry) except + nogil 

