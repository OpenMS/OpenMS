from libcpp cimport bool
from libcpp.string cimport string as libcpp_string
from libcpp.string cimport string as libcpp_utf8_string # triggers input conversion provider for std string
from libcpp.string cimport string as libcpp_utf8_output_string #triggers output conversion provider for std string
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/DATASTRUCTURES/ParamValue.h>" namespace "OpenMS":

    cdef cppclass ParamValue:
        # wrap-doc:
            #  Class to hold strings, numeric values, vectors of strings and vectors of numeric values using the stl types
            #  
            #  - To choose one of these types, just use the appropriate constructor
            #  - Automatic conversion is supported and throws Exceptions in case of invalid conversions
            #  - An empty object is created with the default constructor

         ParamValue() except + nogil 
         ParamValue(ParamValue &) except + nogil 
         ParamValue(char *) except + nogil 
         ParamValue(const libcpp_utf8_string&) except + nogil 
         ParamValue(int) except + nogil 
         ParamValue(double) except + nogil 
         ParamValue(libcpp_vector[ libcpp_utf8_string ]) except + nogil 
         ParamValue(libcpp_vector[ int ]) except + nogil 
         ParamValue(libcpp_vector[ double ]) except + nogil 

         #conversion ops, different declarations as in c++ !
         int operator()(ParamValue) except + nogil  #wrap-cast:toInt
         libcpp_utf8_output_string operator()(ParamValue) except + nogil  #wrap-cast:toString
         double operator()(ParamValue) except + nogil  #wrap-cast:toDouble
         libcpp_vector[ libcpp_string ] toStringVector() except + nogil  # wrap-doc:Explicitly convert ParamValue to string vector
         libcpp_vector[ double ] toDoubleVector() except + nogil  # wrap-doc:Explicitly convert ParamValue to DoubleList
         libcpp_vector[ int ] toIntVector() except + nogil  # wrap-doc:Explicitly convert ParamValue to IntList
         bool toBool() except + nogil  # wrap-doc:Converts the strings 'true' and 'false' to a bool

         ValueType valueType() except + nogil 

         int isEmpty() except + nogil  # wrap-doc:Test if the value is empty

cdef extern from "<OpenMS/DATASTRUCTURES/ParamValue.h>" namespace "OpenMS::ParamValue":

    cdef enum ValueType "OpenMS::ParamValue::ValueType":
        STRING_VALUE # string value
        INT_VALUE # integer value
        DOUBLE_VALUE # double value
        STRING_LIST # string list
        INT_LIST # integer list
        DOUBLE_LIST # double list
        EMPTY_VALUE # empty value
