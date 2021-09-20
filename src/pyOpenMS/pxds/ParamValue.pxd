from libcpp cimport bool
from libcpp.string cimport string as libcpp_utf8_string # triggers input conversion provider for std string
from libcpp.string cimport string as libcpp_utf8_output_string #triggers output conversion provider for std string
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/DATASTRUCTURES/ParamValue.h>" namespace "OpenMS":

    cdef cppclass ParamValue:
        # wrap-doc:
            #   Class to hold strings, numeric values, vectors of strings and vectors of numeric values using the stl types
            #   -----
            #   - To choose one of these types, just use the appropriate constructor
            #   - Automatic conversion is supported and throws Exceptions in case of invalid conversions
            #   - An empty object is created with the default constructor

         ParamValue() nogil except +
         ParamValue(ParamValue &) nogil except +
         ParamValue(char *) nogil except +
         ParamValue(const libcpp_utf8_string&) nogil except +
         ParamValue(int) nogil except +
         ParamValue(double) nogil except +
         ParamValue(libcpp_vector[ libcpp_utf8_string ]) nogil except +
         ParamValue(libcpp_vector[ int ]) nogil except +
         ParamValue(libcpp_vector[ double ]) nogil except +

         #conversion ops, different declarations as in c++ !
         int operator()(ParamValue) nogil except + #wrap-cast:toInt
         libcpp_utf8_output_string operator()(ParamValue) nogil except + #wrap-cast:toString
         double operator()(ParamValue) nogil except + #wrap-cast:toDouble
         libcpp_vector[ libcpp_utf8_string ] toStringVector() nogil except + # wrap-doc:Explicitly convert ParamValue to string vector
         libcpp_vector[ double ] toDoubleVector() nogil except + # wrap-doc:Explicitly convert ParamValue to DoubleList
         libcpp_vector[ int ] toIntVector() nogil except + # wrap-doc:Explicitly convert ParamValue to IntList
         bool toBool() nogil except + # wrap-doc:Converts the strings 'true' and 'false' to a bool

         ValueType valueType() nogil except +

         int isEmpty() nogil except + # wrap-doc:Test if the value is empty

cdef extern from "<OpenMS/DATASTRUCTURES/ParamValue.h>" namespace "OpenMS::ParamValue":

    cdef enum ValueType "OpenMS::ParamValue::ValueType":
        STRING_VALUE # string value
        INT_VALUE # integer value
        DOUBLE_VALUE # double value
        STRING_LIST # string list
        INT_LIST # integer list
        DOUBLE_LIST # double list
        EMPTY_VALUE # empty value
