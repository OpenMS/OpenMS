from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from StringList cimport *
from IntList cimport *
from DoubleList cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DataValue.h>" namespace "OpenMS":

    cdef cppclass DataValue:
         DataValue() nogil except +
         DataValue(DataValue) nogil except + # wrap-ignore
         DataValue(char *) nogil except +
         DataValue(int)    nogil except +
         DataValue(double)    nogil except +
         DataValue(StringList)     nogil except +
         DataValue(IntList)  nogil except +
         DataValue(DoubleList)  nogil except +

         #conversion ops, different declarations as in c++ !
         int operator()(int)    nogil except +  #wrap-cast:toInt
         libcpp_string operator()(DataValue)  nogil except +  #wrap-cast:toString
         double operator()(DataValue)  nogil except +         #wrap-cast:toDouble
         StringList toStringList() nogil except +
         libcpp_vector[ double ] toDoubleList() nogil except +
         libcpp_vector[ int ] toIntList() nogil except +

         DataType valueType() nogil except +
         int isEmpty() nogil except +

         # QString toQString()
         String toString() nogil except +
         bool toBool() nogil except +
         bool hasUnit() nogil except +
         String  getUnit() nogil except +
         void setUnit(const String & unit) nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/DataValue.h>" namespace "OpenMS::DataValue":

    cdef enum DataType:
         STRING_VALUE, INT_VALUE, DOUBLE_VALUE, STRING_LIST, INT_LIST, \
         DOUBLE_LIST, EMPTY_VALUE
