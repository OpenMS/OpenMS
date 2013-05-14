from libcpp.string cimport string as libcpp_string
from libcpp.vector cimport vector as libcpp_vector
from StringList cimport *
from String cimport *
from IntList cimport *
from DoubleList cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DataValue.h>" namespace "OpenMS":

    cdef cppclass DataValue:
         DataValue()
         DataValue(DataValue) nogil except + # wrap-ignore
         DataValue(libcpp_vector[libcpp_string]) nogil except +
         DataValue(libcpp_vector[int])  nogil except +
         DataValue(libcpp_vector[float])  nogil except +
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
         StringList operator()(DataValue) nogil except +      #wrap-cast:toStringList
         DoubleList operator()(DataValue) nogil except +      #wrap-cast:toDoubleList
         IntList operator()(DataValue)  nogil except +        #wrap-cast:toIntList

         DataType valueType() nogil except +
         int isEmpty() nogil except +

         # QString toQString()
         String toString()
         bool toBool()
         bool hasUnit()
         String  getUnit()
         void setUnit(String & unit)

cdef extern from "<OpenMS/DATASTRUCTURES/DataValue.h>" namespace "OpenMS::DataValue":

    cdef enum DataType:
         STRING_VALUE, INT_VALUE, DOUBLE_VALUE, STRING_LIST, INT_LIST, \
         DOUBLE_LIST, EMPTY_VALUE
