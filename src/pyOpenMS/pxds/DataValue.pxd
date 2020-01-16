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
         DataValue(const String&) nogil except +
         DataValue(int) nogil except +
         DataValue(double) nogil except +
         DataValue(StringList) nogil except +
         DataValue(IntList) nogil except +
         DataValue(DoubleList) nogil except +

         #conversion ops, different declarations as in c++ !
         int operator()(int) nogil except + #wrap-cast:toInt
         String operator()(DataValue) nogil except + #wrap-cast:toString
         double operator()(DataValue) nogil except + #wrap-cast:toDouble
         StringList toStringList() nogil except +
         libcpp_vector[ double ] toDoubleList() nogil except +
         libcpp_vector[ int ] toIntList() nogil except +
         String toString() nogil except +
         bool toBool() nogil except +

         DataType valueType() nogil except +

         int isEmpty() nogil except +

         UnitType getUnitType() nogil except +
         void setUnitType(UnitType u) nogil except +

         bool hasUnit() nogil except +
         int getUnit() nogil except +
         void setUnit(int unit_id) nogil except +

cdef extern from "<OpenMS/DATASTRUCTURES/DataValue.h>" namespace "OpenMS::DataValue":

    cdef enum DataType "OpenMS::DataValue::DataType":
        STRING_VALUE # string value
        INT_VALUE # integer value
        DOUBLE_VALUE # double value
        STRING_LIST # string list
        INT_LIST # integer list
        DOUBLE_LIST # double list
        EMPTY_VALUE # empty value

    cdef enum UnitType "OpenMS::DataValue::UnitType":
        UNIT_ONTOLOGY # unit.ontology UO:
        MS_ONTOLOGY # ms.ontology MS:
        OTHER # undefined ontology

