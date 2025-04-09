from libcpp.vector cimport vector as libcpp_vector
from ParamValue cimport *
from String cimport *
from StringList cimport *
from IntList cimport *
from DoubleList cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/DataValue.h>" namespace "OpenMS":

    cdef cppclass DataValue:
         DataValue() except + nogil 
         DataValue(DataValue &) except + nogil 
         DataValue(char *) except + nogil 
         DataValue(const String&) except + nogil 
         DataValue(int) except + nogil 
         DataValue(double) except + nogil 
         DataValue(StringList) except + nogil 
         DataValue(IntList) except + nogil 
         DataValue(DoubleList) except + nogil 
         DataValue(ParamValue) except + nogil 

         #conversion ops, different declarations as in c++ !
         int operator()(DataValue) except + nogil  #wrap-cast:toInt
         String operator()(DataValue) except + nogil  #wrap-cast:toString
         double operator()(DataValue) except + nogil  #wrap-cast:toDouble
         StringList toStringList() except + nogil 
         libcpp_vector[ double ] toDoubleList() except + nogil 
         libcpp_vector[ int ] toIntList() except + nogil 
         String toString() except + nogil 
         bool toBool() except + nogil 

         DataType valueType() except + nogil 

         int isEmpty() except + nogil 

         UnitType getUnitType() except + nogil 
         void setUnitType(UnitType u) except + nogil 

         bool hasUnit() except + nogil 
         int getUnit() except + nogil 
         void setUnit(int unit_id) except + nogil 

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

