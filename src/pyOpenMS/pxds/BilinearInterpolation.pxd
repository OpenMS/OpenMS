from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from Matrix cimport *

cdef extern from "<OpenMS/ML/INTERPOLATION/BilinearInterpolation.h>" namespace "OpenMS::Math":
    
    cdef cppclass BilinearInterpolation[KeyType,ValueType]:
        # wrap-instances:
        #  BilinearInterpolation := BilinearInterpolation[double, double]

        BilinearInterpolation() except + nogil 
        BilinearInterpolation(BilinearInterpolation &) except + nogil 
        ValueType value(KeyType arg_pos_0, KeyType arg_pos_1) except + nogil 
        void addValue(KeyType arg_pos_0, KeyType arg_pos_1, ValueType arg_value) except + nogil  # wrap-doc:Performs bilinear resampling.  The arg_value is split up and added to the data points around arg_pos.  ("forward resampling")

        Matrix[ValueType] getData() except + nogil 
        void setData(Matrix[ValueType] & data) except + nogil  # wrap-doc:Assigns data to the internal random access container storing the data. SourceContainer must be assignable to ContainerType

        bool empty() except + nogil 

        KeyType key2index_0(KeyType pos) except + nogil  # wrap-doc:The transformation from "outside" to "inside" coordinates
        KeyType index2key_0(KeyType pos) except + nogil  # wrap-doc:The transformation from "inside" to "outside" coordinates
        KeyType key2index_1(KeyType pos) except + nogil  # wrap-doc:The transformation from "outside" to "inside" coordinates
        KeyType index2key_1(KeyType pos) except + nogil  # wrap-doc:The transformation from "inside" to "outside" coordinates

        KeyType getScale_0() except + nogil 
        void setScale_0(KeyType & scale) except + nogil 
        KeyType getScale_1() except + nogil 
        void setScale_1(KeyType & scale) except + nogil 

        KeyType getOffset_0() except + nogil  # wrap-doc:Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)"
        void setOffset_0(KeyType & offset) except + nogil  
        KeyType getOffset_1() except + nogil  # wrap-doc:Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)"
        void setOffset_1(KeyType & offset) except + nogil 

        void setMapping_0(KeyType & scale, KeyType & inside, KeyType & outside) except + nogil 
        void setMapping_0(KeyType & inside_low, KeyType & outside_low, KeyType & inside_high, KeyType & outside_high) except + nogil 
        void setMapping_1(KeyType & scale, KeyType & inside, KeyType & outside) except + nogil 
        void setMapping_1(KeyType & inside_low, KeyType & outside_low, KeyType & inside_high, KeyType & outside_high) except + nogil 

        KeyType  getInsideReferencePoint_0() except + nogil 
        KeyType  getInsideReferencePoint_1() except + nogil 

        KeyType  getOutsideReferencePoint_0() except + nogil 
        KeyType  getOutsideReferencePoint_1() except + nogil 

        KeyType supportMin_0() except + nogil  # wrap-doc:Lower boundary of the support, in "outside" coordinates
        KeyType supportMin_1() except + nogil  # wrap-doc:Lower boundary of the support, in "outside" coordinates

        KeyType supportMax_0() except + nogil  # wrap-doc:Upper boundary of the support, in "outside" coordinates
        KeyType supportMax_1() except + nogil  # wrap-doc:Upper boundary of the support, in "outside" coordinates

        # BilinearInterpolation(KeyType scale, KeyType offset) except + nogil 
