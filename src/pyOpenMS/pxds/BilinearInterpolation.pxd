from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from Types cimport *
from Matrix cimport *

cdef extern from "<OpenMS/MATH/MISC/BilinearInterpolation.h>" namespace "OpenMS::Math":
    
    cdef cppclass BilinearInterpolation[KeyType,ValueType]:
        # wrap-instances:
        #   BilinearInterpolation := BilinearInterpolation[double, double]

        BilinearInterpolation() nogil except +
        BilinearInterpolation(BilinearInterpolation &) nogil except +
        ValueType value(KeyType arg_pos_0, KeyType arg_pos_1) nogil except +
        void addValue(KeyType arg_pos_0, KeyType arg_pos_1, ValueType arg_value) nogil except + # wrap-doc:Performs bilinear resampling.  The arg_value is split up and added to the data points around arg_pos.  ("forward resampling")

        Matrix[ValueType] getData() nogil except +
        void setData(Matrix[ValueType] & data) nogil except + # wrap-doc:Assigns data to the internal random access container storing the data. SourceContainer must be assignable to ContainerType

        bool empty() nogil except +

        KeyType key2index_0(KeyType pos) nogil except + # wrap-doc:The transformation from "outside" to "inside" coordinates
        KeyType index2key_0(KeyType pos) nogil except + # wrap-doc:The transformation from "inside" to "outside" coordinates
        KeyType key2index_1(KeyType pos) nogil except + # wrap-doc:The transformation from "outside" to "inside" coordinates
        KeyType index2key_1(KeyType pos) nogil except + # wrap-doc:The transformation from "inside" to "outside" coordinates

        KeyType getScale_0() nogil except +
        void setScale_0(KeyType & scale) nogil except +
        KeyType getScale_1() nogil except +
        void setScale_1(KeyType & scale) nogil except +

        KeyType getOffset_0() nogil except + # wrap-doc:Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)"
        void setOffset_0(KeyType & offset) nogil except + 
        KeyType getOffset_1() nogil except + # wrap-doc:Accessor.  "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)"
        void setOffset_1(KeyType & offset) nogil except +

        void setMapping_0(KeyType & scale, KeyType & inside, KeyType & outside) nogil except +
        void setMapping_0(KeyType & inside_low, KeyType & outside_low, KeyType & inside_high, KeyType & outside_high) nogil except +
        void setMapping_1(KeyType & scale, KeyType & inside, KeyType & outside) nogil except +
        void setMapping_1(KeyType & inside_low, KeyType & outside_low, KeyType & inside_high, KeyType & outside_high) nogil except +

        KeyType  getInsideReferencePoint_0() nogil except +
        KeyType  getInsideReferencePoint_1() nogil except +

        KeyType  getOutsideReferencePoint_0() nogil except +
        KeyType  getOutsideReferencePoint_1() nogil except +

        KeyType supportMin_0() nogil except + # wrap-doc:Lower boundary of the support, in "outside" coordinates
        KeyType supportMin_1() nogil except + # wrap-doc:Lower boundary of the support, in "outside" coordinates

        KeyType supportMax_0() nogil except + # wrap-doc:Upper boundary of the support, in "outside" coordinates
        KeyType supportMax_1() nogil except + # wrap-doc:Upper boundary of the support, in "outside" coordinates

        # BilinearInterpolation(KeyType scale, KeyType offset) nogil except +
