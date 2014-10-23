from Types cimport *
from Feature cimport *
from Peak2D cimport *
from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/FeatureHandle.h>" namespace "OpenMS":
    
    cdef cppclass FeatureHandle(Peak2D,UniqueIdInterface) :
        # wrap-inherits:
        #  Peak2D
        #  UniqueIdInterface

        FeatureHandle() nogil except +
        FeatureHandle(FeatureHandle) nogil except +

        FeatureHandle(UInt64 map_index, Peak2D & point, UInt64 element_index) nogil except +
        # FeatureHandle(UInt64 map_index, BaseFeature & feature) nogil except +
        # FeatureHandleMutable_  asMutable() nogil except +

        UInt64 getMapIndex() nogil except +
        void setMapIndex(UInt64 i) nogil except +
        void setCharge(Int charge) nogil except +
        Int getCharge() nogil except +
        void setWidth(float width) nogil except +
        float getWidth() nogil except +

        bool operator==(FeatureHandle & i) nogil except +
        bool operator!=(FeatureHandle & i) nogil except +

