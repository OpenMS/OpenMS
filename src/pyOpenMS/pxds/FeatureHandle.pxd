from Types cimport *
from Feature cimport *
from Peak2D cimport *
from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/FeatureHandle.h>" namespace "OpenMS":
    
    cdef cppclass FeatureHandle(Peak2D,UniqueIdInterface) :
        # wrap-inherits:
        #  Peak2D
        #  UniqueIdInterface

        FeatureHandle() nogil except + # wrap-doc:Representation of a Peak2D, RichPeak2D or Feature
        FeatureHandle(FeatureHandle &) nogil except +

        FeatureHandle(UInt64 map_index, Peak2D & point, UInt64 element_index) nogil except +
        # FeatureHandle(UInt64 map_index, BaseFeature & feature) nogil except +
        # FeatureHandleMutable_  asMutable() nogil except +

        UInt64 getMapIndex() nogil except + # wrap-doc:Returns the map index
        void setMapIndex(UInt64 i) nogil except + # wrap-doc:Sets the map index
        void setCharge(Int charge) nogil except + # wrap-doc:Sets the charge
        Int getCharge() nogil except + # wrap-doc:Returns the charge
        void setWidth(float width) nogil except + # wrap-doc:Sets the width (FWHM)
        float getWidth() nogil except + # wrap-doc:Returns the width (FWHM)

        bool operator==(FeatureHandle & i) nogil except +
        bool operator!=(FeatureHandle & i) nogil except +

