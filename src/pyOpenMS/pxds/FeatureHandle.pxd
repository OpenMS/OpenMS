from Types cimport *
from Feature cimport *
from Peak2D cimport *
from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/FeatureHandle.h>" namespace "OpenMS":
    
    cdef cppclass FeatureHandle(Peak2D,UniqueIdInterface) :
        # wrap-inherits:
        #  Peak2D
        #  UniqueIdInterface

        FeatureHandle() except + nogil  # wrap-doc:Representation of a Peak2D, RichPeak2D or Feature
        FeatureHandle(FeatureHandle &) except + nogil 

        FeatureHandle(UInt64 map_index, Peak2D & point, UInt64 element_index) except + nogil 
        # FeatureHandle(UInt64 map_index, BaseFeature & feature) except + nogil 
        # FeatureHandleMutable_  asMutable() except + nogil 

        UInt64 getMapIndex() except + nogil  # wrap-doc:Returns the map index
        void setMapIndex(UInt64 i) except + nogil  # wrap-doc:Sets the map index
        void setCharge(Int charge) except + nogil  # wrap-doc:Sets the charge
        Int getCharge() except + nogil  # wrap-doc:Returns the charge
        void setWidth(float width) except + nogil  # wrap-doc:Sets the width (FWHM)
        float getWidth() except + nogil  # wrap-doc:Returns the width (FWHM)

        bool operator==(FeatureHandle & i) except + nogil 
        bool operator!=(FeatureHandle & i) except + nogil 

