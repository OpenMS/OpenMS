from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from DRange cimport *

cdef extern from "<OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>" namespace "OpenMS":

    cdef cppclass FeatureFileOptions:

        FeatureFileOptions() nogil except +
        FeatureFileOptions(FeatureFileOptions) nogil except +

        void setMetadataOnly(bool) nogil except +
        bool getMetadataOnly()     nogil except +

        void setSizeOnly(bool) nogil except +
        bool getSizeOnly()     nogil except +

        void setLoadConvexHull(bool) nogil except +
        bool getLoadConvexHull()     nogil except +

        void setLoadSubordinates(bool) nogil except +
        bool getLoadSubordinates()     nogil except +

        void setRTRange(DRange1 & range_) nogil except +
        bool hasRTRange() nogil except +
        DRange1 getRTRange() nogil except +
        void setMZRange(DRange1 & range_) nogil except +
        bool hasMZRange() nogil except +
        DRange1 getMZRange() nogil except +
        void setIntensityRange(DRange1 & range_) nogil except +
        bool hasIntensityRange() nogil except +
        DRange1 getIntensityRange() nogil except +

