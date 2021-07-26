from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from DRange cimport *

cdef extern from "<OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>" namespace "OpenMS":

    cdef cppclass FeatureFileOptions:

        FeatureFileOptions() nogil except + # wrap-doc:Options for loading files containing features.
        FeatureFileOptions(FeatureFileOptions &) nogil except + # compiler

        void setMetadataOnly(bool) nogil except + # wrap-doc:Sets whether or not to load only meta data
        bool getMetadataOnly()     nogil except + # wrap-doc:Returns whether or not to load only meta data

        void setSizeOnly(bool) nogil except + # wrap-doc:Sets whether or not to load only feature count
        bool getSizeOnly()     nogil except + # wrap-doc:Returns whether or not to load only meta data

        void setLoadConvexHull(bool) nogil except + # wrap-doc:Sets whether or not to load convex hull
        bool getLoadConvexHull()     nogil except + # wrap-doc:Returns whether or not to load convex hull

        void setLoadSubordinates(bool) nogil except + # wrap-doc:Sets whether or not load subordinates
        bool getLoadSubordinates()     nogil except + # wrap-doc:Returns whether or not to load subordinates

        void setRTRange(DRange1 & range_) nogil except + # wrap-doc:Restricts the range of RT values for peaks to load
        bool hasRTRange() nogil except + # wrap-doc:Returns true if an RT range has been set
        DRange1 getRTRange() nogil except + # wrap-doc:Returns the RT range
        void setMZRange(DRange1 & range_) nogil except + # wrap-doc:Restricts the range of MZ values for peaks to load
        bool hasMZRange() nogil except + # wrap-doc:Returns true if an MZ range has been set
        DRange1 getMZRange() nogil except + # wrap-doc:Returns the MZ range
        void setIntensityRange(DRange1 & range_) nogil except + # wrap-doc:Restricts the range of intensity values for peaks to load
        bool hasIntensityRange() nogil except + # wrap-doc:Returns true if an intensity range has been set
        DRange1 getIntensityRange() nogil except + # wrap-doc:Returns the intensity range

