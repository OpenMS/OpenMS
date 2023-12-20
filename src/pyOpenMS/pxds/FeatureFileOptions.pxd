from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from DRange cimport *

cdef extern from "<OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>" namespace "OpenMS":

    cdef cppclass FeatureFileOptions:

        FeatureFileOptions() except + nogil  # wrap-doc:Options for loading files containing features.
        FeatureFileOptions(FeatureFileOptions &) except + nogil  # compiler

        void setMetadataOnly(bool) except + nogil  # wrap-doc:Sets whether or not to load only meta data
        bool getMetadataOnly()     except + nogil  # wrap-doc:Returns whether or not to load only meta data

        void setSizeOnly(bool) except + nogil  # wrap-doc:Sets whether or not to load only feature count
        bool getSizeOnly()     except + nogil  # wrap-doc:Returns whether or not to load only meta data

        void setLoadConvexHull(bool) except + nogil  # wrap-doc:Sets whether or not to load convex hull
        bool getLoadConvexHull()     except + nogil  # wrap-doc:Returns whether or not to load convex hull

        void setLoadSubordinates(bool) except + nogil  # wrap-doc:Sets whether or not load subordinates
        bool getLoadSubordinates()     except + nogil  # wrap-doc:Returns whether or not to load subordinates

        void setRTRange(DRange1 & range_) except + nogil  # wrap-doc:Restricts the range of RT values for peaks to load
        bool hasRTRange() except + nogil  # wrap-doc:Returns true if an RT range has been set
        DRange1 getRTRange() except + nogil  # wrap-doc:Returns the RT range
        void setMZRange(DRange1 & range_) except + nogil  # wrap-doc:Restricts the range of MZ values for peaks to load
        bool hasMZRange() except + nogil  # wrap-doc:Returns true if an MZ range has been set
        DRange1 getMZRange() except + nogil  # wrap-doc:Returns the MZ range
        void setIntensityRange(DRange1 & range_) except + nogil  # wrap-doc:Restricts the range of intensity values for peaks to load
        bool hasIntensityRange() except + nogil  # wrap-doc:Returns true if an intensity range has been set
        DRange1 getIntensityRange() except + nogil  # wrap-doc:Returns the intensity range

