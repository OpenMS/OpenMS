from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *

cdef extern from "<OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>" namespace "OpenMS":

    cdef cppclass FeatureFileOptions:

        FeatureFileOptions() nogil except +
        FeatureFileOptions(FeatureFileOptions) nogil except + # wrap-ignore
        
        void setMetadataOnly(bool) nogil except +
        bool getMetadataOnly()     nogil except +

        #void setSizeOnly(bool) nogil except +
        #bool getSizeOnly()     nogil except +

        void setLoadConvexHull(bool) nogil except +
        bool getLoadConvexHull()     nogil except +

        void setLoadSubordinates(bool) nogil except +
        bool getLoadSubordinates()     nogil except +

