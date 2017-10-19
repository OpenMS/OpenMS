from Types cimport *
from KDTreeFeatureMaps cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>" namespace "OpenMS":
    
    cdef cppclass KDTreeFeatureNode "OpenMS::KDTreeFeatureNode":
        # KDTreeFeatureNode() nogil except +
        KDTreeFeatureNode(KDTreeFeatureNode) nogil except +
        KDTreeFeatureNode(KDTreeFeatureMaps * data, Size idx) nogil except +
        double operator[](Size i) nogil except +
        Size getIndex() nogil except +

