from Types cimport *
from KDTreeFeatureMaps cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>" namespace "OpenMS":
    
    cdef cppclass KDTreeFeatureNode "OpenMS::KDTreeFeatureNode":
    # wrap-doc:
        # A node of the kD-tree with pointer to corresponding data and index
        
        # KDTreeFeatureNode() nogil except +
        KDTreeFeatureNode(KDTreeFeatureNode &) nogil except +
        KDTreeFeatureNode(KDTreeFeatureMaps * data, Size idx) nogil except +
        double operator[](Size i) nogil except +
        Size getIndex() nogil except + # wrap-doc:Returns index of corresponding feature in data_

