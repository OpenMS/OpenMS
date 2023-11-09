from Types cimport *
from KDTreeFeatureMaps cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>" namespace "OpenMS":
    
    cdef cppclass KDTreeFeatureNode "OpenMS::KDTreeFeatureNode":
    # wrap-doc:
        # A node of the kD-tree with pointer to corresponding data and index
        
        # KDTreeFeatureNode() except + nogil 
        KDTreeFeatureNode(KDTreeFeatureNode &) except + nogil 
        KDTreeFeatureNode(KDTreeFeatureMaps * data, Size idx) except + nogil 
        double operator[](Size i) except + nogil 
        Size getIndex() except + nogil  # wrap-doc:Returns index of corresponding feature in data_

