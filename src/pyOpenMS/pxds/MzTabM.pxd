from Types cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/FORMAT/MzTabM.h>" namespace "OpenMS":

    cdef cppclass MzTabM:
        # wrap-doc:
                #   Data model of MzTabM files
                #   -----
                #   Please see the official MzTabM specification at https://github.com/HUPO-PSI/mzTab/tree/master/specification_document-releases/2_0-Metabolomics-Release

        MzTabM() nogil except +
        MzTabM(MzTabM &) nogil except + # compiler

        MzTabM exportFeatureMapToMzTabM(FeatureMap feature_map) # wrap-doc:Export FeatureMap with Identifications to MzTabM
