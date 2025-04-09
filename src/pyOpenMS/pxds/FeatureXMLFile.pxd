from Feature cimport *
from FeatureMap cimport *
from String cimport *
from FeatureFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/FeatureXMLFile.h>" namespace "OpenMS":

    cdef cppclass FeatureXMLFile:
        FeatureXMLFile() except + nogil  # wrap-doc:This class provides Input/Output functionality for feature maps

        void load(String, FeatureMap &) except + nogil  # wrap-doc:Loads the file with name `filename` into `map` and calls updateRanges()
        void store(String, FeatureMap &) except + nogil  # wrap-doc:Stores the map `feature_map` in file with name `filename`

        FeatureFileOptions getOptions() except + nogil  # wrap-doc:Access to the options for loading/storing
        void setOptions(FeatureFileOptions) except + nogil  # wrap-doc:Setter for options for loading/storing

        Size loadSize(String path) except + nogil 

