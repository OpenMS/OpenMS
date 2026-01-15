from Feature cimport *
from FeatureMap cimport *
from String cimport *
from FeatureFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/FeatureXMLFile.h>" namespace "OpenMS":

    cdef cppclass FeatureXMLFile:
        # wrap-doc:
        #  File adapter for featureXML files
        #  
        #  Provides methods to load and store feature maps in featureXML format.
        #  FeatureXML files store LC-MS features with their convex hulls, intensities, and metadata.
        #  
        #  Usage:
        #  
        #  .. code-block:: python
        #  
        #    fm = FeatureMap()
        #    FeatureXMLFile().load("test.featureXML", fm)
        #    for feature in fm:
        #      print(feature.getRT(), feature.getMZ(), feature.getIntensity())

        FeatureXMLFile() except + nogil

        void load(String, FeatureMap &) except + nogil  # wrap-doc:Loads the file with name `filename` into `map` and calls updateRanges()
        void store(String, FeatureMap &) except + nogil  # wrap-doc:Stores the map `feature_map` in file with name `filename`

        FeatureFileOptions getOptions() except + nogil  # wrap-doc:Access to the options for loading/storing
        void setOptions(FeatureFileOptions) except + nogil  # wrap-doc:Setter for options for loading/storing

        Size loadSize(String path) except + nogil  # wrap-doc:Counts the number of features in the file without loading the full data

