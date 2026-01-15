from ConsensusMap cimport *
from String cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/ConsensusXMLFile.h>" namespace "OpenMS":

    cdef cppclass ConsensusXMLFile:
        # wrap-doc:
        #  File adapter for consensusXML files
        #
        #  Provides methods to load and store consensus maps in consensusXML format.
        #  ConsensusXML files store consensus features linking corresponding features across
        #  multiple LC-MS runs, typically used for quantification workflows.
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    cm = ConsensusMap()
        #    ConsensusXMLFile().load("test.consensusXML", cm)
        #    for cf in cm:
        #      print(cf.getRT(), cf.getMZ(), cf.getIntensity())

        ConsensusXMLFile() except + nogil

        void load(String, ConsensusMap &) except + nogil  # wrap-doc:Loads a consensus map from file and calls updateRanges
        void store(String, ConsensusMap &) except + nogil  # wrap-doc:Stores a consensus map to file

        PeakFileOptions getOptions() except + nogil  # wrap-doc:Mutable access to the options for loading/storing
