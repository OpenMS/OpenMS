from Types cimport *
from String cimport *
from Feature cimport *
from FeatureMap cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FORMAT/MsInspectFile.h>" namespace "OpenMS":
    
    cdef cppclass MsInspectFile "OpenMS::MsInspectFile":
        MsInspectFile() nogil except +
        MsInspectFile(MsInspectFile &) nogil except + # compiler
        void load(const String & filename, FeatureMap & feature_map) nogil except +
            # wrap-doc:
                #   Loads a MsInspect file into a featureXML
                #   -----
                #   The content of the file is stored in `features`
                #   :raises:
                #     Exception: FileNotFound is thrown if the file could not be opened
                #   :raises:
                #     Exception: ParseError is thrown if an error occurs during parsing

        void store(const String & filename, MSSpectrum & spectrum) nogil except + # wrap-doc:Stores a featureXML as a MsInspect file

