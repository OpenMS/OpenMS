from Types cimport *
from libcpp.pair cimport pair as libcpp_pair
from String cimport *
from MSExperiment cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FORMAT/MascotGenericFile.h>" namespace "OpenMS":
    
    cdef cppclass MascotGenericFile(ProgressLogger, DefaultParamHandler) :
        # wrap-inherits:
        #  ProgressLogger
        #  DefaultParamHandler
        MascotGenericFile() nogil except +
        MascotGenericFile(MascotGenericFile &) nogil except +
        void store(const String & filename, MSExperiment & experiment) nogil except +
        # NAMESPACE # void store(std::ostream & os, const String & filename, MSExperiment & experiment) nogil except +
        void load(const String & filename, MSExperiment & exp) nogil except +
            # wrap-doc:
                #   Loads a Mascot Generic File into a PeakMap
                #   -----
                #   :param filename: File name which the map should be read from
                #   :param exp: The map which is filled with the data from the given file
                #   :raises:
                #     Exception: FileNotFound is thrown if the given file could not be found

        libcpp_pair[ String, String ] getHTTPPeakListEnclosure(const String & filename) nogil except +
            # wrap-doc:
                #   Enclosing Strings of the peak list body for HTTP submission
                #   -----
                #   Can be used to embed custom content into HTTP submission (when writing only the MGF header in HTTP format and then
                #   adding the peaks (in whatever format, e.g. mzXML) enclosed in this body
                #   The `filename` can later be found in the Mascot response
      
        void updateMembers_() nogil except + # wrap-doc:Docu in base class

