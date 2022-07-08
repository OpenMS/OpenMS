from DefaultParamHandler cimport *
from MSExperiment cimport *
from Param cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/MSPGenericFile.h>" namespace "OpenMS":
    cdef cppclass MSPGenericFile(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        MSPGenericFile() nogil except +
        MSPGenericFile(MSPGenericFile &) nogil except + # compiler
        MSPGenericFile(const String& filename, MSExperiment& library) nogil except +

        void load(const String& filename, MSExperiment& library) nogil except +
            # wrap-doc:
                #   Load the file's data and metadata, and save it into an `MSExperiment`
                #   -----
                #   :param filename: Path to the MSP input file
                #   :param library: The variable into which the extracted information will be saved
                #   :raises:
                #     Exception: FileNotFound If the file could not be found

        void store(const String& filename, const MSExperiment& library) nogil except +
            # wrap-doc:
                #   Save data and metadata into a file
                #   -----
                #   :param filename: Path to the MSP input file
                #   :param library: The variable from which extracted information will be saved
                #   :raises:
                #     Exception: FileNotWritable If the file is not writable

        void getDefaultParameters(Param & params) nogil except + # wrap-doc:Returns the class' default parameters

