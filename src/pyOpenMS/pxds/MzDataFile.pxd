
from MSExperiment  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from StringList cimport *
from ProgressLogger cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/MzDataFile.h>" namespace "OpenMS":

    cdef cppclass MzDataFile(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        MzDataFile() nogil except + # wrap-doc:File adapter for MzData files
        MzDataFile(MzDataFile &) nogil except +

        void load(const String& filename, MSExperiment & map) nogil except +
            # wrap-doc:
            #   Loads a map from a MzData file
            #   -----
            #   :param filename: Directory of the file with the file name
            #   :param map: It has to be a MSExperiment or have the same interface
            #   -----
            #   :raises:
            #     Exception: FileNotFound is thrown if the file could not be opened
            #   :raises:
            #     Exception: ParseError is thrown if an error occurs during parsing

        void store(const String& filename, MSExperiment & map) nogil except +
            # wrap-doc:
            #   Stores a map in a MzData file
            #   -----
            #   :param filename: Directory of the file with the file name
            #   :param map: It has to be a MSExperiment or have the same interface
            #   -----
            #   :raises:
            #     Exception: UnableToCreateFile is thrown if the file could not be created

        PeakFileOptions getOptions() nogil except + # wrap-doc:Returns the options for loading/storing
        void setOptions(PeakFileOptions) nogil except + # wrap-doc:Sets options for loading/storing

        bool isSemanticallyValid(const String& filename, StringList & errors, StringList & warnings) nogil except +
            # wrap-doc:
            #   Checks if a file is valid with respect to the mapping file and the controlled vocabulary
            #   -----
            #   :param filename: File name of the file to be checked
            #   :param errors: Errors during the validation are returned in this output parameter
            #   :param warnings: Warnings during the validation are returned in this output parameter
            #   -----
            #   :raises:
            #     Exception: FileNotFound is thrown if the file could not be opened

