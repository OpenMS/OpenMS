from Types cimport *

from PeakFileOptions cimport *
from MzMLFile cimport *
from OnDiscMSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/IndexedMzMLFileLoader.h>" namespace "OpenMS":

    cdef cppclass IndexedMzMLFileLoader:

        IndexedMzMLFileLoader() except + nogil  # wrap-doc:A class to load an indexedmzML file
 
        bool load(String, OnDiscMSExperiment &) except + nogil 
            # wrap-doc:
                #  Load a file\n
                #  
                #  Tries to parse the file, success needs to be checked with the return value
                #
                #  
                #  :param filename: Filename determines where the file is located
                #  :param exp: Object which will contain the data after the call
                #  :return: Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed)

        void store(String, OnDiscMSExperiment &) except + nogil 
            # wrap-doc:
                #  Store a file from an on-disc data-structure
                #  
                #  
                #  :param filename: Filename determines where the file will be stored 
                #  :param exp: MS data to be stored

        void store(String, MSExperiment &) except + nogil 
            # wrap-doc:
                #  Store a file from an in-memory data-structure
                #  
                #  
                #  :param filename: Filename determines where the file will be stored 
                #  :param exp: MS data to be stored

        PeakFileOptions getOptions() except + nogil  # wrap-doc:Returns the options for loading/storing
        void setOptions(PeakFileOptions) except + nogil  # wrap-doc:Returns the options for loading/storing

