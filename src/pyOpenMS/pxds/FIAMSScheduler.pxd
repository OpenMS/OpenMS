from Types cimport *
from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/FIAMSScheduler.h>" namespace "OpenMS":

    cdef cppclass FIAMSScheduler "OpenMS::FIAMSScheduler":
        #
        # wrap-doc:
        #     ADD PYTHON DOCUMENTATION HERE
        #

        FIAMSScheduler() nogil except + # wrap-doc:Scheduler for FIA-MS data batches. Works with FIAMSDataProcessor
        FIAMSScheduler(FIAMSScheduler &) nogil except +

        FIAMSScheduler(String filename, String base_dir, bool load_cached_) nogil except +
        void run() nogil except + # wrap-doc:Run the FIA-MS data analysis for the batch defined in the @filename_
        # libcpp_vector[ libcpp_map[ String, String ] ] getSamples() nogil except +
        String getBaseDir() nogil except + # wrap-doc:Returns the base directory for the relevant paths from the csv file
