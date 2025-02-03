from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from ProgressLogger cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *
from libcpp cimport bool
from ChromatogramExtractorAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>" namespace "OpenMS":

    cdef cppclass ChromatogramExtractor(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger
    
        ChromatogramExtractor() except + nogil  # compiler
        ChromatogramExtractor(ChromatogramExtractor &) except + nogil  # compiler

        void extractChromatograms(
            shared_ptr[ SpectrumAccessOpenMS ] input,
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output,
            libcpp_vector[ ExtractionCoordinates ] extraction_coordinates,
            double mz_extraction_window,
            bool ppm,
            double im_extraction_window,
            String filter) except + nogil

# COMMENT: wrap static methods
cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>" namespace "OpenMS::ChromatogramExtractor":
        
        # static members
        void prepare_coordinates(
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output_chromatograms,
            libcpp_vector[ ExtractionCoordinates ] & extraction_coordinates,
            TargetedExperiment & targeted,
            double rt_extraction_window,
            bool ms1,
            int ms1_isotopes) except + nogil # wrap-attach:ChromatogramExtractor
        



