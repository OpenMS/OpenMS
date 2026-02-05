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
        # wrap-doc:
        #  The ChromatogramExtractor extracts chromatograms (intensity vs
        #  retention time) from mass spectrometry data. * * This class provides
        #  functionality to extract chromatographic traces from mass spectrometry
        #  data * based on specified coordinates (m/z, retention time, and
        #  optionally ion mobility values). * * The extractor supports two main
        #  interfaces: * 1. Legacy interface: Takes a TargetedExperiment object
        #  containing transitions and extracts * chromatograms at the m/z values
        #  specified in those transitions. * 2. Modern interface: Takes a set of
        #  ExtractionCoordinates that specify the exact coordinates * for
        #  extraction. This provides more flexibility and control over the
        #  extraction process. * The prepare_coordinates() helper function can
        #  generate these coordinates for common * MS1 and MS2 extraction
        #  scenarios. * * Key features: * - Supports both MS1 and MS2 level
        #  extractions * - Configurable extraction window sizes in m/z dimension
        #  (absolute or ppm) * - Multiple filter types available (Bartlett,
        #  tophat) for signal processing * - Handles ion mobility data when
        #  available * - Optimized for SWATH/DIA (Data Independent Acquisition)
        #  experiments * - Progress logging capabilities through ProgressLogger
        #  base class * * For MS2 extractions, the input data is expected to come
        #  from a SWATH/DIA experiment * where precursor ions are fragmented in
        #  wide isolation windows, allowing for extraction * of fragment ion
        #  chromatograms. * * @see ChromatogramExtractorAlgorithm For the
        #  underlying extraction algorithm implementation * @see
        #  ExtractionCoordinates For the coordinate specification format
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

        # static members
        @staticmethod
        void prepare_coordinates(
            libcpp_vector[ shared_ptr[OSChromatogram] ] & output_chromatograms,
            libcpp_vector[ ExtractionCoordinates ] & extraction_coordinates,
            TargetedExperiment & targeted,
            double rt_extraction_window,
            bool ms1,
            int ms1_isotopes) except + nogil

