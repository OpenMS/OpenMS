from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from ProgressLogger cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>" namespace "OpenMS":

    cdef cppclass ChromatogramExtractor(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        ChromatogramExtractor()                  nogil except +
        ChromatogramExtractor(ChromatogramExtractor)   nogil except + 

        void extractChromatograms(MSExperiment & input,
                                  MSExperiment & output, 
                                  TargetedExperiment & transition_exp,
                                  double extract_window,
                                  bool ppm,
                                  TransformationDescription trafo,
                                  double rt_extraction_window,
                                  String filter) nogil except +

        # TODO immutable types by reference
        # void extract_value_tophat(MSSpectrum input, double mz,
        #  Size peak_idx, double integrated_intensity, double extract_window, bool ppm)
        # void extract_value_bartlett(MSSpectrum input, double mz,
        #  Size peak_idx, double integrated_intensity, double extract_window, bool ppm)
    
