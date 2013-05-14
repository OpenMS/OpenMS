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
        ChromatogramExtractor(ChromatogramExtractor)   nogil except + #wrap-ignore

        void extractChromatograms(MSExperiment[Peak1D, ChromatogramPeak] & input,
                                  MSExperiment[Peak1D, ChromatogramPeak] & output, 
                                  TargetedExperiment & transition_exp, 
                                  double extract_window,
                                  bool ppm,
                                  TransformationDescription trafo,
                                  double rt_extraction_window,
                                  String filter)


