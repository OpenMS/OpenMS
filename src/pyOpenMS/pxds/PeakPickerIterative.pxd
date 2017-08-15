from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerIterative.h>" namespace "OpenMS":

    cdef cppclass PeakPickerIterative(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        PeakPickerIterative() nogil except +
        PeakPickerIterative(PeakPickerIterative) nogil except + #wrap-ignore

        void pick(MSSpectrum & input,
                  MSSpectrum & output
                 ) nogil except +

        void pickExperiment(MSExperiment & input,
                            MSExperiment & output
                           ) nogil except +



# cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerIterative.h>" namespace "OpenMS":
# 
#     cdef cppclass PeakCandidate:
#         PeakCandidate() nogil except +
#         PeakCandidate(PeakCandidate) nogil except +
#         int index
#         double peak_apex_intensity
# 
#         double integrated_intensity
#         double leftWidth
#         double rightWidth
#         float mz
# 
