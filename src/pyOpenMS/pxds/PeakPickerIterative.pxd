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
        PeakPickerIterative(PeakPickerIterative &) nogil except + # compiler

        void pick(MSSpectrum & input,
                  MSSpectrum & output
                 ) nogil except +
            # wrap-doc:
                #   This will pick one single spectrum. The PeakPickerHiRes is used to
                #   generate seeds, these seeds are then used to re-center the mass and
                #   compute peak width and integrated intensity of the peak
                #   -----
                #   Finally, other peaks that would fall within the primary peak are
                #   discarded
                #   -----
                #   The output are the remaining peaks

        void pickExperiment(MSExperiment & input,
                            MSExperiment & output
                           ) nogil except +



# cdef extern from "<OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerIterative.h>" namespace "OpenMS":
# 
#     cdef cppclass PeakCandidate:
#         PeakCandidate() nogil except +
#         PeakCandidate(PeakCandidate) nogil except + # wrap-ignore
#         int index
#         double peak_apex_intensity
# 
#         double integrated_intensity
#         double leftWidth
#         double rightWidth
#         float mz
# 
