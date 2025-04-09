from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h>" namespace "OpenMS":

    cdef cppclass PeakPickerIterative(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger

        PeakPickerIterative() except + nogil 
        PeakPickerIterative(PeakPickerIterative &) except + nogil  # compiler

        void pick(MSSpectrum & input,
                  MSSpectrum & output
                 ) except + nogil 
            # wrap-doc:
                #  This will pick one single spectrum. The PeakPickerHiRes is used to
                #  generate seeds, these seeds are then used to re-center the mass and
                #  compute peak width and integrated intensity of the peak
                #  
                #  Finally, other peaks that would fall within the primary peak are
                #  discarded
                #  
                #  The output are the remaining peaks

        void pickExperiment(MSExperiment & input,
                            MSExperiment & output
                           ) except + nogil 



# cdef extern from "<OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h>" namespace "OpenMS":
# 
#    cdef cppclass PeakCandidate:
#        PeakCandidate() except + nogil 
#        PeakCandidate(PeakCandidate) except + nogil  # wrap-ignore
#        int index
#        double peak_apex_intensity
# 
#        double integrated_intensity
#        double leftWidth
#        double rightWidth
#        float mz
# 
