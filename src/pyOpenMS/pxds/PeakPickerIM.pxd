from libcpp cimport bool
from Types cimport *
from MSSpectrum cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>" namespace "OpenMS":

    cdef cppclass PeakPickerIM(DefaultParamHandler):
        # wrap-doc:
        #  Peak picking algorithm for ion mobility data

        PeakPickerIM() except + nogil
        PeakPickerIM(PeakPickerIM &) except + nogil
        void pickIMTraces(MSSpectrum& s) except + nogil  # wrap-doc:Use trace detection for IM peak picking.
        void pickIMCluster(MSSpectrum& s) except + nogil  # wrap-doc:Use clustering for IM peak picking.
        void pickIMElutionProfiles(MSSpectrum& s) except + nogil  # wrap-doc:Use elution profile detection for IM peak picking.
