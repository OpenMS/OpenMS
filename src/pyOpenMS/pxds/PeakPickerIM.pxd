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

cdef extern from "<OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>" namespace "OpenMS::PeakPickerIM":
    
    # static members
    void pickIMTraces(MSSpectrum& s) except + nogil   # wrap-attach:PeakPickerIM wrap-doc:use trace detection for IM peak picking.
