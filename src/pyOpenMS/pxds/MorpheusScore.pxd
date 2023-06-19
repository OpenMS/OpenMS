from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/MorpheusScore.h>" namespace "OpenMS":
    
    cdef cppclass MorpheusScore "OpenMS::MorpheusScore":
        MorpheusScore() nogil except + # compiler
        MorpheusScore(MorpheusScore &) nogil except + # compiler
        MorpheusScore_Result compute(double fragment_mass_tolerance,
                                     bool fragment_mass_tolerance_unit_ppm,
                                     const MSSpectrum & exp_spectrum,
                                     const MSSpectrum & theo_spectrum) nogil except + # wrap-doc:Returns Morpheus Score

cdef extern from "<OpenMS/ANALYSIS/RNPXL/MorpheusScore.h>" namespace "OpenMS::MorpheusScore":
    
    cdef cppclass MorpheusScore_Result "OpenMS::MorpheusScore::Result":
        MorpheusScore_Result() nogil except + # compiler
        MorpheusScore_Result(MorpheusScore_Result &) nogil except + # compiler
        Size matches
        Size n_peaks
        float score
        float MIC
        float TIC
        float err

