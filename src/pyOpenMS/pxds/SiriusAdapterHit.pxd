from Types cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterHit "OpenMS::SiriusMzTabWriter::SiriusAdapterHit":
        SiriusAdapterHit() nogil except +
        SiriusAdapterHit(SiriusAdapterHit &) nogil except + # compiler

        String formula
        String adduct
        String precursor_formula
        int rank
        double iso_score
        double tree_score
        double sirius_score
        int explainedpeaks
        double explainedintensity
        double median_mass_error_fragment_peaks_ppm
        double median_absolute_mass_error_fragment_peaks_ppm
        double mass_error_precursor_ppm
