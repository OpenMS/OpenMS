from Types cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterHit "OpenMS::SiriusMzTabWriter::SiriusAdapterHit":
        SiriusAdapterHit() nogil except +
        SiriusAdapterHit(SiriusAdapterHit) nogil except + #wrap-ignore

        String formula
        String adduct
        String precursor_formula
        int rank
        double rank_score
        double iso_score
        double tree_score
        double sirius_score
        int explainedpeaks
        double explainedintensity

