from Types cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>" namespace "OpenMS::SiriusMzTabWriter":
    
    cdef cppclass SiriusAdapterHit "OpenMS::SiriusMzTabWriter::SiriusAdapterHit":
        SiriusAdapterHit() nogil except +
        SiriusAdapterHit(SiriusAdapterHit) nogil except + #wrap-ignore

        String formula
        String adduct
        int rank
        double score
        double treescore
        double isoscore
        int explainedpeaks
        double explainedintensity
