from Types cimport *
from MSExperiment cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from SpectrumMetaDataLookup cimport *

cdef extern from "<OpenMS/FORMAT/PercolatorOutfile.h>" namespace "OpenMS":
    
    cdef cppclass PercolatorOutfile "OpenMS::PercolatorOutfile":
        PercolatorOutfile() nogil except +
        PercolatorOutfile(PercolatorOutfile) nogil except + #wrap-ignore

        # libcpp_string score_type_names()
        PercolatorOutfile_ScoreType getScoreType(String score_type_name) nogil except +

        void load(const String & filename, ProteinIdentification & proteins,
                  libcpp_vector[ PeptideIdentification ] & peptides,
                  SpectrumMetaDataLookup & lookup, 
                  PercolatorOutfile_ScoreType output_score) nogil except +

cdef extern from "<OpenMS/FORMAT/PercolatorOutfile.h>" namespace "OpenMS::PercolatorOutfile":
    cdef enum PercolatorOutfile_ScoreType "OpenMS::PercolatorOutfile::ScoreType":
        #wrap-attach:
        #    PercolatorOutfile
        QVALUE
        POSTERRPROB
        SCORE
        SIZE_OF_SCORETYPE

