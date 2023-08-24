from Types cimport *
from MSExperiment cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *
from SpectrumMetaDataLookup cimport *

cdef extern from "<OpenMS/FORMAT/PercolatorOutfile.h>" namespace "OpenMS":
    
    cdef cppclass PercolatorOutfile "OpenMS::PercolatorOutfile":
        # wrap-doc:
            #  Class for reading Percolator tab-delimited output files
            #  
            #  For PSM-level output, the file extension should be ".psms"

        PercolatorOutfile() except + nogil 
        PercolatorOutfile(PercolatorOutfile &) except + nogil  # compiler

        # libcpp_string score_type_names()
        PercolatorOutfile_ScoreType getScoreType(String score_type_name) except + nogil  # wrap-doc:Returns a score type given its name

        void load(const String & filename, ProteinIdentification & proteins,
                  libcpp_vector[ PeptideIdentification ] & peptides,
                  SpectrumMetaDataLookup & lookup, 
                  PercolatorOutfile_ScoreType output_score) except + nogil  # wrap-doc:Loads a Percolator output file

cdef extern from "<OpenMS/FORMAT/PercolatorOutfile.h>" namespace "OpenMS::PercolatorOutfile":
    cdef enum PercolatorOutfile_ScoreType "OpenMS::PercolatorOutfile::ScoreType":
        #wrap-attach:
        #   PercolatorOutfile
        QVALUE
        POSTERRPROB
        SCORE
        SIZE_OF_SCORETYPE

