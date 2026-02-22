from libcpp.vector cimport vector

from DefaultParamHandler cimport DefaultParamHandler
from PeptideIdentificationList cimport PeptideIdentificationList
from MSExperiment cimport MSExperiment
from FASTAFile cimport FASTAFile
from ProteinIdentification cimport ProteinIdentification
from Types cimport Size, UInt
from String cimport String


cdef extern from "<OpenMS/QC/DBSuitability.h>" namespace "OpenMS":

    cdef cppclass DBSuitability(DefaultParamHandler):

        cdef cppclass SuitabilityData:
            Size num_top_novo
            Size num_top_db
            Size num_interest
            Size num_re_ranked
            double cut_off
            double suitability
            double suitability_no_rerank
            double suitability_corr_no_rerank

        DBSuitability() except +

        void compute(PeptideIdentificationList&& pep_ids,
                     const MSExperiment& exp,
                     const vector[FASTAFile.FASTAEntry]& original_fasta,
                     const vector[FASTAFile.FASTAEntry]& novo_fasta,
                     const ProteinIdentification.SearchParameters& search_params) except +

        const vector[SuitabilityData]& getResults() const except +