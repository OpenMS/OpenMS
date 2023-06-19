from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from TargetedExperiment cimport *
from FASTAFile cimport *
from IntList cimport *
from PeptideIdentification cimport *
from FeatureMap cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/InclusionExclusionList.h>" namespace "OpenMS":

    cdef cppclass InclusionExclusionList(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        InclusionExclusionList() nogil except + # wrap-doc:Provides functionality for writing inclusion or exclusion lists
        InclusionExclusionList(InclusionExclusionList &) nogil except +
        void writeTargets(libcpp_vector[ FASTAEntry ] & fasta_entries, const String & out_path, IntList & charges, const String rt_model_path) nogil except + # wrap-doc:Writes inclusion or exclusion list of tryptic peptides of the given proteins (tab-delimited)
        void writeTargets(FeatureMap & map_, const String & out_path) nogil except + # wrap-doc:Writes inclusion or exclusion list of given feature map
        void writeTargets(libcpp_vector[ PeptideIdentification ] & pep_ids, const String & out_path, IntList & charges) nogil except + # wrap-doc:Writes inclusion or exclusion list of given peptide ids (tab-delimited)
