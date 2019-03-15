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
        InclusionExclusionList() nogil except +
        InclusionExclusionList(InclusionExclusionList) nogil except + #wrap-ignore
        void writeTargets(libcpp_vector[ FASTAEntry ] & fasta_entries, const String & out_path, IntList & charges, const String rt_model_path) nogil except +
        void writeTargets(FeatureMap & map_, const String & out_path) nogil except +
        void writeTargets(libcpp_vector[ PeptideIdentification ] & pep_ids, const String & out_path, IntList & charges) nogil except +

