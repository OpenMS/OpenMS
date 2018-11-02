from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from FileHandler cimport *
from DataValue cimport *
from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>" namespace "OpenMS":
    
    cdef cppclass PercolatorFeatureSetHelper "OpenMS::PercolatorFeatureSetHelper":

        PercolatorFeatureSetHelper(PercolatorFeatureSetHelper) nogil except + #wrap-ignore

        void concatMULTISEPeptideIds(libcpp_vector[ PeptideIdentification ] & all_peptide_ids, libcpp_vector[ PeptideIdentification ] & new_peptide_ids, String search_engine) nogil except +
        void mergeMULTISEPeptideIds(libcpp_vector[ PeptideIdentification ] & all_peptide_ids, libcpp_vector[ PeptideIdentification ] & new_peptide_ids, String search_engine) nogil except +
        void mergeMULTISEProteinIds(libcpp_vector[ ProteinIdentification ] & all_protein_ids, libcpp_vector[ ProteinIdentification ] & new_protein_ids) nogil except +
        void addMSGFFeatures(libcpp_vector[ PeptideIdentification ] & peptide_ids, StringList & feature_set) nogil except +
        void addXTANDEMFeatures(libcpp_vector[ PeptideIdentification ] & peptide_ids, StringList & feature_set) nogil except +
        void addCOMETFeatures(libcpp_vector[ PeptideIdentification ] & peptide_ids, StringList & feature_set) nogil except +
        void addMASCOTFeatures(libcpp_vector[ PeptideIdentification ] & peptide_ids, StringList & feature_set) nogil except +
        void addMULTISEFeatures(libcpp_vector[ PeptideIdentification ] & peptide_ids, StringList & search_engines_used, StringList & feature_set, bool complete_only, bool limits_imputation) nogil except +
        void addCONCATSEFeatures(libcpp_vector[ PeptideIdentification ] & peptide_id_list, StringList & search_engines_used, StringList & feature_set) nogil except +
        void checkExtraFeatures(libcpp_vector[ PeptideHit ] & psms, StringList & extra_features) nogil except +

