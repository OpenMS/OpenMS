from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from CVTermList cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>" namespace "OpenMS":
    
    cdef cppclass IncludeExcludeTarget :
        IncludeExcludeTarget() nogil except +
        IncludeExcludeTarget(IncludeExcludeTarget) nogil except +
        void setName(String & name) nogil except +
        String  getName() nogil except +
        void setPeptideRef(String & peptide_ref) nogil except +
        String  getPeptideRef() nogil except +
        void setCompoundRef(String & compound_ref) nogil except +
        String  getCompoundRef() nogil except +
        void setPrecursorMZ(DoubleReal mz) nogil except +
        DoubleReal getPrecursorMZ() nogil except +
        void setPrecursorCVTermList(CVTermList & list_) nogil except +
        void addPrecursorCVTerm(CVTerm & cv_term) nogil except +
        CVTermList  getPrecursorCVTermList() nogil except +
        void setProductMZ(DoubleReal mz) nogil except +
        DoubleReal getProductMZ() nogil except +
        void setProductCVTermList(CVTermList & list_) nogil except +
        void addProductCVTerm(CVTerm & cv_term) nogil except +
        CVTermList  getProductCVTermList() nogil except +
        void setInterpretations(libcpp_vector[ CVTermList ] & interpretations) nogil except +
        libcpp_vector[ CVTermList ]  getInterpretations() nogil except +
        void addInterpretation(CVTermList & interpretation) nogil except +
        void setConfigurations(libcpp_vector[ Configuration ] & configuration) nogil except +
        libcpp_vector[ Configuration ]  getConfigurations() nogil except +
        void addConfiguration(Configuration & configuration) nogil except +
        void setPrediction(CVTermList & prediction) nogil except +
        void addPredictionTerm(CVTerm & prediction) nogil except +
        CVTermList  getPrediction() nogil except +
        void setRetentionTime(RetentionTime rt) nogil except +
        RetentionTime  getRetentionTime() nogil except +
        bool operator==(IncludeExcludeTarget & rhs) nogil except +
        bool operator!=(IncludeExcludeTarget & rhs) nogil except +

        # from CVTermsList
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except +
        void replaceCVTerm(CVTerm & term)               nogil except +

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                            String accession
                           ) nogil except +

        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map
                           ) nogil except +

        Map[String, libcpp_vector[CVTerm] ] getCVTerms()
        void addCVTerm(CVTerm & term)                   nogil except +

        bool hasCVTerm(String accession)  nogil except +
        bool empty()                      nogil except +
