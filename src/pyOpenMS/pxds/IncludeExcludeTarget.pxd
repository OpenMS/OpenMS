from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from CVTermList cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>" namespace "OpenMS":
    
    cdef cppclass IncludeExcludeTarget :
        IncludeExcludeTarget() nogil except + # wrap-doc:This class stores a SRM/MRM transition
        IncludeExcludeTarget(IncludeExcludeTarget &) nogil except + # TODO
        void setName(const String & name) nogil except + # TODO
        String  getName() nogil except + # TODO
        void setPeptideRef(const String & peptide_ref) nogil except + # TODO
        String  getPeptideRef() nogil except + # TODO
        void setCompoundRef(const String & compound_ref) nogil except + # TODO
        String  getCompoundRef() nogil except + # TODO
        void setPrecursorMZ(double mz) nogil except + # TODO
        double getPrecursorMZ() nogil except + # TODO
        void setPrecursorCVTermList(CVTermList & list_) nogil except + # TODO
        void addPrecursorCVTerm(CVTerm & cv_term) nogil except + # TODO
        CVTermList  getPrecursorCVTermList() nogil except + # TODO
        void setProductMZ(double mz) nogil except + # TODO
        double getProductMZ() nogil except + # TODO
        void setProductCVTermList(CVTermList & list_) nogil except + # TODO
        void addProductCVTerm(CVTerm & cv_term) nogil except + # TODO
        CVTermList  getProductCVTermList() nogil except + # TODO
        void setInterpretations(libcpp_vector[ CVTermList ] & interpretations) nogil except + # TODO
        libcpp_vector[ CVTermList ]  getInterpretations() nogil except + # TODO
        void addInterpretation(CVTermList & interpretation) nogil except + # TODO
        void setConfigurations(libcpp_vector[ Configuration ] & configuration) nogil except + # TODO
        libcpp_vector[ Configuration ]  getConfigurations() nogil except + # TODO
        void addConfiguration(Configuration & configuration) nogil except + # TODO
        void setPrediction(CVTermList & prediction) nogil except + # TODO
        void addPredictionTerm(CVTerm & prediction) nogil except + # TODO
        CVTermList  getPrediction() nogil except + # TODO
        void setRetentionTime(RetentionTime rt) nogil except + # TODO
        RetentionTime  getRetentionTime() nogil except + # TODO
        bool operator==(IncludeExcludeTarget & rhs) nogil except + # TODO
        bool operator!=(IncludeExcludeTarget & rhs) nogil except + # TODO

        # from CVTermsList
        void setCVTerms(libcpp_vector[CVTerm] & terms)  nogil except + # TODO
        void replaceCVTerm(CVTerm & term)               nogil except + # TODO

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                            String accession
                           ) nogil except + # TODO

        void replaceCVTerms(Map[String, libcpp_vector[CVTerm] ] cv_term_map
                           ) nogil except + # TODO

        Map[String, libcpp_vector[CVTerm] ] getCVTerms() nogil except + # TODO
        void addCVTerm(CVTerm & term)                   nogil except + # TODO

        bool hasCVTerm(String accession)  nogil except + # TODO
        bool empty()                      nogil except + # TODO
