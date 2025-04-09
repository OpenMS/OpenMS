from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from CVTermList cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>" namespace "OpenMS":
    
    cdef cppclass IncludeExcludeTarget :
        IncludeExcludeTarget() except + nogil  # wrap-doc:This class stores a SRM/MRM transition
        IncludeExcludeTarget(IncludeExcludeTarget &) except + nogil  # TODO
        void setName(const String & name) except + nogil  # TODO
        String  getName() except + nogil  # TODO
        void setPeptideRef(const String & peptide_ref) except + nogil  # TODO
        String  getPeptideRef() except + nogil  # TODO
        void setCompoundRef(const String & compound_ref) except + nogil  # TODO
        String  getCompoundRef() except + nogil  # TODO
        void setPrecursorMZ(double mz) except + nogil  # TODO
        double getPrecursorMZ() except + nogil  # TODO
        void setPrecursorCVTermList(CVTermList & list_) except + nogil  # TODO
        void addPrecursorCVTerm(CVTerm & cv_term) except + nogil  # TODO
        CVTermList  getPrecursorCVTermList() except + nogil  # TODO
        void setProductMZ(double mz) except + nogil  # TODO
        double getProductMZ() except + nogil  # TODO
        void setProductCVTermList(CVTermList & list_) except + nogil  # TODO
        void addProductCVTerm(CVTerm & cv_term) except + nogil  # TODO
        CVTermList  getProductCVTermList() except + nogil  # TODO
        void setInterpretations(libcpp_vector[ CVTermList ] & interpretations) except + nogil  # TODO
        libcpp_vector[ CVTermList ]  getInterpretations() except + nogil  # TODO
        void addInterpretation(CVTermList & interpretation) except + nogil  # TODO
        void setConfigurations(libcpp_vector[ Configuration ] & configuration) except + nogil  # TODO
        libcpp_vector[ Configuration ]  getConfigurations() except + nogil  # TODO
        void addConfiguration(Configuration & configuration) except + nogil  # TODO
        void setPrediction(CVTermList & prediction) except + nogil  # TODO
        void addPredictionTerm(CVTerm & prediction) except + nogil  # TODO
        CVTermList  getPrediction() except + nogil  # TODO
        void setRetentionTime(RetentionTime rt) except + nogil  # TODO
        RetentionTime  getRetentionTime() except + nogil  # TODO
        bool operator==(IncludeExcludeTarget & rhs) except + nogil  # TODO
        bool operator!=(IncludeExcludeTarget & rhs) except + nogil  # TODO

        # from CVTermsList
        void setCVTerms(libcpp_vector[CVTerm] & terms)  except + nogil  # TODO
        void replaceCVTerm(CVTerm & term)               except + nogil  # TODO

        void replaceCVTerms(libcpp_vector[CVTerm] cv_terms,
                            String accession
                           ) except + nogil  # TODO

        void replaceCVTerms(libcpp_map[String, libcpp_vector[CVTerm] ] cv_term_map
                           ) except + nogil  # TODO

        libcpp_map[String, libcpp_vector[CVTerm] ] getCVTerms() except + nogil  # TODO
        void addCVTerm(CVTerm & term)                   except + nogil  # TODO

        bool hasCVTerm(String accession)  except + nogil  # TODO
        bool empty()                      except + nogil  # TODO
