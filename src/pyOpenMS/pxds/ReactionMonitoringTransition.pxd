from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS":

    cdef cppclass ReactionMonitoringTransition(CVTermList):
        # wrap-inherits:
        #   CVTermList
        # wrap-doc:
        #  This class stores a SRM/MRM transition
        #  
        #  This class is capable of representing a <Transition> tag in a TraML
        #  document completely and contains all associated information
        #  
        #  The default values for precursor m/z is 0.0 which indicates that it is
        #  uninitialized

        ReactionMonitoringTransition() except + nogil 
        ReactionMonitoringTransition(ReactionMonitoringTransition &)   except + nogil 
        String getName() except + nogil 
        String getNativeID() except + nogil 
        String getPeptideRef() except + nogil 
        void setName(String name) except + nogil 
        void setNativeID(String name) except + nogil 
        void setPeptideRef(String peptide_ref) except + nogil 

        double getProductMZ() except + nogil 
        void setProductMZ(double) except + nogil 

        double getPrecursorMZ() except + nogil  # wrap-doc:Returns the precursor mz (Q1 value)
        void setPrecursorMZ(double) except + nogil  # wrap-doc:Sets the precursor mz (Q1 value)
      
        DecoyTransitionType getDecoyTransitionType() except + nogil  # wrap-doc:Returns the type of transition (target or decoy)

        void setCompoundRef(const String & compound_ref) except + nogil 
        String  getCompoundRef() except + nogil 

        bool hasPrecursorCVTerms() except + nogil  # wrap-doc:Returns true if precursor CV Terms exist (means it is safe to call getPrecursorCVTermList)
        void setPrecursorCVTermList(CVTermList & list_) except + nogil  # wrap-doc:Sets a list of precursor CV Terms
        void addPrecursorCVTerm(CVTerm & cv_term) except + nogil  # wrap-doc:Adds precursor CV Term
        CVTermList getPrecursorCVTermList() except + nogil  # wrap-doc:Obtains the list of CV Terms for the precursor

        void addProductCVTerm(CVTerm & cv_term) except + nogil 

        libcpp_vector[ TraMLProduct ]  getIntermediateProducts() except + nogil 
        void addIntermediateProduct(TraMLProduct product) except + nogil 
        void setIntermediateProducts(libcpp_vector[ TraMLProduct ] & products) except + nogil 

        void setProduct(TraMLProduct product) except + nogil 
        TraMLProduct getProduct() except + nogil 
        void setRetentionTime(RetentionTime rt) except + nogil 
        RetentionTime getRetentionTime() except + nogil 

        void setPrediction(Prediction & prediction) except + nogil  # wrap-doc:Sets prediction
        void addPredictionTerm(CVTerm & prediction) except + nogil  # wrap-doc:Adds prediction term
        bool hasPrediction() except + nogil  # wrap-doc:Returns true if a Prediction object exists (means it is safe to call getPrediction)
        Prediction getPrediction() except + nogil  # wrap-doc:Obtains the Prediction object 

        void setDecoyTransitionType(DecoyTransitionType & d) except + nogil  # wrap-doc:Sets the type of transition (target or decoy)

        double getLibraryIntensity() except + nogil  # wrap-doc:Returns the library intensity (ion count or normalized ion count from a spectral library)
        void setLibraryIntensity(double intensity) except + nogil  # wrap-doc:Sets the library intensity (ion count or normalized ion count from a spectral library)

        int getProductChargeState() except + nogil  # wrap-doc:Returns the charge state of the product
        bool isProductChargeStateSet() except + nogil  # wrap-doc:Returns true if charge state of product is already set

        bool isDetectingTransition() except + nogil 
        void setDetectingTransition(bool val) except + nogil 

        bool isIdentifyingTransition() except + nogil 
        void setIdentifyingTransition(bool val) except + nogil 

        bool isQuantifyingTransition() except + nogil 
        void setQuantifyingTransition(bool val) except + nogil 

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS::ReactionMonitoringTransition":

    cdef enum DecoyTransitionType:

        UNKNOWN, TARGET, DECOY
