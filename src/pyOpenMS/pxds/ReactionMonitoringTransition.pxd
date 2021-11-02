from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS":

    cdef cppclass ReactionMonitoringTransition(CVTermList):
        # wrap-inherits:
        #    CVTermList
        # wrap-doc:
        #   This class stores a SRM/MRM transition
        #   -----
        #   This class is capable of representing a <Transition> tag in a TraML
        #   document completely and contains all associated information
        #   -----
        #   The default values for precursor m/z is 0.0 which indicates that it is
        #   uninitialized

        ReactionMonitoringTransition() nogil except +
        ReactionMonitoringTransition(ReactionMonitoringTransition &)   nogil except +
        String getName() nogil except +
        String getNativeID() nogil except +
        String getPeptideRef() nogil except +
        void setName(String name) nogil except +
        void setNativeID(String name) nogil except +
        void setPeptideRef(String peptide_ref) nogil except +

        double getProductMZ() nogil except +
        void setProductMZ(double) nogil except +

        double getPrecursorMZ() nogil except + # wrap-doc:Returns the precursor mz (Q1 value)
        void setPrecursorMZ(double) nogil except + # wrap-doc:Sets the precursor mz (Q1 value)
      
        DecoyTransitionType getDecoyTransitionType() nogil except + # wrap-doc:Returns the type of transition (target or decoy)

        void setCompoundRef(const String & compound_ref)nogil except +
        String  getCompoundRef()nogil except +

        bool hasPrecursorCVTerms() nogil except + # wrap-doc:Returns true if precursor CV Terms exist (means it is safe to call getPrecursorCVTermList)
        void setPrecursorCVTermList(CVTermList & list_)nogil except + # wrap-doc:Sets a list of precursor CV Terms
        void addPrecursorCVTerm(CVTerm & cv_term)nogil except + # wrap-doc:Adds precursor CV Term
        CVTermList getPrecursorCVTermList()nogil except + # wrap-doc:Obtains the list of CV Terms for the precursor

        void addProductCVTerm(CVTerm & cv_term)nogil except +

        libcpp_vector[ TraMLProduct ]  getIntermediateProducts()nogil except +
        void addIntermediateProduct(TraMLProduct product)nogil except +
        void setIntermediateProducts(libcpp_vector[ TraMLProduct ] & products)nogil except +

        void setProduct(TraMLProduct product)nogil except +
        TraMLProduct getProduct()nogil except +
        void setRetentionTime(RetentionTime rt)nogil except +
        RetentionTime getRetentionTime()nogil except +

        void setPrediction(Prediction & prediction)nogil except + # wrap-doc:Sets prediction
        void addPredictionTerm(CVTerm & prediction)nogil except + # wrap-doc:Adds prediction term
        bool hasPrediction() nogil except + # wrap-doc:Returns true if a Prediction object exists (means it is safe to call getPrediction)
        Prediction getPrediction()nogil except + # wrap-doc:Obtains the Prediction object 

        void setDecoyTransitionType(DecoyTransitionType & d)nogil except + # wrap-doc:Sets the type of transition (target or decoy)

        double getLibraryIntensity()nogil except + # wrap-doc:Returns the library intensity (ion count or normalized ion count from a spectral library)
        void setLibraryIntensity(double intensity)nogil except + # wrap-doc:Sets the library intensity (ion count or normalized ion count from a spectral library)

        int getProductChargeState() nogil except + # wrap-doc:Returns the charge state of the product
        bool isProductChargeStateSet() nogil except + # wrap-doc:Returns true if charge state of product is already set

        bool isDetectingTransition() nogil except +
        void setDetectingTransition(bool val) nogil except +

        bool isIdentifyingTransition() nogil except +
        void setIdentifyingTransition(bool val) nogil except +

        bool isQuantifyingTransition() nogil except +
        void setQuantifyingTransition(bool val) nogil except +

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS::ReactionMonitoringTransition":

    cdef enum DecoyTransitionType:

        UNKNOWN, TARGET, DECOY
