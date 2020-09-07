from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS":

    cdef cppclass ReactionMonitoringTransition(CVTermList):
        # wrap-inherits:
        #    CVTermList

        ReactionMonitoringTransition()                  nogil except +
        ReactionMonitoringTransition(ReactionMonitoringTransition)   nogil except + #wrap-ignore
        String getName()                           nogil except +
        String getNativeID()                       nogil except +
        String getPeptideRef()                     nogil except +
        void setName(String name)                  nogil except +
        void setNativeID(String name)              nogil except +
        void setPeptideRef(String peptide_ref)     nogil except +

        double getProductMZ()         nogil except +
        void setProductMZ(double)     nogil except +

        double getPrecursorMZ()       nogil except +
        void setPrecursorMZ(double)   nogil except +
      
        DecoyTransitionType getDecoyTransitionType() nogil except +

        void setCompoundRef(const String & compound_ref)nogil except +
        String  getCompoundRef()nogil except +

        bool hasPrecursorCVTerms() nogil except +
        void setPrecursorCVTermList(CVTermList & list_)nogil except +
        void addPrecursorCVTerm(CVTerm & cv_term)nogil except +
        CVTermList getPrecursorCVTermList()nogil except +

        void addProductCVTerm(CVTerm & cv_term)nogil except +

        libcpp_vector[ TraMLProduct ]  getIntermediateProducts()nogil except +
        void addIntermediateProduct(TraMLProduct product)nogil except +
        void setIntermediateProducts(libcpp_vector[ TraMLProduct ] & products)nogil except +

        void setProduct(TraMLProduct product)nogil except +
        TraMLProduct getProduct()nogil except +
        void setRetentionTime(RetentionTime rt)nogil except +
        RetentionTime getRetentionTime()nogil except +

        void setPrediction(Prediction & prediction)nogil except +
        void addPredictionTerm(CVTerm & prediction)nogil except +
        bool hasPrediction() nogil except +
        Prediction getPrediction()nogil except +

        void setDecoyTransitionType(DecoyTransitionType & d)nogil except +

        double getLibraryIntensity()nogil except +
        void setLibraryIntensity(double intensity)nogil except +

        int getProductChargeState() nogil except +
        bool isProductChargeStateSet() nogil except +

        bool isDetectingTransition() nogil except +
        void setDetectingTransition(bool val) nogil except +

        bool isIdentifyingTransition() nogil except +
        void setIdentifyingTransition(bool val) nogil except +

        bool isQuantifyingTransition() nogil except +
        void setQuantifyingTransition(bool val) nogil except +

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS::ReactionMonitoringTransition":

    cdef enum DecoyTransitionType:

        UNKNOWN, TARGET, DECOY

