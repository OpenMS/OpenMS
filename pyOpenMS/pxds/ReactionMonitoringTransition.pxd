from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from String cimport *
from TargetedExperimentHelper cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS":

    cdef cppclass ReactionMonitoringTransition:

        ReactionMonitoringTransition()                  nogil except +
        ReactionMonitoringTransition(ReactionMonitoringTransition)   nogil except + #wrap-ignore
        String getName()                           nogil except +
        String getNativeID()                       nogil except +
        String getPeptideRef()                     nogil except +
        void setName(String name)                  nogil except +
        void setNativeID(String name)              nogil except +
        void setPeptideRef(String peptide_ref)     nogil except +

        DoubleReal getProductMZ()         nogil except +
        void setProductMZ(DoubleReal)     nogil except +

        DoubleReal getPrecursorMZ()       nogil except +
        void setPrecursorMZ(DoubleReal)   nogil except +
      
        DecoyTransitionType getDecoyTransitionType() nogil except +

        void setCompoundRef(String & compound_ref)nogil except +
        String  getCompoundRef()nogil except +
        void setPrecursorCVTermList(CVTermList & list_)nogil except +
        void addPrecursorCVTerm(CVTerm & cv_term)nogil except +
        CVTermList  getPrecursorCVTermList()nogil except +
        void addProductCVTerm(CVTerm & cv_term)nogil except +
        libcpp_vector[ TraMLProduct ]  getIntermediateProducts()nogil except +
        void addIntermediateProduct(TraMLProduct product)nogil except +
        void setIntermediateProducts(libcpp_vector[ TraMLProduct ] & products)nogil except +
        void setProduct(TraMLProduct product)nogil except +
        TraMLProduct  getProduct()nogil except +
        void setRetentionTime(RetentionTime rt)nogil except +
        RetentionTime  getRetentionTime()nogil except +
        void setPrediction(Prediction & prediction)nogil except +
        void addPredictionTerm(CVTerm & prediction)nogil except +
        Prediction  getPrediction()nogil except +
        void setDecoyTransitionType(DecoyTransitionType & d)nogil except +
        DoubleReal getLibraryIntensity()nogil except +
        void setLibraryIntensity(DoubleReal intensity)nogil except +

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS::ReactionMonitoringTransition":

    cdef enum DecoyTransitionType:

        UNKNOWN, TARGET, DECOY

