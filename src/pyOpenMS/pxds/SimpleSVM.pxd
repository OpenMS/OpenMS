from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from String cimport *

ctypedef libcpp_map[String, libcpp_vector[double] ] PredictorMap

cdef extern from "<OpenMS/ANALYSIS/SVM/SimpleSVM.h>" namespace "OpenMS":
    
    cdef cppclass SimpleSVM(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #   Simple interface to support vector machines for classification (via LIBSVM)
        
        SimpleSVM() nogil except + 
        SimpleSVM(SimpleSVM &) nogil except + # compiler

        # Currently nested things inside std::map doesn't work
        # void setup(libcpp_map[String, libcpp_vector[double] ]& predictors, libcpp_map[ Size, Int ] & labels) nogil except +
        void predict(libcpp_vector[ SVMPrediction ] & predictions, libcpp_vector[ size_t ] indexes) nogil except +
            # wrap-doc:
                #   Predict class labels (and probabilities)
                #   -----
                #   :param predictions: Output vector of prediction results (same order as ``)
                #   :param indexes: Vector of observation indexes for which predictions are desired. If empty (default), predictions are made for all observations

        void getFeatureWeights(libcpp_map[ String, double ] & feature_weights) nogil except +
            # wrap-doc:
                #   Returns the weights used for features (predictors) in the SVM model
                #   -----
                #   Currently only supported for two-class classification
                #   If a linear kernel is used, the weights are informative for ranking features

        void writeXvalResults(const String & path) nogil except +


cdef extern from "<OpenMS/ANALYSIS/SVM/SimpleSVM.h>" namespace "OpenMS::SimpleSVM":
    
    cdef cppclass SVMPrediction "OpenMS::SimpleSVM::Prediction":

        SVMPrediction() nogil except + # compiler
        SVMPrediction(SVMPrediction &) nogil except + # compiler

        Int label
        libcpp_map[int, double ] probabilities
