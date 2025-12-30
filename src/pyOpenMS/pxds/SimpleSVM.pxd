from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool
from String cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ML/SVM/SimpleSVM.h>" namespace "OpenMS":
    
    cdef cppclass SimpleSVM(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        SimpleSVM() except + nogil 
        SimpleSVM(SimpleSVM &) except + nogil  # wrap-ignore

        void setup(libcpp_map[String, libcpp_vector[double]] & predictors,
                   libcpp_map[size_t, double] & outcomes,
                   bool classification) except + nogil 
            # wrap-doc:
                #  Load data and train a model
                #  
                #  
                #  :param predictors: Mapping from predictor name to vector of predictor values (for different observations). All vectors should have the same length; values will be changed by scaling
                #  :param outcomes: Mapping from observation index to class label or regression value in the training set
                #  :param classification: True (default) if SVM classification should be used, SVR otherwise

        void predict(libcpp_vector[SimpleSVM_Prediction] & predictions,
                     libcpp_vector[size_t] indexes) except + nogil 
            # wrap-doc:
                #  Predict class labels or regression values (and probabilities)
                #  
                #  
                #  :param predictions: Output vector of prediction results (same order as indexes)
                #  :param indexes: Vector of observation indexes for which predictions are desired. If empty (default), predictions are made for all observations

        void predict(libcpp_map[String, libcpp_vector[double]] & predictors,
                     libcpp_vector[SimpleSVM_Prediction] & predictions) except + nogil 
            # wrap-doc:
                #  Predict class labels or regression values (and probabilities)
                #  
                #  
                #  :param predictors: Mapping from predictor name to vector of predictor values (for different observations). All vectors should have the same length; values will be changed by scaling applied to training data in setup
                #  :param predictions: Output vector of prediction results (same order as indexes)

        void getFeatureWeights(libcpp_map[String, double] & feature_weights) except + nogil 
            # wrap-doc:
                #  Get the weights used for features (predictors) in the SVM model
                #  
                #  Currently only supported for two-class classification
                #  If a linear kernel is used, the weights are informative for ranking features

        void writeXvalResults(String path) except + nogil 
            # wrap-doc:Write cross-validation (parameter optimization) results to a CSV file

        libcpp_map[String, libcpp_pair[double, double]] getScaling() except + nogil 
            # wrap-doc:Get data range of predictors before scaling to [0, 1]


cdef extern from "<OpenMS/ML/SVM/SimpleSVM.h>" namespace "OpenMS::SimpleSVM":

    cdef cppclass SimpleSVM_Prediction "OpenMS::SimpleSVM::Prediction":
        SimpleSVM_Prediction() except + nogil 
        SimpleSVM_Prediction(SimpleSVM_Prediction &) except + nogil  # compiler

        double outcome
        libcpp_map[double, double] probabilities