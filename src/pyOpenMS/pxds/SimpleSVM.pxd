from Types cimport *
from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *

# typedef std::map<String, std::vector<double> > PredictorMap;
cdef extern from "<OpenMS/ANALYSIS/SVM/SimpleSVM.h>" namespace "OpenMS":
    
    cdef cppclass SimpleSVM(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        SimpleSVM() nogil except +
        SimpleSVM(SimpleSVM) nogil except + #wrap-ignore
        # void setup(PredictorMap & predictors, libcpp_map[ Size, Int ] & labels) nogil except +
        void predict(libcpp_vector[ SVMPrediction ] & predictions, libcpp_vector[ size_t ] indexes) nogil except +
        void getFeatureWeights(libcpp_map[ String, double ] & feature_weights) nogil except +
        void writeXvalResults(String & path) nogil except +


cdef extern from "<OpenMS/ANALYSIS/SVM/SimpleSVM.h>" namespace "OpenMS::SimpleSVM":
    
    cdef cppclass SVMPrediction "OpenMS::SimpleSVM::Prediction":

        SVMPrediction() nogil except +
        SVMPrediction(SVMPrediction) nogil except + #wrap-ignore

        Int label
        libcpp_map[int, double ] probabilities

