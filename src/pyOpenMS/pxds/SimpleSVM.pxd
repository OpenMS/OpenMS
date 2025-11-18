from libcpp.map cimport map as libcpp_map
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from DefaultParamHandler cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/ML/SVM/SimpleSVM.h>" namespace "OpenMS":

    ctypedef libcpp_map[String, libcpp_vector[double]] SimpleSVM_PredictorMap "OpenMS::SimpleSVM::PredictorMap"
    ctypedef libcpp_map[String, libcpp_pair[double, double]] SimpleSVM_ScaleMap "OpenMS::SimpleSVM::ScaleMap"

    cdef cppclass SimpleSVM_Prediction "OpenMS::SimpleSVM::Prediction":
        SimpleSVM_Prediction() except + nogil
        SimpleSVM_Prediction(SimpleSVM_Prediction&) except + nogil
        double outcome
        libcpp_map[double, double] probabilities

    cdef cppclass SimpleSVM(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        SimpleSVM() except + nogil
        SimpleSVM(SimpleSVM&) except + nogil
        void setup(SimpleSVM_PredictorMap& predictors,
                   const libcpp_map[Size, double]& outcomes,
                   bool classification) except + nogil
        void predict(libcpp_vector[SimpleSVM_Prediction]& predictions,
                     libcpp_vector[Size] indexes) const except + nogil
        void predict(SimpleSVM_PredictorMap& predictors,
                     libcpp_vector[SimpleSVM_Prediction]& predictions) const except + nogil
        void getFeatureWeights(libcpp_map[String, double]& feature_weights) const except + nogil
        void writeXvalResults(const String& path) const except + nogil
        SimpleSVM_ScaleMap getScaling() const except + nogil
