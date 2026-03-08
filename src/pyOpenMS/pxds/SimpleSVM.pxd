from Types cimport *
from String cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/ML/SVM/SimpleSVM.h>" namespace "OpenMS":

    cdef cppclass SimpleSVM(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        SimpleSVM() except + nogil  # wrap-doc:Simple interface to support vector machines for classification and regression via LIBSVM
        SimpleSVM(SimpleSVM &) except + nogil  # compiler

        void setup(libcpp_map[String, libcpp_vector[double]] & predictors,
                   libcpp_map[size_t, double] & outcomes,
                   bool classification) except + nogil
            # wrap-doc:
            #  Train an SVM/SVR model on the given data.
            #
            #  Predictor values are scaled in-place to [0, 1].
            #  Use outcomes to specify which observations are in the training set.
            #
            #  :param predictors: Mapping from predictor name to per-observation values
            #  :param outcomes: Mapping from observation index to class label or regression value
            #  :param classification: True for SVM classification, False for SVR regression

        void predict(libcpp_vector[SimpleSVM_Prediction] & predictions,
                     libcpp_vector[size_t] indexes) except + nogil
            # wrap-doc:
            #  Predict class labels (or regression values) for observations by index.
            #
            #  :param predictions: Output prediction results, in the same order as indexes
            #  :param indexes: Observation indexes to predict; if empty, predicts all observations

        void predict(libcpp_map[String, libcpp_vector[double]] & predictors,
                     libcpp_vector[SimpleSVM_Prediction] & predictions) except + nogil
            # wrap-doc:
            #  Predict class labels (or regression values) for novel data.
            #
            #  Predictor values are scaled using the ranges from training.
            #
            #  :param predictors: Mapping from predictor name to per-observation values
            #  :param predictions: Output prediction results

        void getFeatureWeights(libcpp_map[String, double] & feature_weights) except + nogil
            # wrap-doc:
            #  Get the weights of predictors in the SVM model.
            #
            #  Only supported for two-class classification with a linear kernel.
            #
            #  :param feature_weights: Output mapping from predictor name to weight

        void saveModel(const String & path) except + nogil
            # wrap-doc:
            #  Save the trained SVM model to a file.
            #
            #  :param path: Output file path
            #  :raises:
            #    Exception: Precondition if no model has been trained
            #    Exception: IOException if the file cannot be written

        void loadModel(const String & path) except + nogil
            # wrap-doc:
            #  Load a pre-trained SVM model from a file.
            #
            #  :param path: Input file path
            #  :raises:
            #    Exception: IOException if the file cannot be read

        double predictSingle(libcpp_vector[double] & scaled_feature_values) except + nogil
            # wrap-doc:
            #  Predict a single observation from pre-scaled feature values.
            #
            #  Does not use probability estimation. setup() is not required;
            #  loadModel() is sufficient.
            #
            #  :param scaled_feature_values: Pre-scaled feature values for one observation
            #  :returns: Predicted class label or regression value
            #  :raises:
            #    Exception: Precondition if no model has been trained or loaded

        void writeXvalResults(const String & path) except + nogil
            # wrap-doc:
            #  Write cross-validation parameter optimization results to a CSV file.
            #
            #  :param path: Output file path

        libcpp_map[String, libcpp_pair[double, double]] getScaling() except + nogil
            # wrap-doc:
            #  Get the data range of each predictor before scaling to [0, 1].

cdef extern from "<OpenMS/ML/SVM/SimpleSVM.h>" namespace "OpenMS::SimpleSVM":

    cdef cppclass SimpleSVM_Prediction "OpenMS::SimpleSVM::Prediction":
        SimpleSVM_Prediction() except + nogil
        SimpleSVM_Prediction(SimpleSVM_Prediction &) except + nogil  # compiler
        double outcome
        libcpp_map[double, double] probabilities
