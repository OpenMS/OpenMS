from libcpp.string cimport string as libcpp_string
from Types cimport *
from ProgressLogger cimport *
from String cimport *
from TextFile cimport *
from File cimport *

cdef extern from "<OpenMS/ANALYSIS/SVM/SVMWrapper.h>" namespace "OpenMS":
    
    cdef cppclass SVMWrapper "OpenMS::SVMWrapper":
        SVMWrapper() nogil except +
        SVMWrapper(SVMWrapper) nogil except + #wrap-ignore

        void setParameter(SVM_parameter_type type_, Int value) nogil except +
        void setParameter(SVM_parameter_type type_, double value) nogil except +
        # Int train(struct svm_problem *problem) nogil except +
        Int train(SVMData &problem) nogil except +
        void saveModel(libcpp_string modelFilename) nogil except +
        void loadModel(libcpp_string modelFilename) nogil except +
        # void predict(struct svm_problem *problem, libcpp_vector[ double ] &predicted_labels) nogil except +
        void predict(SVMData &problem, libcpp_vector[ double ] &results) nogil except +
        Int getIntParameter(SVM_parameter_type type_) nogil except +
        double getDoubleParameter(SVM_parameter_type type_) nogil except +
        # TODO STL map with wrapped key
        # void predict(libcpp_vector[ svm_node * ] &vectors, libcpp_vector[ double ] &predicted_rts) nogil except +
        # double performCrossValidation(svm_problem *problem_ul, SVMData &problem_l, bool is_labeled, libcpp_map[ SVM_parameter_type, double ] &start_values_map, libcpp_map[ SVM_parameter_type, double ] &step_sizes_map, libcpp_map[ SVM_parameter_type, double ] &end_values_map, Size number_of_partitions, Size number_of_runs, libcpp_map[ SVM_parameter_type, double ] &best_parameters, bool additive_step_sizesrue, bool outputalse, String, bool mcc_as_performance_measurealse) nogil except +
        double getSVRProbability() nogil except +
        # void getSignificanceBorders(svm_problem *data, libcpp_pair[ double, double ] &borders, double confidence.95, Size number_of_runs, Size number_of_partitions, double step_size.01, Size max_iterations000000) nogil except +
        void getSignificanceBorders(SVMData &data, libcpp_pair[ double, double ] &sigmas, double confidence, Size number_of_runs, Size number_of_partitions, double step_size, Size max_iterations) nogil except +
        double getPValue(double sigma1, double sigma2, libcpp_pair[ double, double ] point) nogil except +
        # void getDecisionValues(svm_problem *data, libcpp_vector[ double ] &decision_values) nogil except +
        # void scaleData(svm_problem *data, Int max_scale_value1) nogil except +
        # svm_problem * computeKernelMatrix(svm_problem *problem1, svm_problem *problem2) nogil except +
        # svm_problem * computeKernelMatrix(SVMData &problem1, SVMData &problem2) nogil except +
        # void setTrainingSample(svm_problem *training_sample) nogil except +
        void setTrainingSample(SVMData &training_sample) nogil except +
        # void getSVCProbabilities(struct svm_problem *problem, libcpp_vector[ double ] &probabilities, libcpp_vector[ double ] &prediction_labels) nogil except +
        void setWeights(libcpp_vector[ int ] &weight_labels, libcpp_vector[ double ] &weights) nogil except +
        # void createRandomPartitions(svm_problem *problem, Size number, libcpp_vector[ svm_problem * ] &partitions) nogil except +
        void createRandomPartitions(SVMData &problem, Size number, libcpp_vector[ SVMData ] &problems) nogil except +
        # TODO: Mismatch between C++ return type ([u'svm_problem *']) and Python return type (['void']) in function public mergePartitions:
        # svm_problem * mergePartitions(libcpp_vector[ svm_problem * ] &problems, Size except) nogil except +
        void mergePartitions(libcpp_vector[ SVMData ] &problems, Size except_, SVMData &merged_problem) nogil except +
        # void getLabels(svm_problem *problem, libcpp_vector[ double ] &labels) nogil except +
        # double kernelOligo(libcpp_vector[ libcpp_pair[ int, double ] ] &x, libcpp_vector[ libcpp_pair[ int, double ] ] &y, libcpp_vector[ double ] &gauss_table, int max_distance1) nogil except +
        # double kernelOligo(svm_node *x, svm_node *y, libcpp_vector[ double ] &gauss_table, double sigma_square, Size max_distance0) nogil except +
        void calculateGaussTable(Size border_length, double sigma, libcpp_vector[ double ] &gauss_table) nogil except +

cdef extern from "<OpenMS/ANALYSIS/SVM/SVMWrapper.h>" namespace "OpenMS":
    
    cdef cppclass SVMData "OpenMS::SVMData":
        SVMData() nogil except +
        SVMData(SVMData) nogil except + #wrap-ignore

        # TODO nested STL
        # libcpp_vector[ libcpp_vector[ libcpp_pair[ Int, double ] ] ] sequences
        libcpp_vector[ double ] labels
        # TODO nested STL
        # SVMData(libcpp_vector[ libcpp_vector[ libcpp_pair[ Int, double ] ] ] &seqs, libcpp_vector[ double ] &lbls) nogil except +
        bool operator==(SVMData &rhs) nogil except +
        bool store(const String &filename) nogil except +
        bool load(const String &filename) nogil except +


cdef extern from "<OpenMS/ANALYSIS/SVM/SVMWrapper.h>" namespace "OpenMS::SVMWrapper":
    cdef enum SVM_parameter_type "OpenMS::SVMWrapper::SVM_parameter_type":
        #wrap-attach:
        #    SVMWrapper
        SVM_TYPE
        KERNEL_TYPE
        DEGREE
        C
        NU
        P
        GAMMA
        PROBABILITY
        SIGMA
        BORDER_LENGTH

cdef extern from "<OpenMS/ANALYSIS/SVM/SVMWrapper.h>" namespace "OpenMS::SVMWrapper":
    cdef enum SVM_kernel_type "OpenMS::SVMWrapper::SVM_kernel_type":
        #wrap-attach:
        #    SVMWrapper
        OLIGO
        OLIGO_COMBINED

