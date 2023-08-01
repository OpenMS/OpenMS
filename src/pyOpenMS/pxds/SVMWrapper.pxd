from Types cimport *
from ProgressLogger cimport *
from String cimport *
from TextFile cimport *
from File cimport *

cdef extern from "<OpenMS/ANALYSIS/SVM/SVMWrapper.h>" namespace "OpenMS":
    
    cdef cppclass SVMWrapper "OpenMS::SVMWrapper":
        SVMWrapper() except + nogil 
        SVMWrapper(SVMWrapper &) except + nogil  # compiler

        void setParameter(SVM_parameter_type type_, Int value) except + nogil 
        void setParameter(SVM_parameter_type type_, double value) except + nogil 
        # Int train(struct svm_problem *problem) except + nogil 
        Int train(SVMData &problem) except + nogil  # wrap-doc:The svm is trained with the data stored in the 'svm_problem' structure
        void saveModel(String modelFilename) except + nogil  # wrap-doc:The model of the trained svm is saved into 'modelFilename'
        void loadModel(String modelFilename) except + nogil  # wrap-doc:The svm-model is loaded. After this, the svm is ready for prediction
        # void predict(struct svm_problem *problem, libcpp_vector[ double ] &predicted_labels) except + nogil 
        void predict(SVMData &problem, libcpp_vector[ double ] &results) except + nogil  # wrap-doc:The prediction process is started and the results are stored in 'predicted_labels'
        Int getIntParameter(SVM_parameter_type type_) except + nogil 
        double getDoubleParameter(SVM_parameter_type type_) except + nogil 
        # TODO STL map with wrapped key
        # void predict(libcpp_vector[ svm_node * ] &vectors, libcpp_vector[ double ] &predicted_rts) except + nogil 
        # double performCrossValidation(svm_problem *problem_ul, SVMData &problem_l, bool is_labeled, libcpp_map[ SVM_parameter_type, double ] &start_values_map, libcpp_map[ SVM_parameter_type, double ] &step_sizes_map, libcpp_map[ SVM_parameter_type, double ] &end_values_map, Size number_of_partitions, Size number_of_runs, libcpp_map[ SVM_parameter_type, double ] &best_parameters, bool additive_step_sizesrue, bool outputalse, String, bool mcc_as_performance_measurealse) except + nogil 
        double getSVRProbability() except + nogil 
        # void getSignificanceBorders(svm_problem *data, libcpp_pair[ double, double ] &borders, double confidence.95, Size number_of_runs, Size number_of_partitions, double step_size.01, Size max_iterations000000) except + nogil 
        void getSignificanceBorders(SVMData &data, libcpp_pair[ double, double ] &sigmas, double confidence, Size number_of_runs, Size number_of_partitions, double step_size, Size max_iterations) except + nogil 
        double getPValue(double sigma1, double sigma2, libcpp_pair[ double, double ] point) except + nogil 
        # void getDecisionValues(svm_problem *data, libcpp_vector[ double ] &decision_values) except + nogil 
        # void scaleData(svm_problem *data, Int max_scale_value1) except + nogil 
        # svm_problem * computeKernelMatrix(svm_problem *problem1, svm_problem *problem2) except + nogil 
        # svm_problem * computeKernelMatrix(SVMData &problem1, SVMData &problem2) except + nogil 
        # void setTrainingSample(svm_problem *training_sample) except + nogil 
        void setTrainingSample(SVMData &training_sample) except + nogil 
        # void getSVCProbabilities(struct svm_problem *problem, libcpp_vector[ double ] &probabilities, libcpp_vector[ double ] &prediction_labels) except + nogil 
        void setWeights(libcpp_vector[ int ] &weight_labels, libcpp_vector[ double ] &weights) except + nogil 
        # void createRandomPartitions(svm_problem *problem, Size number, libcpp_vector[ svm_problem * ] &partitions) except + nogil 
        void createRandomPartitions(SVMData &problem, Size number, libcpp_vector[ SVMData ] &problems) except + nogil 
        # TODO: Mismatch between C++ return type ([u'svm_problem *']) and Python return type (['void']) in function public mergePartitions:
        # svm_problem * mergePartitions(libcpp_vector[ svm_problem * ] &problems, Size except) except + nogil 
        void mergePartitions(libcpp_vector[ SVMData ] &problems, Size except_, SVMData &merged_problem) except + nogil 
        # void getLabels(svm_problem *problem, libcpp_vector[ double ] &labels) except + nogil 
        # double kernelOligo(libcpp_vector[ libcpp_pair[ int, double ] ] &x, libcpp_vector[ libcpp_pair[ int, double ] ] &y, libcpp_vector[ double ] &gauss_table, int max_distance1) except + nogil 
        # double kernelOligo(svm_node *x, svm_node *y, libcpp_vector[ double ] &gauss_table, double sigma_square, Size max_distance0) except + nogil 
        void calculateGaussTable(Size border_length, double sigma, libcpp_vector[ double ] &gauss_table) except + nogil 

cdef extern from "<OpenMS/ANALYSIS/SVM/SVMWrapper.h>" namespace "OpenMS":
    
    cdef cppclass SVMData "OpenMS::SVMData":
        SVMData() except + nogil 
        SVMData(SVMData &) except + nogil  # compiler

        # TODO nested STL
        # libcpp_vector[ libcpp_vector[ libcpp_pair[ Int, double ] ] ] sequences
        libcpp_vector[ double ] labels
        # TODO nested STL
        # SVMData(libcpp_vector[ libcpp_vector[ libcpp_pair[ Int, double ] ] ] &seqs, libcpp_vector[ double ] &lbls) except + nogil 
        bool operator==(SVMData &rhs) except + nogil 
        bool store(const String &filename) except + nogil 
        bool load(const String &filename) except + nogil 


cdef extern from "<OpenMS/ANALYSIS/SVM/SVMWrapper.h>" namespace "OpenMS::SVMWrapper":
    cdef enum SVM_parameter_type "OpenMS::SVMWrapper::SVM_parameter_type":
        #wrap-attach:
        #   SVMWrapper
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
        #   SVMWrapper
        OLIGO
        OLIGO_COMBINED

