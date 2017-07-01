from Types cimport *
from Param cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>" namespace "OpenMS":

    cdef cppclass TransformationModel:
        # wrap-ignore

        TransformationModel()  nogil except +
        TransformationModel(TransformationModel, Param) nogil except + # wrap-ignore

        Param getParameters() nogil except +

        # double evaluate(double value) nogil except +
        # void getDefaultParameters(Param & params) nogil except +
        
        virtual void weightData(libcpp_vector[ libcpp_pair[double, double] ]& data, Param & params) nogil except +
        virtual void weightData(libcpp_vector[ libcpp_pair[double, double] ]& data, Param & params) nogil except +
        bool checkValidWeight(libcpp_string& weight, libcpp_vector[libcpp_string] ]& valid_weights) nogil except +
        double weightDatum(double& datum, libcpp_string& weight) nogil except +
        double unWeightDatum(double& datum, libcpp_string& weight) nogil except +
        libcpp_string getValidXWeights() nogil except +
        libcpp_string getValidYWeights() nogil except +

