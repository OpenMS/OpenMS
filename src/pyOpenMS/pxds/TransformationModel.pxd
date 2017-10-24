from Types cimport *
from Param cimport *
from String cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as String

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>" namespace "OpenMS":

    cdef cppclass TM_DataPoint "OpenMS::TransformationModel::DataPoint":
        TM_DataPoint() nogil except +
        TM_DataPoint(double, double) nogil except +
        TM_DataPoint(double, double, String&) nogil except +
        TM_DataPoint(TM_DataPoint) nogil except + #wrap-ignore
        bool operator<(TM_DataPoint&) nogil except +
        bool operator==(TM_DataPoint&) nogil except +

        double first
        double second
        String note


    cdef cppclass TransformationModel:
        # wrap-ignore
        # no-pxd-import

        TransformationModel()  nogil except +
        TransformationModel(TransformationModel, Param) nogil except + # wrap-ignore

        Param getParameters() nogil except +

        # double evaluate(double value) nogil except +
        # void getDefaultParameters(Param & params) nogil except +
        
        void weightData(libcpp_vector[ libcpp_pair[double, double] ]& data, Param & params) nogil except +
        void weightData(libcpp_vector[ libcpp_pair[double, double] ]& data, Param & params) nogil except +
        bool checkValidWeight(String& weight, libcpp_vector[String]& valid_weights) nogil except +
        double weightDatum(double& datum, String& weight) nogil except +
        double unWeightDatum(double& datum, String& weight) nogil except +
        # String getValidXWeights() nogil except +
        # String getValidYWeights() nogil except +
