from Types cimport *
from Param cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>" namespace "OpenMS":

    cdef cppclass TM_DataPoint "OpenMS::TransformationModel::DataPoint":
        TM_DataPoint() nogil except +
        TM_DataPoint(double, double) nogil except +
        TM_DataPoint(double, double, const String&) nogil except +
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
        
        void weightData(libcpp_vector[TM_DataPoint]& data) nogil except +
        void weightData(libcpp_vector[TM_DataPoint]& data) nogil except +
        bool checkValidWeight(const String& weight, libcpp_vector[String]& valid_weights) nogil except +
        double weightDatum(double& datum, const String& weight) nogil except +
        double unWeightDatum(double& datum, const String& weight) nogil except +
        libcpp_vector[ String ] getValidXWeights() nogil except +
        libcpp_vector[ String ] getValidYWeights() nogil except +

        void unWeightData(libcpp_vector[TM_DataPoint] & data) nogil except +
        double checkDatumRange(const double & datum, const double & datum_min, const double & datum_max) nogil except +

