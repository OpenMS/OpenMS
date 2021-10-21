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

        TransformationModel() nogil except +
        # private 
        TransformationModel(TransformationModel, Param) nogil except + # wrap-ignore

        Param getParameters() nogil except +

        # double evaluate(double value) nogil except +
        # void getDefaultParameters(Param & params) nogil except +
        
        void weightData(libcpp_vector[TM_DataPoint]& data) nogil except + # wrap-doc:Weight the data by the given weight function
        bool checkValidWeight(const String& weight, libcpp_vector[String]& valid_weights) nogil except + # wrap-doc:Check for a valid weighting function string
        double weightDatum(double& datum, const String& weight) nogil except + # wrap-doc:Weight the data according to the weighting function
        double unWeightDatum(double& datum, const String& weight) nogil except + # wrap-doc:Apply the reverse of the weighting function to the data
        libcpp_vector[ String ] getValidXWeights() nogil except + # wrap-doc:Returns a list of valid x weight function stringss
        libcpp_vector[ String ] getValidYWeights() nogil except + # wrap-doc:Returns a list of valid y weight function strings

        void unWeightData(libcpp_vector[TM_DataPoint] & data) nogil except + # wrap-doc:Unweight the data by the given weight function
        double checkDatumRange(const double & datum, const double & datum_min, const double & datum_max) nogil except + # wrap-doc:Check that the datum is within the valid min and max bounds
