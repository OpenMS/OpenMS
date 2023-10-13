from Types cimport *
from Param cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>" namespace "OpenMS":

    cdef cppclass TM_DataPoint "OpenMS::TransformationModel::DataPoint":
        TM_DataPoint() except + nogil 
        TM_DataPoint(double, double) except + nogil 
        TM_DataPoint(double, double, const String&) except + nogil 
        TM_DataPoint(TM_DataPoint) except + nogil  #wrap-ignore
        bool operator<(TM_DataPoint&) except + nogil 
        bool operator==(TM_DataPoint&) except + nogil 

        double first
        double second
        String note


    cdef cppclass TransformationModel:
        # wrap-ignore
        # no-pxd-import

        TransformationModel() except + nogil 
        # private 
        TransformationModel(TransformationModel, Param) except + nogil  # wrap-ignore

        Param getParameters() except + nogil 

        # double evaluate(double value) except + nogil 
        # void getDefaultParameters(Param & params) except + nogil 
        
        void weightData(libcpp_vector[TM_DataPoint]& data) except + nogil  # wrap-doc:Weight the data by the given weight function
        bool checkValidWeight(const String& weight, libcpp_vector[String]& valid_weights) except + nogil  # wrap-doc:Check for a valid weighting function string
        double weightDatum(double& datum, const String& weight) except + nogil  # wrap-doc:Weight the data according to the weighting function
        double unWeightDatum(double& datum, const String& weight) except + nogil  # wrap-doc:Apply the reverse of the weighting function to the data
        libcpp_vector[ String ] getValidXWeights() except + nogil  # wrap-doc:Returns a list of valid x weight function stringss
        libcpp_vector[ String ] getValidYWeights() except + nogil  # wrap-doc:Returns a list of valid y weight function strings

        void unWeightData(libcpp_vector[TM_DataPoint] & data) except + nogil  # wrap-doc:Unweight the data by the given weight function
        double checkDatumRange(const double & datum, const double & datum_min, const double & datum_max) except + nogil  # wrap-doc:Check that the datum is within the valid min and max bounds
