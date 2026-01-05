from Types cimport *
from TransformationModel cimport *
from BSpline2d cimport *

# ctypedef libcpp_vector[ libcpp_pair[double, double] ] DataPoints;

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>" namespace "OpenMS":

    cdef cppclass TransformationModelBSpline(TransformationModel) :
        # wrap-inherits:
        #  TransformationModel

        # copy constructor of 'TransformationModelBSpline' is implicitly deleted because base class 'OpenMS::TransformationModel' has an inaccessible copy constructor public TransformationModel
        TransformationModelBSpline(TransformationModelBSpline) except + nogil  # wrap-ignore

        TransformationModelBSpline(libcpp_vector[TM_DataPoint]& data, Param& params) except + nogil
        double evaluate(double value) except + nogil  # wrap-doc:Evaluates the model at the given values

        @staticmethod
        void getDefaultParameters(Param& params) except + nogil  # wrap-doc:Gets the default parameters
