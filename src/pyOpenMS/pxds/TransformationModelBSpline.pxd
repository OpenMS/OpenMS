from Types cimport *
from TransformationModel cimport *
from BSpline2d cimport *

# ctypedef libcpp_vector[ libcpp_pair[double, double] ] DataPoints;

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>" namespace "OpenMS":

    cdef cppclass TransformationModelBSpline(TransformationModel) :
        # wrap-inherits:
        #  TransformationModel

        TransformationModelBSpline(TransformationModelBSpline) nogil except + #wrap-ignore

        TransformationModelBSpline(libcpp_vector[TM_DataPoint]& data, Param& params) nogil except +
        void getDefaultParameters(Param &) nogil except +
        double evaluate(double value) nogil except +

# COMMENT: wrap static methods
#cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>" namespace "OpenMS::TransformationModelBSpline":
#    void getDefaultParameters(Param& params) nogil except +  # wrap-attach:TransformationModelBSpline
