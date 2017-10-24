from Types cimport *
from Param cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>" namespace "OpenMS":


    cdef cppclass TransformationModelLinear(TransformationModel):
        # wrap-inherits:
        #   TransformationModel

        TransformationModelLinear(TransformationModelLinear) nogil except + #wrap-ignore
        TransformationModelLinear(libcpp_vector[TM_DataPoint]& data, Param& params) nogil except +

        void getDefaultParameters(Param&)

        double evaluate(double value) nogil except +
        # void getParameters(double & slope, double & intercept, libcpp_string& x_weight, libcpp_string& y_weight, double & x_datum_min, double & x_datum_max, double & y_datum_min, double & y_datum_max) nogil except +
        void invert() nogil except +

