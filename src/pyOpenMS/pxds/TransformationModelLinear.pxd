from Types cimport *
from Param cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

from TransformationModel cimport *

# ctypedef libcpp_vector[ libcpp_pair[double, double] ] DataPoints;

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>" namespace "OpenMS":


    cdef cppclass TransformationModelLinear(TransformationModel):
        # wrap-inherits:
        #   TransformationModel

        TransformationModelLinear(TransformationModelLinear) nogil except + #wrap-ignore
        TransformationModelLinear(libcpp_vector[ libcpp_pair[double, double] ]& data, Param & params) nogil except +

        void getDefaultParameters(Param &)
        
        double evaluate(double value) nogil except +
        void getParameters(double & slope, double & intercept) nogil except +
        void invert() nogil except +

