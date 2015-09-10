from Types cimport *
from Param cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>" namespace "OpenMS":

    cdef cppclass TransformationModelLinear(TransformationModel):
        # wrap-inherits:
        #   TransformationModel

        void getDefaultParameters(Param &)

        # TransformationModelLinear(DataPoints & data, Param & params) nogil except +
        # double evaluate(double value) nogil except +
        # void getParameters(double & slope, double & intercept) nogil except +
        void invert() nogil except +

