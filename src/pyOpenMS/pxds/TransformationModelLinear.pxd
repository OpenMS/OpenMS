from Types cimport *
from Param cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

from TransformationModel cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>" namespace "OpenMS":

    cdef cppclass TransformationModelLinear(TransformationModel):
        # wrap-inherits:
        #   TransformationModel
        pass
        void getDefaultParameters(Param &)
