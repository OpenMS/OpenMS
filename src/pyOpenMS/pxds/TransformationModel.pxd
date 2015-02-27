from Types cimport *
from Param cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>" namespace "OpenMS":

    cdef cppclass TransformationModel:
        # wrap-ignore

        TransformationModel()  nogil except +
        TransformationModel(TransformationModel, Param) nogil except + # wrap-ignore

        Param getParameters() nogil except +

        # double evaluate(double value) nogil except +
        # void getDefaultParameters(Param & params) nogil except +

