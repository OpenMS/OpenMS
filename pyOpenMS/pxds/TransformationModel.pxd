from Types cimport *
from Param cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>" namespace "OpenMS":


    cdef cppclass TransformationModel:
        # wrap-ignore
        
        TransformationModel()  nogil except +
        TransformationModel(TransformationModel, Param) nogil except + # wrap-ignore
        void getParameters(Param &)

    cdef cppclass TransformationModelLinear(TransformationModel):
        # wrap-inherits:
        #   TransformationModel
        pass

    cdef cppclass TransformationModelBSpline(TransformationModel):
        # wrap-inherits:
        #   TransformationModel
        pass
        
    cdef cppclass TransformationModelInterpolated(TransformationModel):
        # wrap-inherits:
        #   TransformationModel
        pass
        
