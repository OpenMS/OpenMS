from Types cimport *
from TransformationModel cimport *
from BSpline2d cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>" namespace "OpenMS":
    
    cdef cppclass TransformationModelBSpline(TransformationModel) :
        # wrap-inherits:
        #  TransformationModel
        TransformationModelBSpline() nogil except +
        TransformationModelBSpline(TransformationModelBSpline) nogil except + #wrap-ignore
        TransformationModelBSpline(DataPoints & data, Param & params) nogil except +
        double evaluate(double value) nogil except +
        void getDefaultParameters(Param & params) nogil except +

