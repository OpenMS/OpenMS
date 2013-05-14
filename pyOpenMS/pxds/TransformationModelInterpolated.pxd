from Param cimport *

# class TransformationModelInterpolated declared in TransformationModel.pxd

# keep TransformationModelInterpolated, TransformationModelLinear and
# TransformationModelBSpline in separate files. Else autowrap can not
# distinguish the getDefaultParameters() static methods

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>" namespace "OpenMS::TransformationModelInterpolated":

    void getDefaultParameters(Param &) # wrap-attach:TransformationModelInterpolated
