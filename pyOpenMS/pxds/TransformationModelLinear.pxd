from Param cimport *

# class TransformationModelLinear declared in TransformationModel.pxd

# keep TransformationModelInterpolated and TransformationModelLinear
# in separate files. Else autowrap can not
# distinguish the getDefaultParameters() static methods

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>" namespace "OpenMS::TransformationModelLinear":

    void getDefaultParameters(Param &) # wrap-attach:TransformationModelLinear
