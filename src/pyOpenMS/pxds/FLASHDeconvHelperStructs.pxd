from libcpp.vector cimport vector as libcpp_vector



cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>" namespace "OpenMS":

    cdef cppclass MassFeature_FDHS "OpenMS::FLASHDeconvHelperStructs::MassFeature":

        # wrap-inherits:

        # default constructor
        MassFeature_FDHS() except + nogil
        # copy constructor
        MassFeature_FDHS(MassFeature_FDHS &) except + nogil


