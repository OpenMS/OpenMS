from libcpp.vector cimport vector as libcpp_vector
from String cimport *



cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>" namespace "OpenMS":

    cdef cppclass MassFeature_FDHS "OpenMS::FLASHDeconvHelperStructs::MassFeature":

        # wrap-inherits:

        # default constructor
        MassFeature_FDHS() except + nogil
        # copy constructor
        MassFeature_FDHS(MassFeature_FDHS &) except + nogil

    cdef cppclass Tag_FDHS "OpenMS::FLASHDeconvHelperStructs::Tag":

        # wrap-inherits:

        # default constructor
        Tag_FDHS(String seq, double n_mass, double c_mass, int charge, double score, libcpp_vector[double] & mzs) except + nogil
        # copy constructor
        Tag_FDHS(Tag_FDHS &) except + nogil




