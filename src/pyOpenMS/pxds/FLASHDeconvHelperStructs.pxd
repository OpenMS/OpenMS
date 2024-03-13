from libcpp.vector cimport vector as libcpp_vector
from String cimport *



cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>" namespace "OpenMS":

    cdef cppclass MassFeature_FDHS "OpenMS::FLASHDeconvHelperStructs::MassFeature":

        # wrap-inherits:

        # default constructor
        MassFeature_FDHS() except + nogil
        # copy constructor
        MassFeature_FDHS(MassFeature_FDHS &) except + nogil

    cdef cppclass Tag "OpenMS::FLASHDeconvHelperStructs::Tag":

        # wrap-inherits:

        # default constructor
        Tag(String seq, double n_mass, double c_mass, libcpp_vector[int]& scores, libcpp_vector[double] & mzs) except + nogil
        # copy constructor
        Tag(Tag &) except + nogil

        String getSequence() except + nogil
        libcpp_vector[double] getMzs() except + nogil
        double getNtermMass() except + nogil
        double getCtermMass() except + nogil
        int getScore() except + nogil
        int getScore(int pos) except + nogil


