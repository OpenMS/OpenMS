from String cimport *
from Param cimport *

cdef extern from "<OpenMS/APPLICATIONS/MapAlignerBase.h>" namespace "OpenMS":

    cdef cppclass MapAlignerBase:
        pass

cdef extern from "<OpenMS/APPLICATIONS/MapAlignerBase.h>" namespace "OpenMS::MapAlignerBase":

    Param getModelDefaults(const String& default_model) except + nogil  # wrap-attach:MapAlignerBase
