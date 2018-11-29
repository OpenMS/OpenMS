from Types cimport *
from String cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/FragmentAnnotation.h>" namespace "OpenMS":

    cdef cppclass FragmentAnnotation:
        FragmentAnnotation() nogil except +
        FragmentAnnotation(FragmentAnnotation) nogil except +
                
        void extractFragmentAnnotationMapping(String& path_to_sirius_workspace,
                                              MSSpectrum& msspectrum_to_fill,
                                              bool use_exact_mass) nogil except +