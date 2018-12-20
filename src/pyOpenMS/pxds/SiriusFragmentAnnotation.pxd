from Types cimport *
from String cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>" namespace "OpenMS":

    cdef cppclass SiriusFragmentAnnotation:
        SiriusFragmentAnnotation() nogil except +
        SiriusFragmentAnnotation(SiriusFragmentAnnotation) nogil except +
                
        void  extractSiriusFragmentAnnotationMapping(String& path_to_sirius_workspace,
                                                     MSSpectrum& msspectrum_to_fill,
                                                     bool use_exact_mass) nogil except +