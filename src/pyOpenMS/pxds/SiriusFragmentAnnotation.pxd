from Types cimport *
from String cimport *
from MSSpectrum cimport *


cdef extern from "<OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>" namespace "OpenMS":

    cdef cppclass SiriusFragmentAnnotation:
        SiriusFragmentAnnotation() except + nogil 
        SiriusFragmentAnnotation(SiriusFragmentAnnotation &) except + nogil  # compiler

        libcpp_vector[MSSpectrum] extractAnnotationsFromSiriusFile(
                                                    String& path_to_sirius_workspace,
                                                    Size max_rank,
                                                    bool decoy,
                                                    bool use_exact_mass) except + nogil 
        libcpp_vector[MSSpectrum] extractSiriusAnnotationsTgtOnly(
                                                    libcpp_vector[String]& sirius_workspace_subdirs,
                                                    double score_threshold,
                                                    bool use_exact_mass,
                                                    bool resolve) except + nogil 
        libcpp_vector[ SiriusFragmentAnnotation_SiriusTargetDecoySpectra ] extractAndResolveSiriusAnnotations(libcpp_vector[ String ]& sirius_workspace_subdirs,
                                                                                                              double score_threshold,
                                                                                                              bool use_exact_mass) except + nogil 

    cdef cppclass SiriusFragmentAnnotation_SiriusTargetDecoySpectra "OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra":

      SiriusFragmentAnnotation_SiriusTargetDecoySpectra() except + nogil 
      SiriusFragmentAnnotation_SiriusTargetDecoySpectra(SiriusFragmentAnnotation_SiriusTargetDecoySpectra &) except + nogil  # compiler
