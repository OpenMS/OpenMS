from Types cimport *
from TargetedExperiment cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h>" namespace "OpenMS":

  cdef cppclass MetaboTargetedTargetDecoy "OpenMS::MetaboTargetedTargetDecoy":

      MetaboTargetedTargetDecoy() nogil except +
      MetaboTargetedTargetDecoy(MetaboTargetedTargetDecoy) nogil except + #wrap-ignore

      libcpp_vector[ MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] constructTargetDecoyMassMapping(TargetedExperiment& t_exp) nogil except +

      void resolveOverlappingTargetDecoyMassesByIndividualMassShift(TargetedExperiment& t_exp, libcpp_vector[ MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping ]& mappings, double& mass_to_add) nogil except +

      void generateMissingDecoysByMassShift(TargetedExperiment& t_exp, libcpp_vector[ MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping ]& mappings, double& mass_to_add) nogil except +

  cdef cppclass MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping "OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping":

      MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping() nogil except +
      MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping(MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping) nogil except +
