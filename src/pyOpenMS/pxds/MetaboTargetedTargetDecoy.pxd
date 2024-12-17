from Types cimport *
from TargetedExperiment cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h>" namespace "OpenMS":

  cdef cppclass MetaboTargetedTargetDecoy "OpenMS::MetaboTargetedTargetDecoy":

      MetaboTargetedTargetDecoy() except + nogil  # wrap-doc:Resolve overlapping fragments and missing decoys for experimental specific decoy generation in targeted/pseudo targeted metabolomics
      MetaboTargetedTargetDecoy(MetaboTargetedTargetDecoy &) except + nogil  # compiler

      libcpp_vector[ MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] constructTargetDecoyMassMapping(TargetedExperiment& t_exp) except + nogil 
        # wrap-doc:
                #  Constructs a mass mapping of targets and decoys using the unique m_id identifier
                #  
                #  
                #  :param t_exp: TransitionExperiment holds compound and transition information used for the mapping

      void resolveOverlappingTargetDecoyMassesByDecoyMassShift(TargetedExperiment& t_exp, libcpp_vector[ MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping ]& mappings, double& mass_to_add, double& mz_tol, String& mz_tol_unit) except + nogil 
        # wrap-doc:
                #  Resolves overlapping target and decoy transition masses by adding a specifiable mass (e.g. CH2) to the overlapping decoy fragment
                #  
                #  
                #  :param t_exp: TransitionExperiment holds compound and transition information
                #  :param mappings: Map of identifier to target and decoy masses
                #  :param mass_to_add: (e.g. CH2)
                #  :param mz_tol: m/z tolerarance for target and decoy transition masses to be considered overlapping
                #  :param mz_tol_unit: m/z tolerance unit

      void generateMissingDecoysByMassShift(TargetedExperiment& t_exp, libcpp_vector[ MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping ]& mappings, double& mass_to_add) except + nogil 
        # wrap-doc:
                #  Generate a decoy for targets where fragmentation tree re-rooting was not possible, by adding a specifiable mass to the target fragments
                #  
                #  
                #  :param t_exp: TransitionExperiment holds compound and transition information
                #  :param mappings: Map of identifier to target and decoy masses
                #  :param mass_to_add: The maximum number of transitions required per assay

  cdef cppclass MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping "OpenMS::MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping":

      MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping() except + nogil 
      MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping(MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping &) except + nogil 
