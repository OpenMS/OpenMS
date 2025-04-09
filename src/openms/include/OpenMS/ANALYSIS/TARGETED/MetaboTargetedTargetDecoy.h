// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>


namespace OpenMS
{

  /**
  @brief Resolve overlapping fragments and missing decoys for experimental specific decoy generation
   in targeted/pseudo targeted metabolomics.
  */

  class OPENMS_DLLAPI MetaboTargetedTargetDecoy
  {
    public:
    /**
     @brief MetaboTargetDecoyMassMapping introduces a mapping of target and decoy masses and their respective compound reference
     using an identifier
    */
    class MetaboTargetDecoyMassMapping
    {
    public:

      String identifier; ///> unique identifier (e.g. m_id)
      String target_compound_ref; ///> identifier which allows to reference back to the target (e.g. target_transitions_id)
      String decoy_compound_ref; ///> identifier which allows to reference back to the decoy (e.g. decoy_transitions_id)
      std::vector<double> target_product_masses; ///> masses of target transitions
      std::vector<double> decoy_product_masses; ///> masses of decoy transitions
    };

    /**
    @brief Constructs a mass mapping of targets and decoys using the unique m_id identifier.

    @param t_exp TransitionExperiment holds compound and transition information used for the mapping.

    */
    static std::vector<MetaboTargetDecoyMassMapping> constructTargetDecoyMassMapping(const TargetedExperiment& t_exp);

    /**
    @brief Resolves overlapping target and decoy transition masses by adding a specifiable mass (e.g. CH2) to the overlapping decoy fragment.

    @param t_exp TransitionExperiment holds compound and transition information
    @param mappings map of identifier to target and decoy masses
    @param mass_to_add (e.g. CH2)
    @param mz_tol m/z tolerarance for target and decoy transition masses to be considered overlapping
    @param mz_tol_unit m/z tolerance unit ("ppm" or "Da")
    */
    static void resolveOverlappingTargetDecoyMassesByDecoyMassShift(TargetedExperiment& t_exp, std::vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping>& mappings, const double& mass_to_add, const double& mz_tol, const String& mz_tol_unit);

    /**
    @brief Generate a decoy for targets where fragmentation tree re-rooting was not possible, by adding a specifiable mass to the target fragments.

    @param t_exp TransitionExperiment holds compound and transition information
    @param mappings map of identifier to target and decoy masses
    @param mass_to_add the maximum number of transitions required per assay

    */
    static void generateMissingDecoysByMassShift(TargetedExperiment& t_exp, std::vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping>& mappings, const double& mass_to_add);

  protected:
    /**
    @brief Generate a TransitionMap based on Compound_Ref and ReactionMonitoringTransitions

    @param t_exp TransitionExperiment holds compound and transition information
    */
    static std::map<String, std::vector<OpenMS::ReactionMonitoringTransition> > constructTransitionsMap_(const TargetedExperiment& t_exp);
  };

} // namespace OpenMS
