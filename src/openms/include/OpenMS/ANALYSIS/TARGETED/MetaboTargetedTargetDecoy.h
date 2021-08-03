// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

    */
    static void resolveOverlappingTargetDecoyMassesByIndividualMassShift(TargetedExperiment& t_exp, std::vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping>& mappings, const double& mass_to_add);

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
