// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LABELING_BASELABELER_H
#define OPENMS_SIMULATION_LABELING_BASELABELER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{

  /**
  @brief Abstract base class for all kinds of labeling techniques
  */
  class OPENMS_DLLAPI BaseLabeler :
    public DefaultParamHandler
  {
public:

    /// constructor
    BaseLabeler();

    /// destructor
    ~BaseLabeler() override;

    /// register all derived classes here (implemented in file BaseLabeler_impl.h)
    static void registerChildren();

    /**
      @brief Returns the default parameters. Re-implement

      Re-implement if you derive a class and have to incorporate sub-algorithm default parameters.
    */
    virtual Param getDefaultParameters() const;

    /**
      @brief Set the random number generator

      Internally a pointer to the RNG is stored.

    */
    virtual void setRnd(SimTypes::MutableSimRandomNumberGeneratorPtr rng);

    /**
      @brief Checks the (simulation) params passed if they are consistent with
      the labeling technique.

      @param param Param object containing the simulation parameters
      @throws Exception::InvalidParameter if the given parameters are not consistent with the labeling technique
      */
    virtual void preCheck(Param& param) const = 0;

    /**
    @name Labeling Hooks
    */
    //@{

    /// Hook to prepare the simulation process
    virtual void setUpHook(SimTypes::FeatureMapSimVector& /* features */) = 0;

    /// Labeling between digestion and rt simulation
    virtual void postDigestHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) = 0;

    /// Labeling after rt simulation
    virtual void postRTHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) = 0;

    /// Labeling after detectability simulation
    virtual void postDetectabilityHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) = 0;

    /// Labeling after ionization
    virtual void postIonizationHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) = 0;

    /// Labeling after raw signal generation
    virtual void postRawMSHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) = 0;

    /// Labeling after Tandem MS (e.g. iTRAQ)
    virtual void postRawTandemMSHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */, SimTypes::MSSimExperiment& /* simulated map */) = 0;

    //@}

    ConsensusMap& getConsensus();

    /**
      @brief Get short description of the labeler (e.g., channels used)

      Used to add a short description to the labeling section within the INI file.
    */
    const String& getDescription() const;

    /**
      @brief to ensure standardized meta value names across labelers for channel intensity

      Use this function to get the name of the meta value which holds intensity for channel @p channel_index

    */
    String getChannelIntensityName(const Size channel_index) const;


protected:
    /**
      @brief Creates an empty FeatureMap with the merged ProteinIdentifications from
      all FeatureMaps contained in @p maps

      @param maps       Vector of FeatureMaps containing the features that will be merged
      @return           A FeatureMap containing all ProteinIdentifications of the input maps
      */
    SimTypes::FeatureMapSim mergeProteinIdentificationsMaps_(const SimTypes::FeatureMapSimVector& maps);

    /**
      @brief join all protein references of two features

      When merging peptides from different channels, the protein accessions should remain intact.
      Usually joining features is based on peptide sequence, so all protein hits should be valid.

      @param target
      @param source
    */
    void mergeProteinAccessions_(Feature& target, const Feature& source) const;

    /**
      @brief Based on the stored consensus recompute the associations for the passed features, assuming
             that the features where derived from the features stored in the consensus.

      @param simulated_features FeatureMap containing features derived from the ones, stored in the
                                consensus
    */
    void recomputeConsensus_(const SimTypes::FeatureMapSim& simulated_features);


    ConsensusMap consensus_;

    SimTypes::MutableSimRandomNumberGeneratorPtr rng_;

    String channel_description_;

  };
} // namespace OpenMS

#endif //#ifndef OPENMS_SIMULATION_LABELING_BASELABELER_H
