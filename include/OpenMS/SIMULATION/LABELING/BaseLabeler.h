// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LABELING_BASELABELER_H
#define OPENMS_SIMULATION_LABELING_BASELABELER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/SIMULATION/SimTypes.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace OpenMS
{

  /**
  @brief Abstract base class for all kinds of labeling techniques
  */
  class OPENMS_DLLAPI BaseLabeler
    : public DefaultParamHandler
  {
  public:

    /// constructor
    BaseLabeler();

    /// destructor
    virtual ~BaseLabeler()
    {
    }

    /// register all derived classes here (implemented in file BaseLabeler_impl.h)
    static void registerChildren();    
        
    /**
      @brief Returns the default parameters. Reimplement

      Reimplement if you derive a class and have to incoorporate sub-algorithm default parameters.
    */
    virtual Param getDefaultParameters() const
    {
      return this->defaults_;
    }

    virtual void setRnd(const SimRandomNumberGenerator& rng)
    {
      rng_ = &rng;
    }

    /**
      @brief Checks the (simulation) params passed if they are consistent with
      the labeling technique.

      @param param Param object containing the simulation parameters
      @throws Exception::InvalidParameter if the given parameters are not consistent with the labeling technique
      */
    virtual void preCheck(Param &) const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /** 
    @name Labeling Hooks
    */
    //@{
    
    /// Hook to prepare the simulation process
    virtual void setUpHook(FeatureMapSimVector & /* features */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /// Labeling between digestion and rt simulation
    virtual void postDigestHook(FeatureMapSimVector & /* features_to_simulate */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /// Labeling after rt simulation
    virtual void postRTHook(FeatureMapSimVector & /* features_to_simulate */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /// Labeling after detectability simulation
    virtual void postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /// Labeling after ionization
    virtual void postIonizationHook(FeatureMapSimVector & /* features_to_simulate */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /// Labeling after raw signal generation
    virtual void postRawMSHook(FeatureMapSimVector & /* features_to_simulate */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /// Labeling after Tandem MS (e.g. iTRAQ)
    virtual void postRawTandemMSHook(FeatureMapSimVector & /* features_to_simulate */, MSSimExperiment & /* simulated map */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    //@}

    const ConsensusMap& getConsensus() const;

    /**
      @brief to ensure standardized metavalue names across labelers for channel intensity

      Use this function to get the name of the metavalue which holds intensity for channel @p channel_index

    */
    String getChannelIntensityName(const Size channel_index) const;


  protected:
    /**
      @brief Creates an empty FeatureMap with the merged ProteinIdentifications from
      all FeatureMaps contained in @p maps

      @param maps       Vector of FeatureMaps containing the features that will be merged
      @return           A FeatureMap containing all ProteinIdentifications of the input maps
      */
    FeatureMapSim mergeProteinIdentificationsMaps_(const FeatureMapSimVector & maps);

    /**
      @brief join all protein references of two features

      When merging peptides from different channels, the protein accessions should remain intact.
      Usually joining features is based on peptide sequence, so all protein hits should be valid.

    */
    void mergeProteinAccessions_(Feature& target, const Feature& source) const;

    ConsensusMap consensus_;

    SimRandomNumberGenerator const * rng_;

    /**
      @brief Based on the stored consensus recompute the associations for the passed features, assuming
             that the features where derived from the features stored in the consensus.

      @param simulated_features FeatureMap containing features derived from the ones, stored in the
                                consensus
    */
    void recomputeConsensus_(const FeatureMapSim & simulated_features);
  };
} // namespace OpenMS

#endif //#ifndef OPENMS_SIMULATION_LABELING_BASELABELER_H
