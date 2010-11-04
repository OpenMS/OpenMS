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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_MSSIM_H
#define OPENMS_SIMULATION_MSSIM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{


  /**
   @brief Central class for simulation of mass spectrometry experiments
   
   This implementation is an extended and rewritten version of the concepts and ideas presented in:<br>
   <p>
    Ole Schulz-Trieglaff, Nico Pfeifer, Clemens Gropl, Oliver Kohlbacher, and Knut Reinert.<br>
    LC-MSsim - A simulation software for liquid chromatography mass spectrometry data.<br>
    <em>BMC Bioinformatics</em> 9:423, 2008.
   </p>
   
   @htmlinclude OpenMS_MSSim.parameters
   
   @see DetectabilitySimulation
   @see DigestSimulation
   @see IonizationSimulation
   @see RawMSSignalSimulation
   @see RTSimulation
   
   @ingroup Simulation
  */
  class OPENMS_DLLAPI MSSim
    : public DefaultParamHandler,
      public ProgressLogger
  {

  public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Default constructor
    MSSim();

    /// Destructor
    virtual ~MSSim();
    //@}

    /**
     @brief General purpose function to simulate a mass spectrometry run
     
     @param rnd_gen GSL random number generator which will be passed to the different classes
     @param peptides List of peptides and abundances that will be simulated
     @param labeling_type
     */   
    void simulate(const SimRandomNumberGenerator & rnd_gen, SampleChannels& peptides, const String& labeling_tpye);
	
    /// Access the simulated experiment
    MSSimExperiment const & getExperiment() const;
    
    /// Access the simulated features
    FeatureMapSim const & getSimulatedFeatures() const;

		/// Access the charge consensus map of simulated features
    ConsensusMap const & getChargeConsensus() const;

		/// Access the contaminants feature mapmap of simulated features
    FeatureMapSim const & getContaminants() const;

    /// Access the labeling consensus map of simulated features
    ConsensusMap const & getLabelingConsensus() const;

    /// Returns the default parameters for simulation including the labeling technique with name @p labeling_name
    Param getParameters(const String& labeling_name) const;
  protected:
		/// handle global params
		void syncParams_(Param& p, bool to_outer);
		
		/// Convert a list of peptides with given abundance values into a FeatureMap
		void createFeatureMap_(const SampleProteins& peptides, FeatureMapSim& features);

  private:
    /// Synchronize members with param class
		void updateMembers_();        
    
    MSSimExperiment experiment_;
    
    FeatureMapSimVector feature_maps_;

		ConsensusMap consensus_map_;

    FeatureMapSim contaminants_map_;

    /// Labeling functionality
    BaseLabeler * labeler_;
  };

}

#endif
