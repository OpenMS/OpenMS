// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Stephan Aiche Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_MSSIM_H
#define OPENMS_SIMULATION_MSSIM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS {

  /**
   @brief Central class for simulation of mass spectrometry experiments
   
   This implementation is based on the concepts and ideas presented in:
   LC-MSsim paper
   
   @ingroup Simulation
  */
  class OPENMS_DLLAPI MSSim
    : public DefaultParamHandler
  {

  public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Default constructor
    MSSim();

    /// Copy constructor
    MSSim(const MSSim& source);

    /// Destructor
    virtual ~MSSim();
    //@}

    MSSim& operator = (const MSSim& source);

    /**
     @brief General purpose function to simulate a mass spectrometry run
     */   
    void simulate(const gsl_rng* rnd_gen, const SamplePeptides& peptides);
	
    /**
     @brief Access the simulated experiment
     */
    MSSimExperiment const & getExperiment() const;
    
    /**
     @brief Access the simulated features
     */
    FeatureMapSim const & getSimulatedFeatures() const;
	protected:
		
		/// convert list of peptides with abundance into a FeatureMap
		void createFeatureMap_(const SamplePeptides& peptides);

    /// generates a MSSimExperiment of correct size
    void createExperiment_(const DoubleReal gradient_time, const DoubleReal rt_sampling_rate);

  private:
    /// set defaults
    void setDefaultParams_();
    
    /// Synchronize members with param class
		void updateMembers_();        
    
    MSSimExperiment experiment_;
    
    FeatureMapSim features_;
  };

}

#endif
