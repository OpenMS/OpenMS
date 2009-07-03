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
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_RTSIMULATION_H
#define OPENMS_SIMULATION_RTSIMULATION_H

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{
  /**
   @brief Simulates/Predicts predict retention times 
   for peptides or peptide separation.
   
   The retention times for the different peptides are determined based on
   a SVM model or are all set to -1 in case of simulations without a LC
   column.
   
   @htmlinclude OpenMS_RTSimulation.parameters
   
   @ingroup Simulation
  */
  class OPENMS_DLLAPI RTSimulation
    : public DefaultParamHandler
  {

  public:
    /** @name Constructors and Destructors
      */
    //@{

    /// Constructor taking a random generator
    RTSimulation(const gsl_rng * random_generator);
    
    /// Copy constructor
    RTSimulation(const RTSimulation& source);

    /// Destructor
    virtual ~RTSimulation();
    //@}

    // Assignment operator
    RTSimulation& operator = (const RTSimulation& source);
    
    /** 
     @brief Predict retention times for given peptide features based on a SVM Model
     
     @param features Feature map for which the retention times will be predicted
     @param features Experiment map which will be build from scratch
     */
    void predictRT(FeatureMapSim & features, MSSimExperiment & experiment);
 
    /**
     @brief Set retention times randomly for given contaminants
     */
    void predictContaminantsRT(FeatureMapSim &);
    
    /**
     @brief Returns true if a RT column was simulated
     */
    bool isRTColumnOn() const;
    
    SimCoordinateType getGradientTime() const;
  private:
    /// Default constructor
    RTSimulation();
    
    /// Set default parameters
    void setDefaultParams_();

    /// Simply set all retention times to -1
    void noRTColumn_(FeatureMapSim &);
    
    /// Predict all retention times based on a svm model
    void predictFeatureRT_(FeatureMapSim &);
  
		/// Size experiment and assign retention time grid
		void createExperiment_(MSSimExperiment & experiment);
    
    // MEMBERS:
    
    // Name of the svm model file
		OpenMS::String rt_model_file_;
    
    /// Total gradient time
    SimCoordinateType gradient_time_;

    /// bin size in rt dimension
    SimCoordinateType rt_sampling_rate_;

    /// LC conditions (noise parameter for EMG)
		DoubleReal distortion_;
    DoubleReal symmetry_up_;
    DoubleReal symmetry_down_;
    
    /// Front part of the LC gradient that will not be directly assigned to guarantee a full elution profile
    static const DoubleReal gradient_front_offset_;
    /// Total part (front + back) of the LC gradient that will not be directly assigned
    static const DoubleReal gradient_total_offset_;
    
  protected:  
		/// Random number generator
		const gsl_rng* rnd_gen_;    
    
    /// Synchronize members with param class
		void updateMembers_();
    
  };

}

#endif
