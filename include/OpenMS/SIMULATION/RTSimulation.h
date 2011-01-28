// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{
  /**
   @brief Simulates/Predicts retention times for peptides or peptide separation.
   
   The retention times for the different peptides are determined based on
   a SVM model or are all set to -1 in case of simulations without a HPLC
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
    RTSimulation(const SimRandomNumberGenerator& random_generator);
    
    /// Copy constructor
    RTSimulation(const RTSimulation& source);

    /// Destructor
    virtual ~RTSimulation();
    //@}

    /// Assignment operator
    RTSimulation& operator = (const RTSimulation& source);
    
    /** 
     @brief Predict retention times for given peptide features based for HPLC or CE
     
     @param features Feature map for which the retention times will be predicted
     */
    void predictRT(FeatureMapSim & features);
 
    /**
     @brief Set retention times randomly for given contaminants
     */
    void predictContaminantsRT(FeatureMapSim &);
    
    /**
     @brief Returns true if a RT column was simulated
     */
    bool isRTColumnOn() const;

		/// Wrapper for the SVM RT Prediction (HPLC) using AASequences
		void wrapSVM(std::vector<AASequence>& peptide_sequences,std::vector<DoubleReal>& predicted_retention_times);

    SimCoordinateType getGradientTime() const;

    /// Size experiment and assign retention time grid
    void createExperiment(MSSimExperiment & experiment);

  private:
    /// Default constructor
    RTSimulation();
    
    /// Set default parameters
    void setDefaultParams_();

    /// Simply set all retention times to -1
    void noRTColumn_(FeatureMapSim &);
    
		/// smoothes the simulated distortion for the elution profiles with a moving average filter of size 3
		void smoothRTDistortion_(MSSimExperiment & experiment);

		/// Wrapper for the Migration time calculation (CE)
    /// @param features will get modified with metavalue "RT_CE_width_factor", describing widening of MT shape
		void calculateMT_(FeatureMapSim & features,std::vector<DoubleReal>& predicted_retention_times);

		void getChargeContribution_(Map< String, double> & q_cterm, 
															  Map< String, double> & q_nterm,
															  Map< String, double> & q_aa_basic,
															  Map< String, double> & q_aa_acidic);
    
		// MEMBERS:

    // Name of the svm model file
		OpenMS::String rt_model_file_;
    
    /// Total gradient time
    SimCoordinateType total_gradient_time_;

    /// gradient ranges

    /// Minimal observed gradient time
    SimCoordinateType gradient_min_;
    /// Maximal observed gradient time
    SimCoordinateType gradient_max_;

    /// bin size in rt dimension
    SimCoordinateType rt_sampling_rate_;

    /// EGH tau value
    DoubleReal egh_tau_location_;
    /// EGH tau scale parameter of the lorentzian variation
    DoubleReal egh_tau_scale_;

    /// EGH sigma value
    DoubleReal egh_variance_location_;
    /// EGH sigma scale parameter of the lorentzian variation
    DoubleReal egh_variance_scale_;

  protected:  
		/// Random number generator
    SimRandomNumberGenerator const * rnd_gen_;
    
    /// Synchronize members with param class
		void updateMembers_();
    
  };

}

#endif
