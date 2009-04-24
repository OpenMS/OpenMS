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

#ifndef OPENMS_SIMULATION_RAWSIGNALSIMULATION_H
#define OPENMS_SIMULATION_RAWSIGNALSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace OpenMS {

  /**
   @brief
   @ingroup Simulation
  */
  class OPENMS_DLLAPI RawSignalSimulation
    : public DefaultParamHandler
  {

  public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Constructor taking a random generator
    RawSignalSimulation(const gsl_rng * random_generator);

    /// Copy constructor
    RawSignalSimulation(const RawSignalSimulation& source);

    /// Destructor
    virtual ~RawSignalSimulation();
    //@}

    RawSignalSimulation& operator = (const RawSignalSimulation& source);

    // TODO: howto add contaminations
    void generateRawSignals(FeatureMapSim &, MSSimExperiment &);

    SimCoordinateType getRTSamplingRate() const;
  private:
    /// Default constructor
    RawSignalSimulation();

    /// Synchronize members with param class
		void updateMembers_();

    ///
    void setDefaultParams_();

    ///
    void addMSSignal(Feature &, MSSimExperiment &);

    void samplePeptideModel_(const ProductModel<2> & pm,
                             const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
                             SimCoordinateType rt_start, SimCoordinateType rt_end,
                             MSSimExperiment &, Feature & activeFeature);

    void chooseElutionProfile_(ProductModel<2>& pm, const SimCoordinateType rt,const double scale);

    
    void addShotNoise_(MSSimExperiment & experiment);
    
    void addBaseLine_(MSSimExperiment & experiment);
    
    void compressSignals_(MSSimExperiment & experiment);
     
    Size compressSignalsRun_(MSSimExperiment & experiment);
    // TODO: the following parameters are imported -> revise
    // TODO: we need to incorporate those parameters into constructors etc.
    
		/// bin sizes
		SimCoordinateType mz_sampling_rate_;
    SimCoordinateType rt_sampling_rate_;

		/// Mean of peak m/z error
		SimCoordinateType mz_error_mean_;
		/// Standard deviation of peak m/z error
		SimCoordinateType mz_error_stddev_;

		/// Mean of peak intensity error
		SimIntensityType intensity_error_mean_;
		/// Standard deviation of peak intensity error
		SimIntensityType intensity_error_stddev_;

		/// Maximum m/z detected by mass analyser
		SimCoordinateType maximal_mz_measurement_limit_;
		/// Minimum m/z detected by mass analyser
		SimCoordinateType minimal_mz_measurement_limit_;

    /// Full width at half maximum of simulated peaks
		SimCoordinateType peak_std_;

    /// Mean intensity scaling
    SimCoordinateType mean_scaling_;

    /// Number of peptide ions
    Size ion_count_;

		/// Remembers which scans were changed after the last call to removeDuplicatePoints_()
		std::vector<bool> changed_scans_;
    
    /// LC conditions (noise parameter for EMG)
		DoubleReal distortion_;
    DoubleReal symmetry_up_;
    DoubleReal symmetry_down_;
    
  protected:
		/// Random number generator
		const gsl_rng* rnd_gen_;
  };

}

#endif
