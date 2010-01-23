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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_RAWTANDEMMSSIGNALSIMULATION_H
#define OPENMS_SIMULATION_RAWTANDEMMSSIGNALSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace OpenMS {

  /**
   @brief Simulates tandem MS signales for a given set of peptides

   Simulates tandem MS signales for a given set of peptides, with charge annotation,
   given detectabilities, predicted retention times and charge values.

   @htmlinclude OpenMS_RawTandemMSSignalSimulation.parameters

   @ingroup Simulation
  */
  class OPENMS_DLLAPI RawTandemMSSignalSimulation
    : public DefaultParamHandler
  {

		typedef ItraqConstants::ChannelInfo ChannelInfo;
		typedef ItraqConstants::ChannelMapType ChannelMapType;
		typedef ItraqConstants::IsotopeMatrices IsotopeMatrices;

  public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Constructor taking a random generator
    RawTandemMSSignalSimulation(const gsl_rng * random_generator);

    /// Copy constructor
    RawTandemMSSignalSimulation(const RawTandemMSSignalSimulation& source);

    /// Destructor
    virtual ~RawTandemMSSignalSimulation();
    //@}

    RawTandemMSSignalSimulation& operator = (const RawTandemMSSignalSimulation& source);

    /**

     */
    void generateRawTandemSignals(FeatureMapSim &, MSSimExperiment &);

  private:
		/// Default constructor (hidden)
    RawTandemMSSignalSimulation();

		/// init object (members, defaults etc)
		void init_();

    /// Synchronize members with param class
		void updateMembers_();

    /// Set default parameters
    void setDefaultParams_();

		/// convert meta information from feature into intensity values for iTRAQ
		Matrix<SimIntensityType> getItraqIntensity_(const Feature & f) const;

		// Members:

		/// set to either ItraqConstants::FOURPLEX or ItraqConstants::EIGHTPLEX
		Int itraq_type_;
		
		/// map the channel-name (eg 114) onto its description and the centroid mass
		/// the channel-name is also the id-string in the mapList section of the ConsensusMap
		ChannelMapType channel_map_;	

		/// Matrixes with isotope correction values (one for each plex-type)
		IsotopeMatrices isotope_corrections_;

  protected:
		/// Random number generator
		const gsl_rng* rnd_gen_;
  };

}

#endif
