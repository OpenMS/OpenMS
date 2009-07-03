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

#ifndef OPENMS_SIMULATION_PTMSIMULATION_H
#define OPENMS_SIMULATION_PTMSIMULATION_H

#include <vector>
#include <utility>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>


namespace OpenMS {

  /**
		@brief Adds PTM to peptides

		This class can add PTMs to a set of peptides.
		Up to one PTM version of a peptide is generated randomly
		according to allowed modifications.

		@htmlinclude OpenMS_PTMSimulation.parameters

		@ingroup Simulation
  */
  class OPENMS_DLLAPI PTMSimulation
    : public DefaultParamHandler
  {

  public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Constructor
    PTMSimulation(const gsl_rng* rnd_gen);

    /// Copy constructor
    PTMSimulation(const PTMSimulation& source);

    /// Destructor
    virtual ~PTMSimulation();
    //@}

		/// Assignment operator
    PTMSimulation& operator = (const PTMSimulation& source); 

    /**
      @brief Adds possible post translational modifications to the peptides

			Up to 'modification_bound' PTMS are added (position independent) according to 
			allowed modifications given in 'potential_modifications'.

			Up to one new feature per input feature is added, i.e. 0 to 1 modified versions of a peptide are possible.

			@param map FeatureMap which is extended by diced PTM's

     */
    void predictPTMs(FeatureMapSim & map);
   
	protected:
		/// see
		void updateMembers_();

		/// valid PTM's
		PTMTable ptms_;

		/// random generator (GSL)
		const gsl_rng* rnd_gen_;

	private:
		/// Default (hidden)
    PTMSimulation();
    
    /// see
    void setDefaultParams_();

  };
}

#endif
