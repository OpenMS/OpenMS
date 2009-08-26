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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_DIGESTSIMULATION_H
#define OPENMS_SIMULATION_DIGESTSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{

  /**
		@brief Simulates protein digestion

		Supports all enzymes supported by EnzymaticDigestion.h
		and additionally incorporates abundance values, which
		are distributed evenly among digestion products of each
		protein.

		@htmlinclude OpenMS_DigestSimulation.parameters

		@ingroup Simulation
  */
  class OPENMS_DLLAPI DigestSimulation
    : public DefaultParamHandler
  {

  public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Default constructor
    DigestSimulation();

    /// Copy constructor
    DigestSimulation(const DigestSimulation& source);

    /// Destructor
    virtual ~DigestSimulation();
    //@}

		/// Assignment operator
		DigestSimulation& operator = (const DigestSimulation& source);


		/**
			@brief Digest a set of proteins into peptides

			Digest proteins to peptides, with protein abundance distributes equally among
			created sibling peptides (this also applies for peptides with missed cleavages).
			Should a peptide be non-unique the abundances of its instances from proteins are summed up.
			
			@param feature_map Input FeatureMap containing the proteins that should be digested as ProteinIdentification
		**/
    void digest(FeatureMapSim & feature_map);    
    
  private:
    /// set defaults
    void setDefaultParams_();

  };

}

#endif
