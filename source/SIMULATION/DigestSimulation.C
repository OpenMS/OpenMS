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

#include<OpenMS/SIMULATION/DigestSimulation.h>

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

namespace OpenMS {

  DigestSimulation::DigestSimulation()
    : DefaultParamHandler("DigestSimulation")
  {
    defaults_.setValue("missed_cleavages",0,"maximum number of missed cleavages");
    defaults_.setValue("min_peptide_length",0,"minimum peptide length after digestion");

		defaultsToParam_();		
	}

  DigestSimulation::DigestSimulation(const DigestSimulation& source)
    : DefaultParamHandler(source)
  {}

  DigestSimulation& DigestSimulation::operator = (const DigestSimulation& source)
  {
    return *this;
  }
  
  DigestSimulation::~DigestSimulation()
  {}
  
  void DigestSimulation::digest(const SampleProteins & proteins, SamplePeptides & peptides)
  {
		// empty return vector
		peptides.clear();

		EnzymaticDigestion digestion;
    digestion.setEnzyme(EnzymaticDigestion::TRYPSIN);

    UInt missed_cleavages  = param_.getValue("missed_cleavages");
    UInt min_peptide_length = param_.getValue("min_peptide_length");

    digestion.setMissedCleavages( missed_cleavages );

		std::vector<AASequence> digestion_products;

    // Iterate through proteins and digest them
		for (SampleProteins::const_iterator protein = proteins.begin();
         protein != proteins.end();
         ++protein)
    {
      digestion.digest(AASequence(protein->first), digestion_products);

			for (std::vector<AASequence>::const_iterator dp_it = digestion_products.begin();
           dp_it != digestion_products.end();
           ++dp_it)
      {
        if (dp_it->size() < min_peptide_length) continue;

				peptides.push_back( std::make_pair<String, SimIntensityType>(dp_it->toUnmodifiedString(), protein->second) );
      }
		}

  }
}
