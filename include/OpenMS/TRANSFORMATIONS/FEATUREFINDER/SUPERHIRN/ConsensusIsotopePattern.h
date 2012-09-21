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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//

#ifndef CONSENS_ISOTOPE_PATTERN_H
#define CONSENS_ISOTOPE_PATTERN_H

#include <OpenMS/CONCEPT/Types.h>

#include <map>
#include <vector>

namespace OpenMS
{

	class OPENMS_DLLAPI ConsensusIsotopePattern
	{

			////////////////////////////////////////////////
			// declaration of the private members:

		private:

			// stores the consensus pattern:
			std::map<double, double> isotopesTrace_;
			std::vector<double> mzIsotopesStDev_;
			std::vector<double> intensIsotopesStDev_;

			// stores the detected patterns by retention time
			std::map<double, std::pair<std::vector<double>, std::vector<double> > > rawIsotopes_;

			////////////////////////////////////////////////
			// declaration of the public members:

		public:

			// class destructor
			~ConsensusIsotopePattern();

			// class constructor
			ConsensusIsotopePattern();
			// class copy constructor
			ConsensusIsotopePattern(const ConsensusIsotopePattern&);
			// class copy constructor
			ConsensusIsotopePattern(const ConsensusIsotopePattern*);

			//////////////////////////////////////////////////
			// overload operators:
			ConsensusIsotopePattern& operator=(const ConsensusIsotopePattern&);
			bool operator==(const ConsensusIsotopePattern&);
			ConsensusIsotopePattern& operator<=(const ConsensusIsotopePattern&);
			ConsensusIsotopePattern& operator>=(const ConsensusIsotopePattern&);
			ConsensusIsotopePattern& operator<(const ConsensusIsotopePattern&);
			ConsensusIsotopePattern& operator>(const ConsensusIsotopePattern&);

			// constructs the consensus pattern:
			void constructConsusPattern();
			// order an isotope trace in the correct cluster:
			void addIsotopeTrace(double, double);
			// condenses the pattern, make average peaks from the traces:
			void condensIsotopePattern(std::pair<std::vector<double>, std::vector<double> >*);

			///////////////////////////////
			// start here all the get / set
			// function to access the
			// variables of the class

			std::map<double, double>::iterator getConsensIsotopeIteratorStart();
			std::map<double, double>::iterator getConsensIsotopeIteratorEnd();

	};

	inline std::map<double, double>::iterator ConsensusIsotopePattern::getConsensIsotopeIteratorStart()
	{
		return isotopesTrace_.begin();
	}

	inline std::map<double, double>::iterator ConsensusIsotopePattern::getConsensIsotopeIteratorEnd()
	{
		return isotopesTrace_.end();
	}

} // ns

#endif

