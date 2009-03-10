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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConsensusFeature.h>

namespace OpenMS
{

  void ConsensusFeature::computeConsensus()
  {
  	// for computing average position and intensity
  	DoubleReal rt=0.0;
  	DoubleReal mz=0.0;
  	DoubleReal intensity=0.0;

  	// The most frequent charge state wins.  Tie breaking prefers smaller charge.
  	std::map<Int,UInt> charge_occ;
  	Int charge_most_frequent = 0;
  	UInt charge_most_frequent_occ = 0;

    for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
    {
    	rt += it->getRT();
    	mz += it->getMZ();
    	intensity += it->getIntensity();
    	const Int it_charge = it->getCharge();
    	const UInt it_charge_occ = ++charge_occ[it_charge];
    	if ( it_charge_occ > charge_most_frequent_occ )
    	{
    	  charge_most_frequent_occ = it_charge_occ;
        charge_most_frequent = it_charge;
    	}
    	else
    	{
        if ( it_charge_occ >= charge_most_frequent_occ && abs(it_charge) < abs(charge_most_frequent) )
        {
          charge_most_frequent = it_charge;
        }
    	}
    }

    // compute the average position and intensity
    setRT(rt / size());
    setMZ(mz / size());
    setIntensity(intensity / size());
    setCharge(charge_most_frequent);
    return;
  }

  std::ostream& operator << (std::ostream& os, const ConsensusFeature& cons)
  {
    os << "---------- CONSENSUS ELEMENT BEGIN -----------------\n";
    os << "Position: " << cons.getPosition()<< std::endl;
    os << "Intensity " << precisionWrapper(cons.getIntensity()) << std::endl;
    os << "Quality " << precisionWrapper(cons.getQuality()) << std::endl;
    os << "Grouped features: " << std::endl;

    for (ConsensusFeature::HandleSetType::const_iterator it = cons.begin(); it != cons.end(); ++it)
    {
      os << " - Map index: " << it->getMapIndex() << std::endl
         << "   Feature index: " << it->getElementIndex() << std::endl
      	 << "   RT: " << precisionWrapper(it->getRT()) << std::endl
      	 << "   m/z: " << precisionWrapper(it->getMZ())  << std::endl
      	 << "   Intensity: " << precisionWrapper(it->getIntensity()) << std::endl;
    }

		os << "Meta information: " << std::endl;
		std::vector< String > keys;
		cons.getKeys(keys);
    for (std::vector< String >::const_iterator it = keys.begin(); it != keys.end(); ++it)
    {
			os << "   " << (*it) << ": " << cons.getMetaValue (*it) << std::endl;
		}
    os << "---------- CONSENSUS ELEMENT END ----------------- " << std::endl;

    return os;
  }
}
