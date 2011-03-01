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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
	ConsensusFeature::ConsensusFeature()
		: BaseFeature(), HandleSetType()
	{
	}

	ConsensusFeature::ConsensusFeature(const ConsensusFeature& rhs)
		: BaseFeature(rhs), HandleSetType(rhs)
	{
	}

	ConsensusFeature::ConsensusFeature(const BaseFeature& feature)
		: BaseFeature(feature), HandleSetType()
	{
	}

	ConsensusFeature::ConsensusFeature(UInt64 map_index, const Peak2D& element, UInt64 element_index)
		: BaseFeature(element), HandleSetType()
	{
		insert(map_index, element, element_index);
	}

	ConsensusFeature::ConsensusFeature(UInt64 map_index, const BaseFeature& element)
		: BaseFeature(element), HandleSetType()
	{
		insert(FeatureHandle(map_index,element));
	}

	ConsensusFeature& ConsensusFeature::operator=(const ConsensusFeature& rhs)
	{
		if (&rhs==this) return *this;

		HandleSetType::operator=(rhs);
		BaseFeature::operator=(rhs);

		return *this;
	}

	ConsensusFeature::~ConsensusFeature()
	{
	}
		
	void ConsensusFeature::insert(const FeatureHandle& handle)
	{
		if (!(HandleSetType::insert(handle).second))
		{
			String key = String("map") + handle.getMapIndex() + "/feature" + handle.getUniqueId();
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set already contained an element with this key.",key) ;
		}
	}

	void ConsensusFeature::insert(const HandleSetType& handle_set)
	{
		for (ConsensusFeature::HandleSetType::const_iterator it = handle_set.begin(); it != handle_set.end(); ++it)
		{
			insert(*it);
		}
	}

	void ConsensusFeature::insert(UInt64 map_index, const Peak2D& element, UInt64 element_index)
	{
		insert(FeatureHandle(map_index,element,element_index));
	}

	void ConsensusFeature::insert(UInt64 map_index, const BaseFeature& element)
	{
		insert(FeatureHandle(map_index,element));
		peptides_.insert(peptides_.end(), element.getPeptideIdentifications().begin(), element.getPeptideIdentifications().end());
	}

	const ConsensusFeature::HandleSetType& ConsensusFeature::getFeatures() const
	{
		return *this;
	}

	DRange<2> ConsensusFeature::getPositionRange() const
	{
		DPosition<2> min = DPosition<2>::maxPositive();
		DPosition<2> max = DPosition<2>::minPositive();
		for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
		{
			if (it->getRT()<min[0]) min[0]=it->getRT();
			if (it->getRT()>max[0]) max[0]=it->getRT();
			if (it->getMZ()<min[1]) min[1]=it->getMZ();
			if (it->getMZ()>max[1]) max[1]=it->getMZ();
		}
		return DRange<2>(min,max);
	}
	
	DRange<1> ConsensusFeature::getIntensityRange() const
	{
		DPosition<1> min = DPosition<1>::maxPositive();
		DPosition<1> max = DPosition<1>::minPositive();
		for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
		{
			if (it->getIntensity()<min[0]) min[0]=it->getIntensity();
			if (it->getIntensity()>max[0]) max[0]=it->getIntensity();
		}
		return DRange<1>(min,max);
	}

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

	void ConsensusFeature::computeMonoisotopicConsensus()
	{
		// for computing average rt position, minimal m/z position and intensity
		DoubleReal rt=0.0;
		DoubleReal mz=std::numeric_limits<DoubleReal>::max();
		DoubleReal intensity=0.0;

		// The most frequent charge state wins.  Tie breaking prefers smaller charge.
		std::map<Int,UInt> charge_occ;
		Int charge_most_frequent = 0;
		UInt charge_most_frequent_occ = 0;

		for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
		{
			rt += it->getRT();
			if (it->getMZ() < mz)
				mz=it->getMZ();
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

		// compute the position and intensity
		setRT(rt / size());
		setMZ(mz);
		setIntensity(intensity / size());
		setCharge(charge_most_frequent);
		return;
	}
  
  void ConsensusFeature::computeDechargeConsensus(const FeatureMap<>& fm, bool intensity_weighted_averaging)
  {
  	// for computing average position and intensity
  	DoubleReal rt=0.0;
  	DoubleReal m=0.0;
  	DoubleReal intensity=0.0;
		
		DoubleReal proton_mass = Constants::PROTON_MASS_U;

		// intensity sum (for weighting)
    for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
    {
			intensity += it->getIntensity();
    }

		// unweighted averaging by default
		// TODO: add outlier removal
		// TODO: split cluster for each channel (in FD.C)
		DoubleReal weighting_factor = 1.0/size();
		
		// RT and Mass
    for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
    {
			Int q = it->getCharge();
			if (q==0) LOG_WARN << "ConsensusFeature::computeDechargeConsensus() WARNING: Feature's charge is 0! This will lead to M=0!\n";
    	DoubleReal adduct_mass;
    	Size index=fm.uniqueIdToIndex(it->getUniqueId());
    	if (index > fm.size()) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__, index, fm.size());
    	if (fm[index].metaValueExists("dc_charge_adduct_mass"))
    	{
    		adduct_mass = (DoubleReal) fm[index].getMetaValue("dc_charge_adduct_mass");
    	}
    	else
    	{
    		adduct_mass = q * proton_mass;
    	}

			if (intensity_weighted_averaging) weighting_factor = it->getIntensity() / intensity;
    	rt += it->getRT() * weighting_factor;
    	m += (it->getMZ() * q - adduct_mass) * weighting_factor;
    }

    // compute the average position and intensity
    setRT(rt);
    setMZ(m);
    setIntensity(intensity);
    setCharge(0);
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
         << "   Feature id: " << it->getUniqueId() << std::endl
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
