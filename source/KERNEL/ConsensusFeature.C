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
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
	ConsensusFeature::ConsensusFeature()
		: RichPeak2D(),
			HandleSetType(),
			quality_(0.0),
			charge_(0),
			peptide_identifications_()
	{
	}

	ConsensusFeature::ConsensusFeature(const ConsensusFeature& rhs)
		: RichPeak2D(rhs),
			HandleSetType(rhs),
			quality_(rhs.quality_),
			charge_(rhs.charge_),
			peptide_identifications_(rhs.peptide_identifications_)
	{
	}

	ConsensusFeature::ConsensusFeature(const RichPeak2D& point)
		: RichPeak2D(point),
			HandleSetType(),
			quality_(0.0),
			charge_(0),
			peptide_identifications_()
	{
	}

	ConsensusFeature::ConsensusFeature(const Peak2D& point)
		: RichPeak2D(point),
			HandleSetType(),
			quality_(0.0),
			charge_(0),
			peptide_identifications_()
	{
	}

	ConsensusFeature::ConsensusFeature(const Feature& feature)
		: RichPeak2D(feature),
			HandleSetType(),
			quality_(0.0),
			charge_(0),
			peptide_identifications_(feature.getPeptideIdentifications())
	{
	}

	ConsensusFeature::ConsensusFeature(Size map_index,	Size element_index, const Peak2D& element)
		: RichPeak2D(element),
			HandleSetType(),
			quality_(0.0),
			charge_(0),
			peptide_identifications_()
	{
		insert(map_index,element_index,element);
	}

	ConsensusFeature::ConsensusFeature(Size map_index,	Size element_index, const Feature& element)
		: RichPeak2D(element),
			HandleSetType(),
			quality_(element.getOverallQuality()),
			charge_(element.getCharge()),
			peptide_identifications_()
	{
		insert(map_index,element_index,element);
	}

	ConsensusFeature::ConsensusFeature(Size map_index,	Size element_index, const ConsensusFeature& element)
		: RichPeak2D(element),
			HandleSetType(),
			quality_(element.getQuality()),
			charge_(element.getCharge()),
			peptide_identifications_()
	{
		insert(map_index,element_index,element);
	}

	ConsensusFeature& ConsensusFeature::operator=(const ConsensusFeature& rhs)
	{
		if (&rhs==this) return *this;

		HandleSetType::operator=(rhs);
		RichPeak2D::operator=(rhs);
		quality_ = rhs.quality_;
		charge_ = rhs.charge_;
		peptide_identifications_ =  rhs.peptide_identifications_;

		return *this;
	}

	ConsensusFeature::~ConsensusFeature()
	{
	}
		

	/**
	@brief Adds an feature handle into the consensus feature

	@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
	*/
	void ConsensusFeature::insert(const FeatureHandle& handle)
	{
		if (!(HandleSetType::insert(handle).second))
		{
			String key = String("map") + handle.getMapIndex() + "/feature" + handle.getElementIndex();
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set already contained an element with this key.",key) ;
		}
	}

	/// Adds all feature handles in @p handle_set to this consensus feature.
	void ConsensusFeature::insert(const HandleSetType& handle_set)
	{
		for (ConsensusFeature::HandleSetType::const_iterator it = handle_set.begin(); it != handle_set.end(); ++it)
		{
			insert(*it);
		}
	}

	/**
		@brief Creates a FeatureHandle and adds it

		@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
	*/
	void ConsensusFeature::insert(Size map_index, Size element_index, const Peak2D& element)
	{
		insert(FeatureHandle(map_index,element_index,element));
	}

	/**
		@brief Creates a FeatureHandle and adds it

		@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
	*/
	void ConsensusFeature::insert(Size map_index, Size element_index, const Feature& element)
	{
		insert(FeatureHandle(map_index,element_index,element));
		peptide_identifications_.insert(peptide_identifications_.end(), element.getPeptideIdentifications().begin(), element.getPeptideIdentifications().end());
	}

	/**
		@brief Creates a FeatureHandle and adds it

		@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
	*/
	void ConsensusFeature::insert(Size map_index, Size element_index, const ConsensusFeature& element)
	{
		insert(FeatureHandle(map_index,element_index,element));
		peptide_identifications_.insert(peptide_identifications_.end(), element.getPeptideIdentifications().begin(), element.getPeptideIdentifications().end());
	}

	/// Non-mutable access to the contained feature handles
	const ConsensusFeature::HandleSetType& ConsensusFeature::getFeatures() const
	{
		return *this;
	}

	/// Returns the quality
	ConsensusFeature::QualityType ConsensusFeature::getQuality() const
	{
		return quality_;
	}
	/// Sets the quality
	void ConsensusFeature::setQuality(ConsensusFeature::QualityType quality)
	{
		quality_ = quality;
	}
	/// Sets the charge
	void ConsensusFeature::setCharge(Int charge)
	{
		charge_ = charge;
	}
	/// Returns the charge
	Int ConsensusFeature::getCharge() const
	{
		return charge_;
	}
	/// Returns the position range of the contained elements
	DRange<2> ConsensusFeature::getPositionRange() const
	{
		DPosition<2> min = DPosition<2>::max();
		DPosition<2> max = DPosition<2>::min();
		for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
		{
			if (it->getRT()<min[0]) min[0]=it->getRT();
			if (it->getRT()>max[0]) max[0]=it->getRT();
			if (it->getMZ()<min[1]) min[1]=it->getMZ();
			if (it->getMZ()>max[1]) max[1]=it->getMZ();
		}
		return DRange<2>(min,max);
	}
	/// Returns the intensity range of the contained elements
	DRange<1> ConsensusFeature::getIntensityRange() const
	{
		DPosition<1> min = DPosition<1>::max();
		DPosition<1> max = DPosition<1>::min();
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
  
  void ConsensusFeature::computeDechargeConsensus(const FeatureMap<>& fm)
  {
  	// for computing average position and intensity
  	DoubleReal rt=0.0;
  	DoubleReal m=0.0;
  	DoubleReal intensity=0.0;
		
		DoubleReal proton_mass = ElementDB::getInstance()->getElement("H")->getMonoWeight();

    for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
    {
    	rt += it->getRT();
			Int q = it->getCharge();
			if (q==0) std::cerr << "ConsensusFeature::computeDechargeConsensus() WARNING: Feature's charge is 0! This will lead to M=0!\n";
    	DoubleReal adduct_mass;
    	if (it->getElementIndex() > fm.size()) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__, it->getElementIndex(), fm.size());
    	if (fm[it->getElementIndex()].metaValueExists("dc_charge_adduct_mass"))
    	{
    		adduct_mass = (DoubleReal) fm[it->getElementIndex()].getMetaValue("dc_charge_adduct_mass");
    	}
    	else
    	{
    		adduct_mass = q * proton_mass;
    	}
    	m += it->getMZ() * q - adduct_mass;
    	intensity += it->getIntensity();
    }

    // compute the average position and intensity
    setRT(rt / size());
    setMZ(m / size());
    setIntensity(intensity);
    setCharge(0);
    return;
  }  

	/// returns a const reference to the PeptideIdentification vector
	const std::vector<PeptideIdentification>& ConsensusFeature::getPeptideIdentifications() const
	{
		return peptide_identifications_;
	}

	/// returns a mutable reference to the PeptideIdentification vector
	std::vector<PeptideIdentification>& ConsensusFeature::getPeptideIdentifications()
	{
		return peptide_identifications_;
	}

	/// sets the PeptideIdentification vector
	void ConsensusFeature::setPeptideIdentifications( const std::vector<PeptideIdentification>& peptide_identifications )
	{
		peptide_identifications_ = peptide_identifications;
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
