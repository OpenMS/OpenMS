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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/SampleTreatment.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	
	const std::string Sample::NamesOfSampleState[] = {"Unknown","solid","liquid","gas","solution","emulsion","suspension"};
	
	Sample::Sample() : 
		MetaInfoInterface(),
		state_(SAMPLENULL), 
		mass_(0.0), 
		volume_(0.0), 
		concentration_(0.0)
	{
	}
	
	Sample::Sample(const Sample& source):
		MetaInfoInterface(source),
	  name_(source.name_),
	  number_(source.number_),
	  comment_(source.comment_),
	  organism_(source.organism_),
	  state_(source.state_),
	  mass_(source.mass_),
	  volume_(source.volume_),
	  concentration_(source.concentration_),
	  subsamples_(source.subsamples_)
	{
		//delete old treatments
		for (list<SampleTreatment*>::iterator it = treatments_.begin(); it!=treatments_.end(); ++it)
		{
			delete(*it);
		}
		treatments_.clear();
		//copy treatments
		for (list<SampleTreatment*>::const_iterator it = source.treatments_.begin(); it!=source.treatments_.end(); ++it)
		{
			treatments_.push_back((*it)->clone());
		}
	}
	
	Sample::~Sample()
	{
		for (list<SampleTreatment*>::iterator it = treatments_.begin(); it!=treatments_.end(); ++it)
		{
			delete(*it);
		}
	}
	
	Sample& Sample::operator = (const Sample& source)
	{
	  if (&source == this) return *this;
	  
    name_ = source.name_;
    number_ = source.number_;
    comment_ = source.comment_;
    organism_ = source.organism_;
    state_ = source.state_;
    mass_ = source.mass_;
    volume_ = source.volume_;
    concentration_ = source.concentration_;
    subsamples_ = source.subsamples_;
    MetaInfoInterface::operator=(source);
		//delete old treatments
		for (list<SampleTreatment*>::iterator it = treatments_.begin(); it!=treatments_.end(); ++it)
		{
			delete(*it);
		}
		treatments_.clear();
		//copy treatments
		for (list<SampleTreatment*>::const_iterator it = source.treatments_.begin(); it!=source.treatments_.end(); ++it)
		{
			treatments_.push_back((*it)->clone());
		}
	  
	  return *this;
	}
	
	bool Sample::operator== (const Sample& rhs) const
	{
		if (
	   name_ != rhs.name_ ||
	   number_ != rhs.number_ ||
	   comment_ != rhs.comment_ ||
	   organism_ != rhs.organism_ ||
	   state_ != rhs.state_ ||
	   mass_ != rhs.mass_ ||
	   volume_ != rhs.volume_ ||
	   concentration_ != rhs.concentration_ ||
	   subsamples_ != rhs.subsamples_ ||
	   MetaInfoInterface::operator!=(rhs) ||
	   treatments_.size()!=rhs.treatments_.size()
	   )
		{
			return false;
	  }
	  
		//treatments
		list<SampleTreatment*>::const_iterator it2 = rhs.treatments_.begin();
		for (list<SampleTreatment*>::const_iterator it = treatments_.begin(); it!=treatments_.end(); ++it,++it2)
		{
			if (*it != *it2) return false;
		}
		return true;	
	}
	
	const String& Sample::getName() const 
	{
	  return name_; 
	}
	
	void Sample::setName(const String& name)
	{
	  name_ = name; 
	}
	
	const String& Sample::getOrganism() const 
	{
	  return organism_; 
	}
	
	void Sample::setOrganism(const String& organism)
	{
	  organism_ = organism; 
	}
	
	const String& Sample::getNumber() const 
	{
	  return number_; 
	}
	
	void Sample::setNumber(const String& number)
	{
	  number_ = number; 
	}
	
	const String& Sample::getComment() const 
	{
	  return comment_; 
	}
	
	void Sample::setComment(const String& comment)
	{
	  comment_ = comment; 
	}
	
	Sample::SampleState Sample::getState() const 
	{
	  return state_; 
	}
	
	void Sample::setState(Sample::SampleState state)
	{
	  state_ = state; 
	}
	
	DoubleReal Sample::getMass() const 
	{
	  return mass_; 
	}
	
	void Sample::setMass(DoubleReal mass)
	{
	  mass_ = mass; 
	}
	
	DoubleReal Sample::getVolume() const 
	{
	  return volume_; 
	}
	
	void Sample::setVolume(DoubleReal volume)
	{
	  volume_ = volume; 
	}
	
	DoubleReal Sample::getConcentration() const 
	{
	  return concentration_; 
	}
	
	void Sample::setConcentration(DoubleReal concentration)
	{
	  concentration_ = concentration; 
	}
	
	const vector<Sample>& Sample::getSubsamples() const
	{
	  return subsamples_;
	}
	
	vector<Sample>& Sample::getSubsamples()
	{
	  return subsamples_;
	}
	
	void Sample::setSubsamples(const vector<Sample>& subsamples)
	{
	  subsamples_ = subsamples;
	}
	
	
	void Sample::addTreatment(const SampleTreatment& treatment,Int before_position)
	{
		if (before_position > Int(treatments_.size()))
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,before_position,treatments_.size());
		}
		list<SampleTreatment*>::iterator it;
		if (before_position>=0)
		{
			it = treatments_.begin();
			for (Int i=0; i<before_position; ++i)
			{
				++it;
			}
		}
		else
		{
			it = treatments_.end();
		}
		SampleTreatment* tmp = treatment.clone();
		treatments_.insert(it,tmp);
	}
	
	const SampleTreatment& Sample::getTreatment(UInt position) const
	{
		if (position >= treatments_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,position,treatments_.size());
		}
		list<SampleTreatment*>::const_iterator it = treatments_.begin();
		for (Size i=0; i<position; ++i)
		{
			++it;
		}
		return **it;
	}
	
	SampleTreatment& Sample::getTreatment(UInt position)
	{
		if (position >= treatments_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,position,treatments_.size());
		}
		list<SampleTreatment*>::const_iterator it = treatments_.begin();
		for (Size i=0; i<position; ++i)
		{
			++it;
		}
		return **it;
	}
		
	void Sample::removeTreatment(UInt position)
	{
		if (position >= treatments_.size())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,position,treatments_.size());
		}
		list<SampleTreatment*>::iterator it = treatments_.begin();
		for (Size i=0; i<position; ++i)
		{
			++it;
		}
		delete(*it);
		treatments_.erase(it);
	}
	
	
	Int Sample::countTreatments() const
	{
		return (Int)treatments_.size();
	}

}


