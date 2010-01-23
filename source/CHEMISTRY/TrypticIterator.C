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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------
 
#include <OpenMS/CHEMISTRY/TrypticIterator.h>
#include <OpenMS/FORMAT/FastaIterator.h>
#include <fstream>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>


namespace OpenMS
{


	typedef std::pair<String,String> FASTAEntry;
	
	///Constructor to intialize algorithm
	TrypticIterator::TrypticIterator ()
		: PepIterator()
	{
		b_ = 0;
		e_ = 0;
		is_at_end_ = false;
		
	}

	///destructor
	TrypticIterator::~TrypticIterator ()
	{
	
	}
	
	TrypticIterator::TrypticIterator (const TrypticIterator & source) :
		PepIterator (source),
		f_file_(source.f_file_),
		actual_pep_(source.actual_pep_),
		is_at_end_(source.is_at_end_),
		f_entry_(source.f_entry_),
		b_(source.b_),
		e_(source.e_)
	{
		f_iterator_=source.f_iterator_;
	}

	FASTAEntry TrypticIterator::operator*()
	{
		if (actual_pep_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		return FASTAEntry(f_entry_.first,actual_pep_);
	}
	
	PepIterator & TrypticIterator::operator++()
	{
		if (actual_pep_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		actual_pep_ = next_();
		
		if (f_iterator_->isAtEnd() && !hasNext_()) 
		{
			is_at_end_ = true;
		}
		return *this;
	}

	PepIterator * TrypticIterator::operator++(int)
	{
		if (actual_pep_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		PepIterator * old = new TrypticIterator(*this);
		actual_pep_ = next_();
		if (f_iterator_->isAtEnd() && !hasNext_()) 
		{
			is_at_end_ = true;
		}
		return old;
	}	

	void TrypticIterator::setFastaFile(const String & f)
	{
		std::fstream fs;
		fs.open(f.c_str());
		if (!fs.is_open())
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, f);
		}
		f_file_ = f;
	}
	
	String TrypticIterator::getFastaFile()
	{
		return (f_file_);
	}

	bool TrypticIterator::begin ()
	{
		if (f_file_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		f_iterator_ = Factory<PepIterator>::create("FastaIterator");
		f_iterator_->setFastaFile(f_file_);
		if (! f_iterator_->begin())
		{
			return false;
			
		}
		
		f_entry_ = **f_iterator_;
		actual_pep_ = next_();
		
		return true;
	}

	std::string TrypticIterator::next_ ()
	{
		
		std::string seq = f_entry_.second;
		
		while (b_ < seq.length())
		{
			bool isInSpec = (e_==0||isDigestingEnd(seq[e_-1],seq[e_])||e_>=(seq.length()-1));
			if (isInSpec)
			{
				if (e_ < seq.length()){
					e_++;
					if (e_>(b_+1)) {
						if (e_==seq.length()) 
						{
							return (seq.substr(b_,e_-b_));
						} else 
						{
							return (seq.substr(b_,e_-b_-1));
						}
					}
				}
			} else {
				e_++;
			}
			if (e_>=seq.length())
			{
				goToNextAA_();
			}
		}
		
		if (f_iterator_->isAtEnd())
		{
			return "";
			
		}
		++(*f_iterator_);
		f_entry_ = **f_iterator_;
		b_ = 0;
		e_ = 0;
		return next_();
	}
	
	
	bool TrypticIterator::hasNext_ ()
	{
		unsigned int bold = b_;
		unsigned int eold = e_;
		std::string res = next_();
		b_ = bold;
		e_ = eold;
		if (res.length()==0)
		{
			return false;
		}
		return true;
	}
	
	void TrypticIterator::goToNextAA_ ()
	{
		std::string seq = f_entry_.second;
		b_++;
		while ((b_ < seq.length())&&!isDigestingEnd (seq[b_-1],seq[b_]))
		{
			b_++;
		}
		e_ = b_;
		
	}
	
	bool TrypticIterator::isAtEnd (){
		return is_at_end_;
	}
	
	bool TrypticIterator::isDigestingEnd(char aa1, char aa2)
	{
		return (aa1 == 'K' || aa1 == 'R') && aa2 != 'P';
	}
		
} // namespace OpenMS	
