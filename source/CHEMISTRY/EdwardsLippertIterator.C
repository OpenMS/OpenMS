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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------
 
#include <OpenMS/CHEMISTRY/EdwardsLippertIterator.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <cstring>
#include <algorithm>
#include <fstream>

using namespace std;

namespace OpenMS
{

/**
@brief comperator for two DoubleReals with a tolerance value
*/
struct FloatsWithTolLess : public binary_function<DoubleReal , DoubleReal, bool>
{
	/**
	@brief constructor
	@param t const reference to the tolerance
	*/
	FloatsWithTolLess(const DoubleReal & t) : tol_(t) {}
	/**
	@brief copy constructor
	*/
	FloatsWithTolLess(const FloatsWithTolLess & rhs ) : tol_(rhs.tol_) {}
	
	/**
	@brief implementation of the '<' operator for two DoubleReals with the tolerance value
	@param f1 first DoubleReal
	@param f2 second DoubleReal
	@return true if first DoubleReal '<' second DoubleReal-tolerance
	*/
	bool operator()( DoubleReal f1, DoubleReal f2) const
	
	{
		return (f1<(f2-tol_));
	}

	protected:
	DoubleReal const & tol_; ///< tolerance value
};


	typedef pair<String,String> FASTAEntry;
	
	///Constructor to intialize algorithm
	EdwardsLippertIterator::EdwardsLippertIterator ()
		: PepIterator()
	{
		ResidueDB* rdb = ResidueDB::getInstance();
		
		char aa[] = "ARNDCEQGHILKMFPSTWYV";
		
		for (Size i = 0; i<255;++i)
		{
			masse_[i]=0;
		}
	
		for (Size i = 0; i<strlen(aa);++i)
		{
			const Residue * r = rdb->getResidue(aa[i]);
			masse_[(int)aa[i]]=r->getAverageWeight();
		}
		
		b_ = 0;
		e_ = 0;
		m_ = 0;
		tol_ = 0.5;
		is_at_end_ = false;
		
	}

	///destructor
	EdwardsLippertIterator::~EdwardsLippertIterator ()
	{
	
	}
	
	EdwardsLippertIterator::EdwardsLippertIterator (const EdwardsLippertIterator & source) :
		PepIterator(source),
		f_file_(source.f_file_),
		actual_pep_(source.actual_pep_),
		spec_(source.spec_),
		tol_(source.tol_),
		is_at_end_(source.is_at_end_),
		f_iterator_(source.f_iterator_),
		f_entry_(source.f_entry_),
		b_(source.b_),
		e_(source.e_),
		m_(source.m_),
		massMax_(source.massMax_)
	{
		for (Size i = 0; i < 256;i++)
		{
			masse_[i]=source.masse_[i];
		}
	}

	FASTAEntry EdwardsLippertIterator::operator*()
	{
		if (actual_pep_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		return FASTAEntry(f_entry_.first,actual_pep_);
	}
	
	PepIterator & EdwardsLippertIterator::operator++()
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

	PepIterator * EdwardsLippertIterator::operator++(int)
	{
		if (actual_pep_=="")
		{
			throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		PepIterator * old = new EdwardsLippertIterator(*this);
		actual_pep_ = next_();
		if (f_iterator_->isAtEnd() && !hasNext_()) 
		{
			is_at_end_ = true;
		}
		return old;
	}	

	void EdwardsLippertIterator::setTolerance (DoubleReal t)
	{
		if (t<0)
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"tolerace must not be negative",(String)t);
		}
		tol_ = t;
	}

	DoubleReal EdwardsLippertIterator::getTolerance ()
	{
		return (tol_);
	}
	
	void EdwardsLippertIterator::setSpectrum (const vector<DoubleReal> & s)
	{
		//check if spectrum is sorted
		for (Size i = 1; i < s.size();++i)
		{
			if (s.at(i-1)>s.at(i))
			{
				throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Spectrum has to be sorted ascendingly","");
			}
		}
		spec_ = s;
		massMax_ = spec_.back();
	}

	const vector<DoubleReal> & EdwardsLippertIterator::getSpectrum ()
	{
		return (spec_);
	}
	
	void EdwardsLippertIterator::setFastaFile(const String & f)
	{
		fstream fs;
		fs.open(f.c_str());
		if (!fs.is_open())
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, f);
		}
		f_file_ = f;
	}
	
	String EdwardsLippertIterator::getFastaFile()
	{
		return (f_file_);
	}

	bool EdwardsLippertIterator::begin ()
	{
		if (f_file_=="" || spec_.size()==0)
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

	string EdwardsLippertIterator::next_ ()
	{
		
		string seq = f_entry_.second;
		
		while (b_ < seq.length())
		{
			
			bool isInSpec = (e_==0 ||e_>=(seq.length()-1)||isDigestingEnd(seq[e_-1],seq[e_])) ? isInSpectrum_(m_) : false;
			//cout<<"b:"<<b_<< " e:"<<e_<< " m:"<<m_<< " seql:"<<seq.length()<<endl;
			if (isInSpec && m_<=massMax_+tol_)
			{
				if (e_<=seq.length()) {
					if (e_<seq.length()) m_+=masse_[(int)seq[e_]];
					e_++;
					if (e_>(b_+1)) {
						return (seq.substr(b_,e_-b_-1));
					}
				}
				
			} else {
				if (m_> massMax_+tol_){
					goToNextAA_();
				} else {
					m_+=masse_[(int)seq[e_]];
					e_++;
				}
			}
			if (e_>seq.length())
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
		m_ = 0;
		return next_();
	}
	
	
	bool EdwardsLippertIterator::hasNext_ ()
	{
		DoubleReal mold = m_;
		unsigned int bold = b_;
		unsigned int eold = e_;
		string res = next_();
		m_ = mold;
		b_ = bold;
		e_ = eold;
		if (res.length()==0)
		{
			return false;
		}
		return true;
	}
	
	void EdwardsLippertIterator::goToNextAA_ ()
	{
		string seq = f_entry_.second;
		m_ = 0;
		b_++;
		
		while ((b_ < seq.length())&&!isDigestingEnd (seq[b_-1],seq[b_]) )
		{
			b_++;
		}
		e_ = b_;
		
	}
	
	bool EdwardsLippertIterator::isAtEnd (){
		return is_at_end_;
	}
	
	bool EdwardsLippertIterator::isDigestingEnd(char, char)
	{
		return true;
	}
		
	
	bool EdwardsLippertIterator::isInSpectrum_ (DoubleReal & mass)
	{
		return (binary_search(spec_.begin(),spec_.end(),mass,FloatsWithTolLess(tol_)));
	}

} // namespace OpenMS

