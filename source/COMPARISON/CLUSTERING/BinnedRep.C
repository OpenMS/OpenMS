// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch  $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <cmath>

using namespace std;

namespace OpenMS
{

  /**
  standard values for bin dimensions: <br>
  size = 1<br>
  spread = 0 <br>
  */
  BinnedRep::BinnedRep() 
    :bins_(),binsize_(1),spread_(0),begin_(0),end_(0),id_(0),retention_(0),parent_m_z_(0),precursorpeakcharge_(0)
  {
  }

  BinnedRep::BinnedRep(double binsize,uint spread ) 
    :bins_(),binsize_(binsize),spread_(spread),begin_(0),end_(0),id_(0),retention_(0),parent_m_z_(0),precursorpeakcharge_(0)
  {
    if (binsize_ <= 0) 
    {
      binsize_ = 1;
      std::cerr << binsize << " is no sensible value for binsize!\n";
    }
  }

	BinnedRep::BinnedRep(const PeakSpectrum& spectrum, double binsize, uint spread)
		: bins_(),
			binsize_(binsize),
			spread_(spread),
			begin_(0),
			end_(0),
			id_(0),
			retention_(0),
			parent_m_z_(0),
			precursorpeakcharge_(0)
	{
		*this << spectrum;
	}
	
  BinnedRep::BinnedRep(const BinnedRep& source)
    :bins_(source.bins_), binsize_(source.binsize_),spread_(source.spread_), begin_(source.begin_), end_(source.end_), id_(source.id_), retention_(source.retention_), parent_m_z_(source.parent_m_z_),precursorpeakcharge_(source.precursorpeakcharge_)
  {
  }

  BinnedRep& BinnedRep::operator = (const BinnedRep& source)
  {
    bins_=source.bins_;
    binsize_=source.binsize_;
    begin_=source.begin_;
    end_=source.end_;
    id_=source.id_;
    retention_=source.retention_;
    parent_m_z_=source.parent_m_z_;
    precursorpeakcharge_ = source.precursorpeakcharge_;
    spread_=source.spread_;
    return *this;
  }

  BinnedRep::~BinnedRep()
  {
  }

  double BinnedRep::operator[] (int index) const
  {
    return bins_.at(index);
  }

  String BinnedRep::str () const
  {
    stringstream ss;
    ss << "# id: " << id_ << " retention time (normalized): " << retention_ 
       << " parent mass: " << parent_m_z_ << " size: " << bins_.size() 
       << " binsize: " << binsize_ << " begin: " << begin_ << " end: " << end_ << "\n";
    for ( BinnedRep::const_iterator it = bins_.begin(); it != bins_.end(); it.hop())
    {
      ss << it.position()*binsize_ + begin_ << "\t" << *it << "\n";
    }
    return ss.str();
  }

  /**
  required for BinnedRepWidget <br>
  */
  void BinnedRep::normalize()
  {
    double max = 0;
    for ( SparseVector<int>::iterator bit = bins_.begin(); bit != bins_.end(); bit.hop())
    {
      if ( (int)*bit > max )
      {
        max = (int)*bit;
      }
    }
    for ( SparseVector<int>::iterator bit = bins_.begin(); bit != bins_.end(); bit.hop())
    {
      *bit = (int)*bit/max;
    }
  }

  void BinnedRep::clear_()
  {
    id_ = 0;
    retention_ = 0;
    parent_m_z_ = 0;
    begin_ = 0;
    end_ = 0;
    bins_.clear();
  }

 	void operator << (BinnedRep& binrep, const PeakSpectrum& spectrum)
  {
    binrep.clear_();
    //set information from PeakSpectrum TODO
    binrep.id_ = spectrum.getPersistenceId();
    binrep.retention_ = spectrum.getRT();
    binrep.parent_m_z_ = spectrum.getPrecursorPeak().getPosition()[0];
    binrep.precursorpeakcharge_ = spectrum.getPrecursorPeak().getCharge();
    
    if (spectrum.getContainer().size())
    {
      //calculate size of needed container 
      binrep.begin_ = ((int)(spectrum.begin()->getMZ()/binrep.binsize_)) * binrep.binsize_;
      binrep.end_ =((int)((spectrum.end()-1)->getMZ()/binrep.binsize_)) * binrep.binsize_;
      int size = (int) ((binrep.end_ - binrep.begin_) / binrep.binsize_ + 3 ); // +3 is safety distance for rounding errors
      binrep.bins_.resize(size);
      
      //iterate over Peaks and throw them into bins
      for (PeakSpectrum::ConstIterator cit = spectrum.begin(); cit != spectrum.end(); ++cit)
      {
        uint position = (uint) ((cit->getMZ()-binrep.begin_) / binrep.binsize_);
        try 
				{
          binrep.bins_[position] =  binrep.bins_.at(position) + cit->getIntensity();
        }
        catch (out_of_range&)
        {
          cerr <<  "out of range in BinnedRep, this should not happen\n"<< endl;
        }

        //insert into neighbouring bins according to spread
        if (binrep.spread_ > 0 )
        {
          for (uint i = 1; i <= binrep.spread_ && position + i < binrep.size(); ++i)
          {
            binrep.bins_[position+i] = binrep.bins_.at(position+i)+ cit->getIntensity();
          }
          for (uint i = 1; i <= binrep.spread_ && position >= i ; ++i)
          {
            binrep.bins_[position-i]  = binrep.bins_.at(position -i)+ cit->getIntensity();
          }
        }
      }
      binrep.normalize(); //some other functions depend on normalized binreps
    }
    return;
  }
} 

