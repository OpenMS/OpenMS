// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

#include <sstream>

using namespace std;

namespace OpenMS
{

  ClusterSpectrum::DifferentSpectra::DifferentSpectra(const char* file, int line, const char* function) throw()
    : Base(file, line, function)
  {
    what_ = "The Spectra are incompatible\nif thrown by the ctor the ids probably are not equal\nif thrown by a CompareFunctor the bins are probably not compatible";

    Exception::globalHandler.setMessage(what_);
  }

  ClusterSpectrum::DifferentSpectra::~DifferentSpectra() throw()
  {
  }

  ClusterSpectrum::WrongRepresentation::WrongRepresentation(const char* file, int line, const char* function, const char* message ) throw() 
    : Base(file, line, function, "ClusterSpectrum::WrongRepresentation",message)
  {
  }

  ClusterSpectrum::WrongRepresentation::~WrongRepresentation() throw()
  {
  }

  // intended for the use in Containers, like the stl
  ClusterSpectrum::ClusterSpectrum()
    : specp_(0),
    	binrepp_(0),
    	binsize_(0),
    	binspread_(0),
    	id_(0),
    	cached_(0),
    	retention_(0),
    	parent_mass_(0),
    	parentioncharge_(0)
  {
  }
  
  ClusterSpectrum::ClusterSpectrum(long id,double size,uint spread)
    : specp_(0),
    	binrepp_(0),
    	binsize_(size), 
    	binspread_(spread), 
    	id_(id),
    	cached_(0),
    	retention_(0),
    	parent_mass_(0),
    	parentioncharge_(0)
  {
  }

  ClusterSpectrum::ClusterSpectrum(const PeakSpectrum& spec, double binsize , uint binspread )
    : specp_(new PeakSpectrum(spec)), 
    	binrepp_(0),
    	binsize_(binsize),
    	binspread_(binspread),
    	id_(spec.getPersistenceId()),
    	cached_(0),
    	retention_(0),
    	parent_mass_(0),
    	parentioncharge_(0)
  {
    updatecache_();
  }
  
  
  ClusterSpectrum::ClusterSpectrum(PeakSpectrum* specp, double binsize, uint binspread)
    : specp_(specp),
    	binrepp_(0),
    	binsize_(binsize), 
    	binspread_(binspread), 
    	id_(specp->getPersistenceId()),
    	cached_(0),
    	retention_(0),
    	parent_mass_(0),
    	parentioncharge_(0)
  {
    updatecache_();
  }

  ClusterSpectrum::ClusterSpectrum(BinnedRep* binrepp)
    : specp_(0), 
    	binrepp_(binrepp), 
    	binsize_(binrepp->getBinSize()), 
    	binspread_(binrepp->getBinSpread()),
    	id_(binrepp->id()), 
    	cached_(0),
    	retention_(0),
    	parent_mass_(0),
    	parentioncharge_(0)
  {
    updatecache_();
  }

  ClusterSpectrum::ClusterSpectrum(PeakSpectrum* specp, BinnedRep* binrepp)
    : specp_(specp),
    	binrepp_(binrepp),
    	binsize_(binrepp->getBinSize()), 
    	binspread_(binrepp->getBinSpread()),
    	id_(binrepp->id()),
    	cached_(0),
    	retention_(0),
    	parent_mass_(0),
    	parentioncharge_(0)
  {
    if (specp->getPersistenceId() != binrepp->id())
    {
      throw DifferentSpectra(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    updatecache_();
  }

  ClusterSpectrum::ClusterSpectrum( const ClusterSpectrum& source)
    : specp_(0),
    	binrepp_(0),
    	binsize_(source.binsize_),
    	binspread_(source.binspread_),
    	id_(source.id_),
    	cached_(source.cached_),
    	retention_(source.retention_),
    	parent_mass_(source.parent_mass_),
    	parentioncharge_(source.parentioncharge_)
  {
    if ( source.specp_ )
    {
      specp_ = new PeakSpectrum(*source.specp_);
    }
    if ( source.binrepp_ )
    {
      binrepp_ = new BinnedRep(*source.binrepp_);
    }
  }

  ClusterSpectrum& ClusterSpectrum::operator=(const ClusterSpectrum& source)
  {
    specp_ = 0;
    binrepp_ = 0;
    if ( source.specp_ )
    {
      specp_ = new PeakSpectrum(*source.specp_);
    }
    if ( source.binrepp_ )
    {
      binrepp_ = new BinnedRep(*source.binrepp_);
    }
    binsize_=source.binsize_;   
    binspread_=source.binspread_;   
    id_=source.id_;
    cached_=source.cached_;
    retention_=source.retention_;
    parent_mass_=source.parent_mass_;
    parentioncharge_=source.parentioncharge_;
    return *this;
  }
  
  ClusterSpectrum::~ClusterSpectrum()
  {
    delete binrepp_;
    delete specp_;
  }

  int ClusterSpectrum::id() const
  {
    return id_;
  }

  void ClusterSpectrum::strip() const
  {
    delete binrepp_;
    binrepp_ = 0;
    delete specp_;
    specp_ = 0;
  }
  
  const BinnedRep& ClusterSpectrum::getBinrep() const
  {
    if (binrepp_ )
    {
      return *binrepp_;
    }
    else if ( specp_ )
    {
      binrepp_ = new BinnedRep(binsize_,binspread_);
      (*binrepp_) << (*specp_);
      return *binrepp_;
    }
    else
    {
      throw WrongRepresentation(__FILE__, __LINE__, __PRETTY_FUNCTION__,"no BinnedRep in ClusterSpectrum and no DBAdapter* given");
    }
  }

  const PeakSpectrum& ClusterSpectrum::getSpec() const
  {
    return const_cast<ClusterSpectrum*>(this)->spec();
  }
  
  /** 
  do not use with python!
  */
  PeakSpectrum& ClusterSpectrum::spec()
  {
    // in case spec is changed, binrep needs to be recreated
    if ( binrepp_ )
    {
      delete binrepp_;
      binrepp_ = 0;
    }

    if (specp_ )
    {
      return *specp_;
    }
    else 
    {
      throw WrongRepresentation(__FILE__, __LINE__, __PRETTY_FUNCTION__,"no PeakSpectrum in ClusterSpectrum and no DBAdapter* given");
    }
  }

  const std::vector<Identification>& ClusterSpectrum::getIdentification() const
  {
    //todo 

    getSpec();

    return specp_->getIdentifications();
  }
  
  //todo not usable without prior spec() or binrep if created from id
  const uint& ClusterSpectrum::getParentionCharge() const
  {
    if ( !cached_ ) updatecache_();
    return parentioncharge_;
  }
  
  const double& ClusterSpectrum::getParentMass() const
  {
    if ( !cached_ ) updatecache_();
    return parent_mass_;
  }

  const double& ClusterSpectrum::getRetention() const
  {
    if ( !cached_ ) updatecache_();
    return retention_;
  }
  
  void ClusterSpectrum::updatecache_() const
  {
    if (specp_)
    {
      retention_ = specp_->getRetentionTime();
      parent_mass_ = specp_->getPrecursorPeak().getPosition()[0];
      parentioncharge_ = specp_->getPrecursorPeak().getCharge();
      cached_ = 1;     
    }
    else if ( binrepp_ )
    {
      retention_ = binrepp_->getRetention();
      parent_mass_ = binrepp_->getParentmz();
      parentioncharge_ = binrepp_->getPrecursorPeakCharge();
      cached_ = 1;     
    }
    else
    {
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"empty ClusterSpectrum","");
    }
  }

  PeptideHit ClusterSpectrum::getTophit() const
  {
    const vector<Identification>& dbsl = getIdentification();
    if ( dbsl.size() > 0 && dbsl.begin()->getPeptideHits().size() > 0 )
    {
    	cout << "Inside"<< endl;
      return *dbsl.begin()->getPeptideHits().begin();
    }
    return PeptideHit();
  }

  void ClusterSpectrum::clearcache_() const 
  {
    cached_ = 0;
  }
}
