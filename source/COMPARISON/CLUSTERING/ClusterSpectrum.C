// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/FORMAT/DBAdapter.h>

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
    	adapterp_(0),
    	binsize_(0),
    	binspread_(0),
    	id_(0),
    	cached_(0),
    	retention_(0),
    	parent_mass_(0),
    	parentioncharge_(0)
  {
  }
  
  ClusterSpectrum::ClusterSpectrum(long id, DBAdapter* adapterp,double size,uint spread)
    : specp_(0),
    	binrepp_(0),
    	adapterp_(adapterp), 
    	binsize_(size), 
    	binspread_(spread), 
    	id_(id),
    	cached_(0),
    	retention_(0),
    	parent_mass_(0),
    	parentioncharge_(0)
  {
  }

  ClusterSpectrum::ClusterSpectrum(const MSSpectrum< DPeak<1> >& spec,DBAdapter* adapterp, double binsize , uint binspread )
    : specp_(new MSSpectrum< DPeak<1> >(spec)), 
    	binrepp_(0),
    	adapterp_(adapterp),
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
  
  
  ClusterSpectrum::ClusterSpectrum(MSSpectrum< DPeak<1> >* specp, DBAdapter* adapterp, double binsize, uint binspread)
    : specp_(specp),
    	binrepp_(0),
    	adapterp_(adapterp), 
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

  ClusterSpectrum::ClusterSpectrum(BinnedRep* binrepp, DBAdapter* adapterp)
    : specp_(0), 
    	binrepp_(binrepp), 
    	adapterp_(adapterp), 
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

  ClusterSpectrum::ClusterSpectrum(MSSpectrum< DPeak<1> >* specp, BinnedRep* binrepp)
    : specp_(specp),
    	binrepp_(binrepp),
    	adapterp_(0),
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
    	adapterp_(source.adapterp_),
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
      specp_ = new MSSpectrum< DPeak<1> >(*source.specp_);
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
      specp_ = new MSSpectrum< DPeak<1> >(*source.specp_);
    }
    if ( source.binrepp_ )
    {
      binrepp_ = new BinnedRep(*source.binrepp_);
    }
    adapterp_=source.adapterp_;   
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
    if (!adapterp_)
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
    //if a DBAdapter is present, the other Representation is created
    else 
    {
      if (binrepp_)
      {
        return *binrepp_;  
      }
      else
      {
        if (binsize_ < 0 )
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"binsize_");
        }
        else 
        {
          if (!specp_)
          {
            getSpec();
          }
          binrepp_ = new BinnedRep(binsize_,binspread_);
          (*binrepp_) << (*specp_);
          // this works if only one ( i.e. binrepp_ or specp_ is needed )
          // otherwise it takes much time
          delete specp_;
          specp_ = 0;
          return *binrepp_;
        }
      }
    }
  }

  const MSSpectrum< DPeak<1> >& ClusterSpectrum::getSpec() const
  {
    return const_cast<ClusterSpectrum*>(this)->spec();
  }
  
  /** 
  do not use with python!
  */
  MSSpectrum< DPeak<1> >& ClusterSpectrum::spec()
  {
    // in case spec is changed, binrep needs to be recreated
    if ( binrepp_ )
    {
      delete binrepp_;
      binrepp_ = 0;
    }
    if (!adapterp_)
    {
      if (specp_ )
      {
        return *specp_;
      }
      else 
      {
        throw WrongRepresentation(__FILE__, __LINE__, __PRETTY_FUNCTION__,"no MSSpectrum< DPeak<1> > in ClusterSpectrum and no DBAdapter* given");
      }
    }
    else 
    {
      if (specp_)
      {
        return *specp_;
      }
      else
      {
        //specp_ = dynamic_cast<MSSpectrum< DPeak<1> >*>(adapterp_->createObject(id_)); //TODO Persistence
        
        specp_->getContainer().sortByPosition();
        if (specp_)
        {
          return *specp_;
        }
        else 
        {
          ostringstream ss;
          ss.str("");
          ss << "id " << id_ << " ";
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,ss.str().c_str());
        }
      }
    }
  }

  const std::vector<Identification>& ClusterSpectrum::getIdentification() const
  {
    //todo 

    getSpec();

/*
    if ( specp_->getIdentification().size() == 0 )
    {
      if (!adapterp_ && ! adapterp )
      {
        cerr << (int) adapterp_ << endl;
        throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"No DBAdapter","this operation needs a DBAdapter*");
      }
      if (!adapterp) adapterp = adapterp_;
      if (specp_->getIdentification().size() == 0 )
      {
				#ifdef OPENMS_HAS_DB
        specp_->loadIdentification(dynamic_cast<DBAdapter*>(adapterp));
				#endif
      }
    }
*/
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
    else if ( adapterp_ )
    {
      getSpec();
      updatecache_();
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
