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
//
#include <OpenMS/FILTERING/TRANSFORMERS/PeakDiffBins.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

using namespace std;

namespace OpenMS
{
  PeakDiffBins::PeakDiffBins()
    : FilterFunctor(),
			mask_() 
  {
		check_defaults_ = false;
		setName(PeakDiffBins::getProductName());
		defaultsToParam_();

    // value from Bioinformatics, Bern 2004
    double mindiff = 1;
    double maxdiff = 187;
    double binsize = 1;
    vector<double> mask;
    for ( double i = mindiff; i < maxdiff ; i += binsize )
    {
      mask.push_back(i);
    }
    setmask(mask);
  }

  PeakDiffBins::PeakDiffBins(const PeakDiffBins& source)
    : FilterFunctor(source),
			mask_(source.mask_) 
  {
  }
  
  PeakDiffBins& PeakDiffBins::operator = (const PeakDiffBins& source)
  {
		if (this != &source)
		{
    	FilterFunctor::operator = (source);
    	mask_ = source.mask_;
		}
    return *this;
  }
  
  PeakDiffBins::~PeakDiffBins()
  {
  }

  /**
  this function violates the interface of FactoryProduct and should only used for testing
  \param newmask borders of new regions
  */
  void PeakDiffBins::setmask(std::vector<double>& newmask)
  {
    mask_.clear();
    for (uint i = 0; i < newmask.size(); ++i)
    {
      mask_.insert(make_pair(newmask[i],i-1));
    }
    //cerr << "mask ";
    //for ( map<double,int>::const_iterator cmit = mask_.begin(); cmit != mask_.end(); ++cmit )
    //{
    //  cerr << cmit->first << " ";
    //}
    //cerr << "\n";
  }
  
  vector<double> PeakDiffBins::operator() ( const ClusterSpectrum& cspec)
  {
    vector<double> result = vector<double>(mask_.size());
    double total = 0;
    //iterate over all peaks
    for (uint i = 0; i < cspec.getSpec().size(); ++i)
    {
      //look for each peakdifference that is in range of aa residuemasses (56/187), if it could be a aa (aamass)
      for (uint j = i; i+j < cspec.getSpec().size(); ++j)
      {
        double diff =  cspec.getSpec().getContainer()[i+j].getPosition()[0] - cspec.getSpec().getContainer()[i].getPosition()[0];
        total += cspec.getSpec().getContainer()[i+j].getIntensity() + cspec.getSpec().getContainer()[i].getIntensity();
        map<double,int>::const_iterator cmit = mask_.upper_bound(diff);
        if (cmit == mask_.begin() || cmit == mask_.end())
        {
          // nop
        }
        else
        {
          result[cmit->second]+= cspec.getSpec().getContainer()[i+j].getIntensity() + cspec.getSpec().getContainer()[i].getIntensity();
        }
      }
    }

    for (uint i = 0; i < result.size(); ++i)
    {
      result[i] /= total;
    }
    
    return result;
  }
}  
