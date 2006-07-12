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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepMutualInformation.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

#include <cmath> 
#include <exception>

using namespace std;

namespace OpenMS
{
  BinnedRepMutualInformation::BinnedRepMutualInformation()
  {
		name_ = BinnedRepMutualInformation::getName();
    defaults_.setValue("intervals", 3);
    usebins_ = true;
		param_ = defaults_;
  }

  BinnedRepMutualInformation::BinnedRepMutualInformation(const BinnedRepMutualInformation& source)
    : CompareFunctor(source)
  {
  }

  BinnedRepMutualInformation& BinnedRepMutualInformation::operator = (const BinnedRepMutualInformation& source)
  {
    CompareFunctor::operator = (source);
		usebins_ = source.usebins_;
    return *this;
  }

  BinnedRepMutualInformation::~BinnedRepMutualInformation()
  {
  }

  double BinnedRepMutualInformation::operator()(const ClusterSpectrum& csa, const ClusterSpectrum& csb)const
  {
    uint intervals = (unsigned int)param_.getValue("intervals");
   
    const BinnedRep& b = csb.getBinrep();
    const BinnedRep& a = csa.getBinrep();
    
    double filterfactor = filter(csa,csb);
    if ( filterfactor < 1e-12) return 0;
    //number of pairs where a falls into interval i and b falls into interval j
    vector<vector<double> > n = vector<vector<double> >(intervals,vector<double>(intervals));
    //marginal frequencies for x and y
    double result = 0;
    vector<double> na = vector<double>(intervals); 
    vector<double> nb = vector<double>(intervals);
    int nab = 0; //number of pairs where signals are present in both binreps
     
    BinnedRep::const_iterator bit = b.begin(); 
    BinnedRep::const_iterator ait = a.begin();
    while (ait != a.end() && bit != b.end())
    {
      //we are at the same position (+- precision)
      if ( fabs( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  < 1e-8)
      {
        uint int_a = (uint)(*ait * intervals);
        if ( int_a == intervals ) int_a--;
        uint int_b = (uint)(*bit * intervals);
        if ( int_b == intervals ) int_b--;
        if( ( *ait > 1e-8 || *bit > 1e-8 ) )
        {
          n.at(int_a).at(int_b)++;
          nab++;
        }
        ait.hop();
        bit.hop();
      }
      //ait lags
      else if ( ( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  < 0 )
      {
        uint int_a = (uint)(*ait * intervals);
        if ( int_a == intervals ) int_a--;
        n.at(int_a).at(0)++;
        nab++;
        ait.hop();
      }
      //bit lags
      else if ( ( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  > 0 )
      {
        uint int_b = (uint)(*bit * intervals);
        if ( int_b == intervals ) int_b--;
        n.at(0).at(int_b)++;
        nab++;
        bit.hop();
      }
      //shouldnt happen 
      else 
      {
        cerr << "precision problems in BinnedRepMutualInformation.C\n"
          << ait.position()*a.getBinSize() + a.min() << " - "
          << bit.position()*b.getBinSize() + b.min() 
          << " is neither < 0 nor > 0 , but it is not between 1e-8 and -1e-8\n"; 
      }
    }
    while ( bit != b.end() ) 
    {
      uint int_b = (uint)(*bit * intervals);
      if ( int_b == intervals ) int_b--;
      n.at(0).at(int_b)++;
      nab++;
      bit.hop();
    }
    while ( ait != a.end() )
    {
      uint int_a = (uint)(*ait * intervals);
      if ( int_a == intervals ) int_a--;
      n.at(int_a).at(0)++;
      nab++;
      ait.hop();
    }
    // if they have no matching pairs, they are uncorrelated
    if ( nab == 0 ) 
    {
      return 0;
    }
    
    for (uint i = 0; i < intervals; ++i)
    {
      double marg_freq_a = 0;
      for (uint j = 0; j < intervals ; ++j)
      {
        marg_freq_a += n[i][j];
      }
      na[i] = marg_freq_a/nab;
    }
    for (uint j = 0; j < intervals; ++j)
    {
      double marg_freq_b = 0;
      for ( uint i = 0; i < intervals; ++i)
      {
        marg_freq_b += n[i][j];
      }
      nb[j] = marg_freq_b/nab;
    }
    for (uint i  = 0; i < intervals; ++i)
    {
      for (uint j = 0; j < intervals; ++j)
      {
        double tempresult = ( n[i][j]/nab ) * log( ( n[i][j]/nab ) / (na[i]*nb[j]) ) / log(2.0f);
        if ( fabs(n[i][j]) > 1e-8 ) result += tempresult;
      }
    }
    return result * filterfactor;
  }

}
