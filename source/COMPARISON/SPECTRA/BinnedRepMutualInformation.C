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
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepMutualInformation.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

#include <cmath> 

using namespace std;

namespace OpenMS
{
  BinnedRepMutualInformation::BinnedRepMutualInformation()
		:	BinnedRepCompareFunctor()
  {
		setName(BinnedRepMutualInformation::getProductName());
    defaults_.setValue("intervals", 3);
		defaultsToParam_();
  }

  BinnedRepMutualInformation::BinnedRepMutualInformation(const BinnedRepMutualInformation& source)
    : BinnedRepCompareFunctor(source)
  {
  }

  BinnedRepMutualInformation& BinnedRepMutualInformation::operator = (const BinnedRepMutualInformation& source)
  {
		if (this != &source)
		{
    	BinnedRepCompareFunctor::operator = (source);
		}
    return *this;
  }

  BinnedRepMutualInformation::~BinnedRepMutualInformation()
  {
  }

	double BinnedRepMutualInformation::operator () (const BinnedRep& a) const
	{
		return operator () (a, a);
	}
	
  double BinnedRepMutualInformation::operator () (const BinnedRep& a, const BinnedRep& b) const
  {
    uint intervals = (unsigned int)param_.getValue("intervals");
   
    //const BinnedRep& b = csb.getBinrep();
    //const BinnedRep& a = csa.getBinrep();
    
    //double filterfactor = filter(csa,csb);
    //if (filterfactor < 1e-12) return 0;
    //number of pairs where a falls into interval i and b falls into interval j
    vector<vector<double> > n = vector<vector<double> >(intervals, vector<double>(intervals, 0.0));
    //marginal frequencies for x and y
    double result = 0;
    vector<double> na = vector<double>(intervals, 0.0);
    vector<double> nb = vector<double>(intervals, 0.0);
    int nab = 0; //number of pairs where signals are present in both binreps
     
    BinnedRep::const_iterator bit = b.begin(); 
    BinnedRep::const_iterator ait = a.begin();
    while (ait != a.end() && bit != b.end())
    {
      //we are at the same position (+- precision)
      if ((ait.position() * a.getBinSize() + a.min()) == (bit.position() * b.getBinSize() + b.min()))
      {
        uint int_a = (int)floor(*ait * intervals + 0.5);
        if (int_a == intervals) int_a--;
        uint int_b = (int)floor(*bit * intervals + 0.5);
        if (int_b == intervals) int_b--;
        if (*ait > 0 || *bit > 0)
        {
          n.at(int_a).at(int_b) += 1;
          nab++;
        }
        ait.hop();
        bit.hop();

      }
      //ait lags
      else if (((ait.position() * a.getBinSize() + a.min()) - (bit.position() * b.getBinSize() + b.min())) < 0)
      {
        uint int_a =(int)floor(*ait * intervals + 0.5);
        if (int_a == intervals) int_a--;
        n.at(int_a).at(0) += 1;
        nab++;
        ait.hop();
      }
      //bit lags
      else if (((ait.position() * a.getBinSize() + a.min()) - (bit.position() * b.getBinSize() + b.min())) > 0)
      {
        uint int_b = (int)floor(*bit * intervals + 0.5);
        if (int_b == intervals) int_b--;
        n.at(0).at(int_b) += 1;
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

    while (bit != b.end()) 
    {
      uint int_b = (int)floor(*bit * intervals + 0.5);
      if (int_b == intervals) int_b--;

      n.at(0).at(int_b) += 1;
      nab++;
      bit.hop();
    }
		
    while (ait != a.end())
    {
      uint int_a = (int)floor(*ait * intervals + 0.5);
      if (int_a == intervals) int_a--;
			
      n.at(int_a).at(0) += 1;
      nab++;
      ait.hop();
    }
	
    for (uint i = 0; i < intervals; ++i)
    {
      double marg_freq_a = 0;
      for (uint j = 0; j < intervals; ++j)
      {
        marg_freq_a += n[i][j];
      }
      na[i] = marg_freq_a/(double)nab;
    }
		
    for (uint j = 0; j < intervals; ++j)
    {
      double marg_freq_b = 0;
      for (uint i = 0; i < intervals; ++i)
      {
        marg_freq_b += n[i][j];
      }
      nb[j] = marg_freq_b/(double)nab;
    }
		
    for (uint i  = 0; i < intervals; ++i)
    {
      for (uint j = 0; j < intervals; ++j)
      {
        double tempresult = (n[i][j]/(double)nab) * log((n[i][j]/(double)nab) / (na[i] * nb[j])) / log(2.0f);
        if (fabs(n[i][j]) > 0) 
				{
					result += tempresult;
				}
      }
    }
    return result;
  }

}
