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
#include <OpenMS/FILTERING/TRANSFORMERS/TradSeqQuality.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

#include <cmath>
//debug
#include <cassert>

using namespace std;

namespace OpenMS{

  TradSeqQuality::TradSeqQuality()
    :FilterFunctor()
  {
		setName(TradSeqQuality::getProductName());
    defaults_.setValue("xcorr_1+", 1.9, "min XCorr for charge state 1");
    defaults_.setValue("xcorr_2+", 2.2, "min XCorr for charge state 2");
    defaults_.setValue("xcorr_3+", 3.75, "min XCorr for charge state 3");

    defaults_.setValue("dCn_1+", 0.08, "min deltaCN for charge state 1");
    defaults_.setValue("dCn_2+", 0.08, "min deltaCN for charge state 2");
    defaults_.setValue("dCn_3+", 0.08, "min deltaCN for charge state 3");
		defaultsToParam_();
  }

  TradSeqQuality::TradSeqQuality(const TradSeqQuality& source)
    :	FilterFunctor(source)
  {
  }
    
  TradSeqQuality& TradSeqQuality::operator = (const TradSeqQuality& source)
  {
		if (this != &source)
		{
    	FilterFunctor::operator=(source);
		}
    return *this;
  }

  TradSeqQuality::~TradSeqQuality()
  {
  }

  /**
   \return a Quality measure based on the often used limits for sequest scores (XCorr and deltaCn) <br>
    if a ProteinIdentification is not fit ( no PeptideHit etc) -1000 is returned
   */
  double TradSeqQuality::operator() (const ClusterSpectrum& cspec)
  {
    //todo really necessary?
    cspec.getSpec();
    
    vector<double> xcorr_c(vector<double>(4));
    xcorr_c[1] = (double)param_.getValue("xcorr_1+");
    xcorr_c[2] = (double)param_.getValue("xcorr_2+");
    xcorr_c[3] = (double)param_.getValue("xcorr_3+");
    
    vector<double> dCn_c(vector<double>(4));
    dCn_c[1] = (double)param_.getValue("dCn_1+");
    dCn_c[2] = (double)param_.getValue("dCn_2+");
    dCn_c[3] = (double)param_.getValue("dCn_3+");
    
    uint charge = cspec.getParentionCharge();

    //vector<double> result;
    
    if (cspec.getPeptideIdentifications().size() == 0 || cspec.getPeptideIdentifications().begin()->getHits().size() == 0)
    {
      //result.push_back(-1000);
      return -1000; //result;
    }
    
    // sort PeptideHits
    const_cast<PeptideIdentification&>(*cspec.getPeptideIdentifications().begin()).sort();
    
    vector<PeptideHit>::const_iterator peph = cspec.getPeptideIdentifications().begin()->getHits().begin();
    
    double xcorr = peph->getScore();

    ++peph;
    
    double deltaCn;
    if ( peph != cspec.getPeptideIdentifications().begin()->getHits().end() ) 
    {
      deltaCn = (xcorr - peph->getScore()) / xcorr;
    }
    else 
    {
      deltaCn = 0;
    }
    
    if (xcorr > xcorr_c[charge] && deltaCn > dCn_c[charge])
    {
      return (xcorr-xcorr_c[charge] + deltaCn - dCn_c[charge]);
    }
    else
    {
      double xc = xcorr - xcorr_c[charge];
      double dc = deltaCn - dCn_c[charge];
      double res = 0;
      if ( xc < 0 ) res += xc;
      if ( dc < 0 ) res += dc;
      //result.push_back(res);
			return res;
    }
    
    return 0;
  }
}
