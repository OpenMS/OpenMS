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
// $Id: TradSeqQuality.C,v 1.5 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/TradSeqQuality.h>

#include <cmath>
//debug
#include <cassert>

using namespace std;

namespace OpenMS{

  const String TradSeqQuality::info_ = "> 0 if XCorr and deltaCN > given values ";

  TradSeqQuality::TradSeqQuality()
    :FilterFunctor()
  {
		name_ = TradSeqQuality::getName();
    defaults_.setValue("xcorr_1+", 1.9);
    defaults_.setValue("xcorr_2+", 2.2);
    defaults_.setValue("xcorr_3+", 3.75);

    defaults_.setValue("dCn_1+", 0.08);
    defaults_.setValue("dCn_2+", 0.08);
    defaults_.setValue("dCn_3+", 0.08);
		param_ = defaults_;
  }

  TradSeqQuality::TradSeqQuality(const TradSeqQuality& source)
    :FilterFunctor(source)
  {
  }
    
  TradSeqQuality& TradSeqQuality::operator=(const TradSeqQuality& source )
  {
    FilterFunctor::operator=(source);
    return *this;
  }

  TradSeqQuality::~TradSeqQuality()
  {
  }

  String TradSeqQuality::info() const
  {
    return info_;
  }

  /**
   \return a Quality measure based on the often used limits for sequest scores (XCorr and deltaCn) <br>
    if a Identification is not fit ( no PeptideHit etc) -1000 is returned
   */
  vector<double> TradSeqQuality::operator() ( const ClusterSpectrum& cspec)
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

    vector<double> result;
    
    if (cspec.getIdentification().size() == 0 || cspec.getIdentification().begin()->getPeptideHits().size() == 0 )
    {
      result.push_back(-1000);
      return result;
    }
    
    // sort PeptideHits
    const_cast<Identification&>(*cspec.getIdentification().begin()).sort();
    
    vector<PeptideHit>::const_iterator peph = cspec.getIdentification().begin()->getPeptideHits().begin();
    
    double xcorr = peph->getScore();

    ++peph;
    
    double deltaCn;
    if ( peph != cspec.getIdentification().begin()->getPeptideHits().end() ) 
    {
      deltaCn = ( xcorr - peph->getScore() ) / xcorr;
    }
    else 
    {
      deltaCn = 0;
    }
    
    if ( xcorr > xcorr_c[charge] && deltaCn > dCn_c[charge] )
    {
      result.push_back(xcorr-xcorr_c[charge] + deltaCn - dCn_c[charge] );
    }
    else
    {
      double xc = xcorr - xcorr_c[charge];
      double dc = deltaCn - dCn_c[charge];
      double res = 0;
      if ( xc < 0 ) res += xc;
      if ( dc < 0 ) res += dc;
      result.push_back(res);
    }
    
    return result;
  }
}
