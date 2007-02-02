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
#include <OpenMS/FILTERING/TRANSFORMERS/KellerQuality.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <cmath>

//debug
#include <cassert>

using namespace std;

namespace OpenMS
{

  KellerQuality::KellerQuality(const KellerQuality& source)
    : FilterFunctor(source)
  {
  }
    
  KellerQuality& KellerQuality::operator = (const KellerQuality& source)
  {
		if (this != &source)
		{
    	FilterFunctor::operator=(source);
		}
    return *this;
  }
  
  KellerQuality::KellerQuality()
  {
		setName(KellerQuality::getProductName());
    // values from private communication from andrew keller

    defaults_.setValue("constant_1+", -0.236);
    defaults_.setValue("constant_2+", -1.498);
    defaults_.setValue("constant_3+", -1.975);

    defaults_.setValue("xcorr_1+", 8.346);
    defaults_.setValue("xcorr_2+", 9.3);
    defaults_.setValue("xcorr_3+", 10.685);

    defaults_.setValue("peplen_fac_1+", 2);
    defaults_.setValue("peplen_fac_2+", 1);
    defaults_.setValue("peplen_fac_3+", 1);
    
    defaults_.setValue("N_L_1+", 1);
    defaults_.setValue("N_L_2+", 2);
    defaults_.setValue("N_L_3+", 4);

    defaults_.setValue("N_C_1+", 1000000000); // very large number i.e. there is no N_C (MAX_DOUBLE etc make problems in FactoryProduct_test
    defaults_.setValue("N_C_2+", 15); 
    defaults_.setValue("N_C_3+", 25);

    defaults_.setValue("dCn_1+", 3.904);
    defaults_.setValue("dCn_2+", 7.317);
    defaults_.setValue("dCn_3+",11.263);
		defaultsToParam_();
  }

  KellerQuality::~KellerQuality()
  {
  }

  double KellerQuality::operator () (const ClusterSpectrum& cspec)
  {
    //todo really necessary?
    cspec.getSpec();
    
    vector<double> constant(vector<double>(4));
    constant[1] = (double)param_.getValue("constant_1+");
    constant[2] = (double)param_.getValue("constant_2+");
    constant[3] = (double)param_.getValue("constant_3+");
  
    vector<double> xcorr_c(vector<double>(4));
    xcorr_c[1] = (double)param_.getValue("xcorr_1+");
    xcorr_c[2] = (double)param_.getValue("xcorr_2+");
    xcorr_c[3] = (double)param_.getValue("xcorr_3+");
    
    vector<double> n_l(vector<double>(4));
    n_l[1] = (double)param_.getValue("N_L_1+");
    n_l[2] = (double)param_.getValue("N_L_2+");
    n_l[3] = (double)param_.getValue("N_L_3+");

    vector<double> peplen_fac(vector<double>(4));
    peplen_fac[1] = (double)param_.getValue("peplen_fac_1+");
    peplen_fac[2] = (double)param_.getValue("peplen_fac_2+");
    peplen_fac[3] = (double)param_.getValue("peplen_fac_3+");

    vector<double> n_c(vector<double>(4));
    n_c[1] = (double)param_.getValue("N_C_1+");
    n_c[2] = (double)param_.getValue("N_C_2+");
    n_c[3] = (double)param_.getValue("N_C_3+");

    vector<double> dCn_c(vector<double>(4));
    dCn_c[1] = (double)param_.getValue("dCn_1+");
    dCn_c[2] = (double)param_.getValue("dCn_2+");
    dCn_c[3] = (double)param_.getValue("dCn_3+");

    uint charge = cspec.getParentionCharge();

    //vector<double> result;
    
    // no identification, filter out
    if (cspec.getIdentification().size() == 0) 
    {
      //result.push_back(-1000);
      return -1000; 
    }
   
    if (cspec.getIdentification().begin()->getPeptideHits().size() == 0)
    {
      //result.push_back(-1000);
      return -1000; //result;
    }
    
    // sort PeptideHits
    const_cast<Identification&>(*cspec.getIdentification().begin()).sort();
    
    vector<PeptideHit>::const_iterator peph = cspec.getIdentification().begin()->getPeptideHits().begin();
    
    double xcorr = peph->getScore();
    uint peplen = peph->getSequence().size();

    ++peph;
    
    double deltaCn(0);
    if (peph != cspec.getIdentification().begin()->getPeptideHits().end())
    {
      deltaCn = (xcorr - peph->getScore()) / xcorr;
    }
    else 
    {
      deltaCn = 0;
    }
    
    double result = constant[charge] + 
										xcorr_c[charge]* (log(xcorr)/log(min(peplen * peplen_fac[charge],n_c[charge])) * n_l[charge]) +
      							dCn_c[charge] * deltaCn;
    return result;
    
    // sprank omitted since it is not in the DB and has a very low coefficient
  }
}
