// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <stdio.h>

namespace OpenMS
{

	InternalCalibration::InternalCalibration()
		:DefaultParamHandler("InternalCalibration"),
		 ProgressLogger()
	{
		defaults_.setValue("mz_tolerance",1.,"Allowed tolerance between peak and reference m/z.");
		defaults_.setMinFloat("mz_tolerance",0.);
		defaults_.setValue("mz_tolerance_unit","Da","Unit for mz_tolerance.");
		defaults_.setValidStrings("mz_tolerance_unit",StringList::create("Da,ppm"));
		defaults_.setValue("hires:percentage",30,"Percentage of spectra a signal has to appear in to be considered as background signal.");
		defaultsToParam_();
	}
	
  InternalCalibration::InternalCalibration(InternalCalibration& obj)
		: DefaultParamHandler(obj),
			ProgressLogger(obj)
  {}
  
  InternalCalibration& InternalCalibration::operator=(const InternalCalibration& obj)
  {
		// take care of self assignments
    if (this == &obj)		return *this;
		DefaultParamHandler::operator=(obj);
    return *this;
  }

	void InternalCalibration::checkReferenceIds_(std::vector<PeptideIdentification>& pep_ids)
	{
		 for(Size p_id = 0; p_id < pep_ids.size();++p_id)
		 {
				 if(pep_ids[p_id].getHits().size() > 1)
				 {
						 throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "InternalCalibration: Your Id-file contains PeptideIdentifications with more than one hit, use the IDFilter to select only the best hits.");
				 }
				 if(!pep_ids[p_id].metaValueExists("RT"))
				 {
						 throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "InternalCalibration: meta data value 'RT' missing for peptide identification!"); 
				 } 
				 if(!pep_ids[p_id].metaValueExists("MZ"))
				 {
						 throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "InternalCalibration: meta data value 'MZ' missing for peptide identification!"); 
				 }
		 }
	}	 


	void InternalCalibration::makeLinearRegression_(std::vector<DoubleReal>& observed_masses, std::vector<DoubleReal>& theoretical_masses)
	{
			if(observed_masses.size() != theoretical_masses.size())
			{
					throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Number of observed and theoretical masses must agree."); 
			}
#ifdef DEBUG_CALIBRATION
			std::ofstream out("calibration_regression.txt");
#endif
			std::vector<DoubleReal> rel_errors(observed_masses.size(),0.);
			// determine rel error in ppm for the two reference masses
			for(Size ref_peak=0; ref_peak < observed_masses.size();++ref_peak)
			{
				rel_errors[ref_peak] = (theoretical_masses[ref_peak]-observed_masses[ref_peak])/theoretical_masses[ref_peak] * 1e6;
#ifdef DEBUG_CALIBRATION
				out << observed_masses[ref_peak] << "\t"<< rel_errors[ref_peak] << "\n";
#endif
			}

			DoubleReal cov00, cov01, cov11, sumsq, slope,intercept;
			// TODO: what exactly is stride?? used 1 here as in the gsl-example :)
			gsl_fit_linear (&(observed_masses[0]), 1, &(rel_errors[0]), 1, observed_masses.size(), &intercept,&slope,&cov00,&cov01,&cov11,&sumsq);
			

			trafo_.setName("linear");
			trafo_.setParam("slope",slope);
			trafo_.setParam("intercept",intercept);

//#ifdef DEBUG_CALIBRATION
  	  printf ("# best fit: Y = %g + %g X\n", intercept, slope);
      printf ("# covariance matrix:\n");
      printf ("# [ %g, %g\n#   %g, %g]\n", 
               cov00, cov01, cov01, cov11);
      printf ("# sumsq = %g\n", sumsq);
//#endif
	}

}

