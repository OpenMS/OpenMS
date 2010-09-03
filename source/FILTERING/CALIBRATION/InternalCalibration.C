// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
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
    defaults_.setValue("rt_tolerance",10,"Allowed tolerance between peak and reference rt.");
		//		defaults_.setValue("hires:percentage",30,"Percentage of spectra a signal has to appear in to be considered as background signal.");
		defaultsToParam_();
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
			std::vector<DoubleReal> rel_errors(observed_masses.size(),0.);
			// determine rel error in ppm for the two reference masses
			for(Size ref_peak=0; ref_peak < observed_masses.size();++ref_peak)
			{
        rel_errors[ref_peak] = (theoretical_masses[ref_peak]-observed_masses[ref_peak])/theoretical_masses[ref_peak] * 1e6;

				out << observed_masses[ref_peak] << "\t"<< rel_errors[ref_peak] << "\n";
        std::cout << observed_masses[ref_peak] <<" "<<theoretical_masses[ref_peak]<<std::endl;
				//				std::cout << observed_masses[ref_peak]<<"\t"<<rel_errors[ref_peak]<<std::endl;
			}
#endif

			DoubleReal cov00, cov01, cov11, sumsq, slope,intercept;
			// TODO: what exactly is stride?? used 1 here as in the gsl-example :)
			gsl_fit_linear (&(observed_masses[0]), 1, &(theoretical_masses[0]), 1, observed_masses.size(), &intercept,&slope,&cov00,&cov01,&cov11,&sumsq);
			trafo_.setName("linear");
			trafo_.setParam("slope",slope);
			trafo_.setParam("intercept",intercept);

#ifdef DEBUG_CALIBRATION
			// 			std::cout <<"\n\n---------------------------------\n\n"<< "after calibration "<<std::endl;
			for(Size i = 0; i < observed_masses.size();++i)
				{
					DoubleReal new_mass = observed_masses[i];
					trafo_.apply(new_mass);

					DoubleReal rel_error = (theoretical_masses[i]-(new_mass))/theoretical_masses[i] * 1e6;
					std::cout << observed_masses[i]<<"\t"<<rel_error<<std::endl;
				}
#endif								
			
#ifdef DEBUG_CALIBRATION
  	  printf ("# best fit: Y = %g + %g X\n", intercept, slope);
      printf ("# covariance matrix:\n");
      printf ("# [ %g, %g\n#   %g, %g]\n", 
               cov00, cov01, cov01, cov11);
      printf ("# sumsq = %g\n", sumsq);
#endif
	}



	void InternalCalibration::calibrateMapGlobally(const FeatureMap<>& feature_map, FeatureMap<>& calibrated_feature_map,
																								 String trafo_file_name)
	{
		// check if the ids 
		checkReferenceIds_(feature_map);
		// first collect theoretical and observed m/z values
		std::vector<DoubleReal> observed_masses;
		std::vector<DoubleReal> theoretical_masses;
		for(Size f = 0; f < feature_map.size();++f)
			{
				// if more than one peptide id exists for this feature we can't use it as reference
				if(feature_map[f].getPeptideIdentifications().size() > 1) continue;
				if(!feature_map[f].getPeptideIdentifications().empty())
					{
						Int charge = feature_map[f].getPeptideIdentifications()[0].getHits()[0].getCharge();
						DoubleReal theo_mass = feature_map[f].getPeptideIdentifications()[0].getHits()[0].getSequence().getMonoWeight(Residue::Full,charge)/(DoubleReal)charge;
						theoretical_masses.push_back(theo_mass);
						observed_masses.push_back(feature_map[f].getMZ());
#ifdef DEBUG_CALIBRATION
						std::cout << feature_map[f].getRT() <<" " <<feature_map[f].getMZ() <<" "<<theo_mass<<std::endl;
						std::cout << feature_map[f].getPeptideIdentifications()[0].getHits().size()<<std::endl;
						std::cout << feature_map[f].getPeptideIdentifications()[0].getHits()[0].getSequence()<<std::endl;
						std::cout << feature_map[f].getPeptideIdentifications()[0].getHits()[0].getCharge()<<std::endl;
#endif
					}
			}
		// then make the linear regression
		makeLinearRegression_(observed_masses,theoretical_masses);
		// apply transformation
    applyTransformation_(feature_map,calibrated_feature_map);
		if(trafo_file_name != "")
			{
				TransformationXMLFile().store(trafo_file_name,trafo_);
			}
	}


	void InternalCalibration::calibrateMapGlobally(const FeatureMap<>& feature_map, FeatureMap<>& calibrated_feature_map,std::vector<PeptideIdentification>& ref_ids,String trafo_file_name)
	{
    checkReferenceIds_(ref_ids);
    
    calibrated_feature_map = feature_map;
    // clear the ids
    for(Size f = 0;f < calibrated_feature_map.size();++f)
    {
			calibrated_feature_map[f].getPeptideIdentifications().clear();
		}

    // map the reference ids onto the features
    IDMapper mapper;
    Param param;
    param.setValue("rt_tolerance",(DoubleReal)param_.getValue("rt_tolerance"));
    param.setValue("mz_tolerance",param_.getValue("mz_tolerance"));
    param.setValue("mz_measure",param_.getValue("mz_tolerance_unit"));
    mapper.setParameters(param);
    std::vector<ProteinIdentification> vec;
    mapper.annotate(calibrated_feature_map,ref_ids,vec);
    
    // calibrate
    calibrateMapGlobally(calibrated_feature_map,calibrated_feature_map,trafo_file_name);

		// copy the old ids
		calibrated_feature_map.setUnassignedPeptideIdentifications(feature_map.getUnassignedPeptideIdentifications());
    for(Size f = 0;f < feature_map.size();++f)
    {
			calibrated_feature_map[f].getPeptideIdentifications().clear();
			if(!feature_map[f].getPeptideIdentifications().empty())
				{
					calibrated_feature_map[f].setPeptideIdentifications(feature_map[f].getPeptideIdentifications());
				}
    }
	}

  void InternalCalibration::applyTransformation_(const FeatureMap<>& feature_map,FeatureMap<>& calibrated_feature_map)
  {
    calibrated_feature_map = feature_map;
		for(Size f = 0; f < feature_map.size();++f)
			{
				DoubleReal mz = feature_map[f].getMZ();
				trafo_.apply(mz);
				calibrated_feature_map[f].setMZ(mz);

				// apply transformation to convex hulls and subordinates
				for(Size s = 0; s < calibrated_feature_map[f].getSubordinates().size();++s)
					{
						// subordinates
						DoubleReal mz = calibrated_feature_map[f].getSubordinates()[s].getMZ();
						trafo_.apply(mz);
						calibrated_feature_map[f].getSubordinates()[s].setMZ(mz);
					}
				for(Size s = 0; s < calibrated_feature_map[f].getConvexHulls().size();++s)
					{
						// convex hulls
						std::vector<DPosition<2> > point_vec = calibrated_feature_map[f].getConvexHulls()[s].getHullPoints();
						calibrated_feature_map[f].getConvexHulls()[s].clear();
						for(Size p = 0; p < point_vec.size(); ++p)
						{
							DoubleReal mz = point_vec[p][1];
							trafo_.apply(mz);
							point_vec[p][1] = mz;
						}
						calibrated_feature_map[f].getConvexHulls()[s].setHullPoints(point_vec);						
					}
			}
  }

	void InternalCalibration::checkReferenceIds_(const FeatureMap<>& feature_map)
	{
		Size num_ids = 0;
		for(Size f = 0; f < feature_map.size();++f)
			{
				if(!feature_map[f].getPeptideIdentifications().empty() && feature_map[f].getPeptideIdentifications()[0].getHits().size() > 1)
					{
						throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "InternalCalibration: Your feature map contains PeptideIdentifications with more than one hit, use the IDFilter to select only the best hits before you map the ids to the feature map.");
					}
				else if(!feature_map[f].getPeptideIdentifications().empty())  ++num_ids;
			}
		if(num_ids < 2)
			{
				throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "InternalCalibration: Your feature map contains less than two PeptideIdentifications, can't perform a linear regression on your data.");
			}
	}

}

