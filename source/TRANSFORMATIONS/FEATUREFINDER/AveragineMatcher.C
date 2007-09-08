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
// $Maintainer: Ole Schulz-Trieglaff $	
// --------------------------------------------------------------------------


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/AveragineMatcher.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LogNormalModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseQuality.h>

#include <iostream>
#include <fstream>
#include <numeric>
#include <math.h>

#include <vector>

using namespace std;

namespace OpenMS
{
	
	namespace Internal
	{
		// Helper struct for AveragineMatcher
		struct ExpFitPolyData
		{
			size_t n;
			string profile;
		};
	}

	using namespace Internal;
	using namespace std;

	//positions and signal values
	vector<double> positionsDC2_;
	vector<double> signalDC2_;

	AveragineMatcher::AveragineMatcher()
		: BaseModelFitter(),
		quality_(0),
	 	model2D_(),
		mz_stat_(),
		rt_stat_(),
		stdev_mz_(0), 
		stdev_rt1_(0), 
		stdev_rt2_(0),
		min_(), 
		max_(),
		counter_(0),
		iso_stdev_first_(0),
		iso_stdev_last_(0),
		iso_stdev_stepsize_(0),
		first_mz_model_(0),
		last_mz_model_(0)
	{
		setName(getProductName());
				
		defaults_.setValue("tolerance_stdev_bounding_box",3.0f,"Bounding box has range [minimim of data, maximum of data] enlarged by tolerance_stdev_bounding_box times the standard deviation of the data");
		defaults_.setValue("intensity_cutoff_factor",0.05f,"Cutoff peaks with a predicted intensity below intensity_cutoff_factor times the maximal intensity of the model");
		defaults_.setValue("feature_intensity_sum",1,"Determines what is reported as feature intensity.\n1: the sum of peak intensities;\n0: the maximum intensity of all peaks");
		
		defaults_.setValue("min_num_peaks:final",5,"Minimum number of peaks left after cutoff. If smaller, feature will be discarded.");
		defaults_.setValue("min_num_peaks:extended",10,"Minimum number of peaks after extension. If smaller, feature will be discarded.");
		defaults_.setSectionDescription("min_num_peaks","Required number of peaks for a feature.");
		
		defaults_.setValue("rt:interpolation_step",0.2f,"Step size in seconds used to interpolate model for RT.");
		defaults_.setValue("rt:max_iteration",500,"Maximum number of iterations for RT fitting.");
		defaults_.setValue("rt:deltaAbsError",0.0001,"Absolute error used by the Levenberg-Marquardt algorithms.");
		defaults_.setValue("rt:deltaRelError",0.0001,"Relative error used by the Levenberg-Marquardt algorithms.");
		defaults_.setValue("rt:profile","EMG","Type of RT model. Currently only 'EMG' is supported.");
		defaults_.setSectionDescription("rt","Model settings in RT dimension.");
		
		defaults_.setValue("mz:interpolation_step",0.03f,"Interpolation step size for m/z.");
		defaults_.setValue("mz:model_type:first",1,"Numeric id of first m/z model fitted (usually indicating the charge state).");
		defaults_.setValue("mz:model_type:last",4,"Numeric id of last m/z model fitted (usually indicating the charge state).");
		defaults_.setSectionDescription("mz","Model settings in m/z dimension.");
		
		defaults_.setValue("quality:type","Correlation","Type of the quality measure used to assess the fit of model vs data ('Correlation','EuclidianDistance','RankCorrelation').");
		defaults_.setValue("quality:minimum",0.65f,"Minimum quality of fit, features below this threshold are discarded.");
		defaults_.setSectionDescription("quality","Fitting quality settings.");
		
		defaults_.setValue("isotope_model:stdev:first",0.04f,"First standard deviation to be considered for isotope model.");
		defaults_.setValue("isotope_model:stdev:last",0.12f,"Last standard deviation to be considered for isotope model.");
		defaults_.setValue("isotope_model:stdev:step",0.04f,"Step size for standard deviations considered for isotope model.");
		defaults_.setSectionDescription("isotope_model:stdev","Instrument resolution settings for m/z dimension.");
		
		defaults_.setValue("isotope_model:averagines:C",0.0443f,"Number of C atoms per Dalton of the mass.");
		defaults_.setValue("isotope_model:averagines:H",0.007f,"Number of H atoms per Dalton of the mass.");
		defaults_.setValue("isotope_model:averagines:N",0.0012f,"Number of N atoms per Dalton of the mass.");
		defaults_.setValue("isotope_model:averagines:O",0.013f,"Number of O atoms per Dalton of the mass.");
		defaults_.setValue("isotope_model:averagines:S",0.00037f,"Number of S atoms per Dalton of the mass.");
		defaults_.setSectionDescription("isotope_model:averagines","Averagines are used to approximate the number of atoms (C,H,N,O,S) which a peptide of a given mass contains.");
		
		defaults_.setValue("isotope_model:isotope:trim_right_cutoff",0.001f,"Cutoff for averagine distribution, trailing isotopes below this relative intensity are not considered.");
		defaults_.setValue("isotope_model:isotope:maximum",100,"Maximum number of isotopes being used for the IsotopeModel.");
		defaults_.setValue("isotope_model:isotope:distance",1.000495f,"Distance between consecutive isotopic peaks.");
		defaults_.setSectionDescription("isotope_model","Settings of the isotope model (m/z).");
		
		
		defaultsToParam_();
	}

	AveragineMatcher::~AveragineMatcher()
	{
		delete quality_;
	}

  AveragineMatcher::AveragineMatcher(const AveragineMatcher& rhs)
    : BaseModelFitter(rhs),
    	quality_(0)
  {
    updateMembers_();
  }
  
  AveragineMatcher& AveragineMatcher::operator= (const AveragineMatcher& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseModelFitter::operator=(rhs);
    
    updateMembers_();
    
    return *this;
  }

	void AveragineMatcher::updateMembers_()
	{
		if (quality_) delete quality_;
		quality_ = Factory<BaseQuality>::create(param_.getValue("quality:type"));
		
		max_iteration_ = param_.getValue("rt:max_iteration");
		eps_abs_        = param_.getValue("rt:deltaAbsError");
		eps_rel_         = param_.getValue("rt:deltaRelError");
		profile_           = (string)param_.getValue("rt:profile");

		interpolation_step_mz_ = param_.getValue("mz:interpolation_step");
		interpolation_step_rt_ = param_.getValue("rt:interpolation_step");
		
		iso_stdev_first_        = param_.getValue("isotope_model:stdev:first");
		iso_stdev_last_        = param_.getValue("isotope_model:stdev:last");
		iso_stdev_stepsize_ = param_.getValue("isotope_model:stdev:step");
		
		first_mz_model_ = (Int) param_.getValue("mz:model_type:first");
		last_mz_model_ = (Int) param_.getValue("mz:model_type:last");
	}
	
	AveragineMatcher::QualityType AveragineMatcher::fit_loop_(const ChargedIndexSet& set, Int& first_mz, Int& last_mz, CoordinateType& sampling_size_mz , ProductModel<2>*& final) 
	{
		QualityType quality = 0.0;
		QualityType max_quality = -std::numeric_limits<QualityType>::max();
		
		for ( float stdev = iso_stdev_first_; stdev <= iso_stdev_last_; stdev += iso_stdev_stepsize_)
		{
				for (Int mz_fit_type = first_mz; mz_fit_type <= last_mz; ++mz_fit_type)
				{
					quality = fit_(set, static_cast<MzFitting>(mz_fit_type), LMAGAUSS, stdev, (UInt) sampling_size_mz );
				
					if (quality > max_quality)
					{
						max_quality = quality;
						final = new ProductModel<2>(model2D_);	// store model
					}
				}
		}		
		return max_quality;
	}
	
	
	Feature AveragineMatcher::fit(const ChargedIndexSet& set) throw (UnableToFit)
	{
		// not enough peaks to fit
		if (set.size() < (UInt)(param_.getValue("min_num_peaks:extended")))
		{
			for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it) 
			{
				traits_->getPeakFlag(*it) =UNUSED;
			}
			
			String mess = String("Skipping feature, IndexSet size too small: ") + set.size();
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,"UnableToFit-IndexSet", mess.c_str());
		}

		quality_->setTraits(traits_);
		mz_lin_int_.getData().clear();		// empty interpolation datastructrure
		
		//QualityType quality = 0.0;
		QualityType max_quality = -std::numeric_limits<QualityType>::max();
		UInt max_peak_scan = 0;

		// Calculate statistics
		mz_stat_.update( IntensityIterator(set.begin(),traits_),IntensityIterator(set.end(),traits_),MzIterator(set.begin(),traits_) );
		rt_stat_.update ( IntensityIterator(set.begin(),traits_),IntensityIterator(set.end(),traits_),RtIterator(set.begin(),traits_) );

		// Calculate bounding box
		CoordinateType min_mz_distance = std::numeric_limits<CoordinateType>::max();
		CoordinateType prev_mz = 0.0;
		IntensityType highest_intensity = 0.0;
		IndexSetIter it=set.begin();
		min_ = max_ = traits_->getPeakPos(*it);
		for ( ++it; it!=set.end(); ++it )
		{
			CoordinateType tmp = traits_->getPeakMz(*it);
			if (min_[MZ] > tmp) min_[MZ] = tmp;
			if (max_[MZ] < tmp) max_[MZ] = tmp;
			
			if (it != set.begin())
			{
				if (fabs(tmp - prev_mz) < min_mz_distance && fabs(tmp - prev_mz)  > 0.0)
				{
					min_mz_distance = fabs(tmp - prev_mz);			
				}		
			}			
			prev_mz = tmp;		
			
			if (traits_->getPeakIntensity(*it) > highest_intensity)
			{
				highest_intensity = traits_->getPeakIntensity(*it);
				max_peak_scan = it->first;
			}
						
			tmp = traits_->getPeakRt(*it);
			if (min_[RT] > tmp) min_[RT] = tmp;
			if (max_[RT] < tmp) max_[RT] = tmp;
		}
		
		CoordinateType sampling_size_mz = 120; //2.0 / min_mz_distance;
		mz_lin_int_.getData().resize((UInt) sampling_size_mz);	
		mz_lin_int_.setMapping( 0, min_[MZ] , sampling_size_mz, max_[MZ]);
		
		for (IndexSetIter it=set.begin(); it!=set.end(); ++it )
		{
			mz_lin_int_.addValue( traits_->getPeakMz(*it),traits_->getPeakIntensity(*it) );			
		}				

		double const tolerance_stdev_box = param_.getValue("tolerance_stdev_bounding_box");
		stdev_mz_ = sqrt ( mz_stat_.variance() ) * tolerance_stdev_box;
		min_[MZ] -= stdev_mz_;
		max_[MZ] += stdev_mz_;

		stdev_rt1_ = sqrt ( rt_stat_.variance1() ) * tolerance_stdev_box;
		stdev_rt2_ = sqrt ( rt_stat_.variance2() ) * tolerance_stdev_box;
		min_[RT] -= stdev_rt1_;
		max_[RT] += stdev_rt2_;

		// create a vector with RT-values & Intensity
		// compute the parameters (intial values) for the EMG & Gauss function and finally,
		// optimize the parameters with Levenberg-Marquardt algorithms
		if (profile_=="LmaGauss" || profile_=="EMG" || profile_=="LogNormal")
		{
			setData(set);
			if (symmetric_==false)
				optimize();
				
			if (gsl_status_!="success") 
			{
				cout << profile_ + " status: " + gsl_status_ << endl;
				//throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,"UnableToFit-BadQuality",String("Skipping feature, " + profile_ + " status: " + gsl_status_));
			}
		}

		// Test charge states and stdevs
		Int first_mz  = first_mz_model_;
		Int last_mz  = last_mz_model_;
		
		// Check charge estimate if charge is not specified by user
		if (set.charge_ != 0)
		{
			first_mz = set.charge_;
			last_mz = set.charge_;
		}
		//cout << "Checking charge state from " << first_mz << " to " << last_mz << endl;
	
		ProductModel<2>* final = 0;	// model  with best correlation		
		fit_loop_(set, first_mz,last_mz,sampling_size_mz,final);
	
		// in this case something went wrong during the modelfitting and we stop.
		if (! final)
		{
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__, "UnableToFit-BadQuality","Zero quality after fitting. Skipping this feature");
			delete final;
		}
		
		// Cutoff low intensities wrt. to averagine model
		IndexSet model_set;
		for (IndexSetIter it=set.begin(); it!=set.end(); ++it) 
		{
			if (mz_model_.getIntensity( traits_->getPeakMz(*it) ) > IntensityType(param_.getValue("intensity_cutoff_factor")) &&
			    rt_model_.getIntensity( traits_->getPeakRt(*it) ) > IntensityType(param_.getValue("intensity_cutoff_factor")) )
			{
					model_set.insert(*it);
 			}
			else		// free dismissed peak by setting the UNUSED flag
			{
					traits_->getPeakFlag(*it) = UNUSED;
			}
		}
		// Print number of selected peaks after cutoff
		cout << " Selected " << model_set.size() << " from " << set.size() << " peaks." << std::endl;

		// not enough peaks left for feature
		if (model_set.size() < (UInt)(param_.getValue("min_num_peaks:final")))
		{				
			delete final;
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,
																	"UnableToFit-FinalSet",
												          String("Skipping feature, IndexSet size after cutoff too small: ") + model_set.size() );
		}
 		
		// normalize data and compute mean position
		IntensityType mz_data_sum  = 0.0;
		for (UInt i=0;i<mz_lin_int_.getData().size();++i)
		{
			mz_data_sum += mz_lin_int_.getData()[i];
		}	
		
		// compute m/z data average
		CoordinateType mz_data_avg = 0.0;
		for (UInt i=0;i<mz_lin_int_.getData().size();++i)
		{
			mz_data_avg += ( mz_lin_int_.getData()[i] / mz_data_sum);	
		}	
		mz_data_avg /= mz_lin_int_.getData().size();
	
		QualityType qual_mz = compute_mz_corr_(mz_data_sum, mz_model_, mz_data_avg);
	
// 		cout << "Quality in m/z : " << qual_mz << endl;
		QualityType qual_rt =	quality_->evaluate(model_set, rt_model_, RT );
		if (isnan(qual_rt) ) qual_rt = -1.0;
// 		cout << "Quality in rt : " << qual_rt << endl;
		
		max_quality = (qual_mz + qual_rt) / 2.0;

		if (max_quality < 0.2)
		{
// 			cout << "max quality is too low: " << max_quality << endl;

			Int fmz = first_mz_model_;
			Int lmz = last_mz_model_;
			
// 			cout << "Checking charge states " << fmz << " to " << lmz << endl;
			fit_loop_(set, fmz,lmz,sampling_size_mz,final);	
			
		}
		
		// fit has too low quality or fit was not possible i.e. because of zero stdev
		if (max_quality < (QualityType)(param_.getValue("quality:minimum")))
		{
			delete final;
			String mess = String("Skipping feature, correlation too small: ") + max_quality;
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,"UnableToFit-Correlation", mess.c_str());
		}
	
		// Build Feature
		// The feature coordinate in rt dimension is given
		// by the centroid of the rt model whereas the coordinate
		// in mz dimension is equal to the monoisotopic peak.
		Feature f;
		f.setModelDescription( ModelDescription<2>(final) );
		f.setOverallQuality(max_quality);
		f.setRT(dynamic_cast<InterpolationModel*>(final->getModel(RT))->getCenter());
		
		CoordinateType mz_guess = mz_model_.getCenter();
		vector< IDX > max_scan;
		for (IndexSetIter it=set.begin(); it!=set.end(); ++it )
		{		
			if (it->first == max_peak_scan)
			{
				max_scan.push_back( *it );
			}
		}
		sort(max_scan.begin(),max_scan.end());
		
		vector<IDX>::iterator iter = lower_bound(max_scan.begin(),max_scan.end(),mz_guess, AveragineMatcher::IndexMzLess::IndexMzLess(traits_));
		
		if ( iter == max_scan.end() ) 
	  {
	  	--iter; // one step back
	  }
		
		IntensityType max_mono_candidate_it       = 0.0;
		CoordinateType max_mono_candidate_mz = traits_->getPeakMz(*iter);
		CoordinateType mz_dist = 0.0;
		
		if ( iter == max_scan.end() ) 
	  {
	  	--iter; // one step back
	  }
		
		// local search for highest peak
		if (iter != max_scan.begin())
		{
			// walk to the left
			for (vector<IDX>::iterator it = iter; it != max_scan.begin() && mz_dist <= 0.25; --it )
			{
// 				cout << "it: " << it->first << " " << it->second << endl;
			
				if (it->first >= traits_->getData().size())
				{
					break;
				}
				
				if (it->second >= traits_->getData()[it->first].size())
				{
					break;
				}				
			
				if (traits_->getPeakIntensity(*it) > max_mono_candidate_it)
				{
					max_mono_candidate_it	   = traits_->getPeakIntensity(*it);
					max_mono_candidate_mz = traits_->getPeakMz(*it);				
				} 		
							
				mz_dist = traits_->getPeakMz(*iter) - traits_->getPeakMz(*it);			
			}
		}
		
		if ( iter != max_scan.end() ) 
		{		
			// walk to the right
			mz_dist = 0.0;
			for (vector<IDX>::iterator it = iter; it != max_scan.end() && mz_dist <= 0.25; ++it )
			{
				if (it->first >= traits_->getData().size())
				{
					break;
				}
				
				if (it->second >= traits_->getData()[it->first].size())
				{
					break;
				}				
				
				if (traits_->getPeakIntensity(*it) > max_mono_candidate_it)
				{
					max_mono_candidate_it	   = traits_->getPeakIntensity(*it);
					max_mono_candidate_mz = traits_->getPeakMz(*it);				
				} 		
				mz_dist = traits_->getPeakMz(*it) - traits_->getPeakMz(*iter);			
			}
		}		
		else
		{
			--iter; // one step back, no need to check left side
	  }
		
// 		cout << "Setting to mz " << 		max_mono_candidate_mz << endl;
// 		cout << "Old estimate " << mz_guess << endl;
		
		f.setMZ( max_mono_candidate_mz /*dynamic_cast<InterpolationModel*>(final->getModel(MZ))->getCenter()*/);
		if (final->getModel(MZ)->getName() == "IsotopeModel")
		{
			f.setCharge(dynamic_cast<IsotopeModel*>(final->getModel(MZ))->getCharge());
		}
		// if we used a simple Gaussian model to fit the feature, we can't say anything about
		// its charge state. The value 0 indicates that charge state is undetermined.
		else 
		{
			f.setCharge(0);		
		}

		Int const intensity_choice = param_.getValue("feature_intensity_sum");
		IntensityType feature_intensity = 0.0;
		
		if (intensity_choice == 1)
		{
			// intensity of the feature is the sum of all included data points
			for (IndexSetIter it=model_set.begin(); it!=model_set.end(); ++it) 
			{
				feature_intensity += traits_->getPeakIntensity(*it);
			}
		}
		else
		{
			// feature intensity is the maximum intensity of all peaks
			for (IndexSetIter it=model_set.begin(); it!=model_set.end(); ++it) 
			{
				if (traits_->getPeakIntensity(*it) > feature_intensity)
					feature_intensity = traits_->getPeakIntensity(*it);
			}	
		} 
		
		f.setIntensity(feature_intensity);
		traits_->addConvexHull(model_set, f);

		std::cout << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << " Feature " << counter_
							<< ": (" << f.getRT()
							<< "," << f.getMZ() << ") Qual.:"
							<< max_quality << "\n";


		f.setQuality(RT, qual_rt);
		f.setQuality(MZ, qual_mz);

		// save meta data in feature for TOPPView
		stringstream meta ;
		meta << "Feature #" << counter_ << ", +"	<< f.getCharge() << ", " << set.size() << "->" << model_set.size() 
				 << ", Corr: (" << max_quality << ","  << f.getQuality(RT) << "," << f.getQuality(MZ) << ")";
		f.setMetaValue(3,String(meta.str()));

		#ifdef DEBUG_FEATUREFINDER
		// write debug output
		CoordinateType rt = f.getRT();
		CoordinateType mz = f.getMZ();
		
		// write feature model 
		String fname = String("model")+ counter_ + "_" + rt + "_" + mz;
		ofstream file(fname.c_str()); 
		for (IndexSetIter it=model_set.begin(); it!=model_set.end(); ++it) 
		{
			DPosition<2> pos = traits_->getPeakPos(*it);
			if ( final->isContained(pos) )
			{
				file << pos[RT] << " " << pos[MZ] << " " << final->getIntensity( traits_->getPeakPos(*it)) << "\n";						
			}
		}
		file.close();
	
		// wrote peaks remaining after model fit
		fname = String("feature") + counter_ +  "_" + rt + "_" + mz;
		ofstream file2(fname.c_str()); 
		for (IndexSetIter it=model_set.begin(); it!=model_set.end(); ++it) 
		{
			DPosition<2> pos = traits_->getPeakPos(*it);
			if ( final->isContained(pos) )
			{
				file2 << pos[RT] << " " << pos[MZ] << " " << traits_->getPeakIntensity(*it) << "\n";						
			}
		}
		file2.close();
		#endif

		++counter_;
		
		delete final;
		
 		return f;
	}

	AveragineMatcher::QualityType AveragineMatcher::fit_(const ChargedIndexSet& set, MzFitting mz_fit, RtFitting /*rt_fit*/,
																	 					  																	  Coordinate isotope_stdev, UInt sampling_size)
	{
			// Build Models
			rt_model_.setInterpolationStep(interpolation_step_rt_);

			Param tmp;
			tmp.setValue("bounding_box:min",min_[RT] );
			tmp.setValue("bounding_box:max",max_[RT] );
			tmp.setValue("statistics:variance",rt_stat_.variance() );
			tmp.setValue("statistics:mean",rt_stat_.mean() );			
			tmp.setValue("emg:height",height_);
			tmp.setValue("emg:width",width_);
			tmp.setValue("emg:symmetry",symmetry_);
			tmp.setValue("emg:retention",retention_);
	
			rt_model_.setParameters(tmp);
			double res = fit_mz_(set, sampling_size, mz_fit,isotope_stdev);
			
			model2D_.setModel(MZ, &mz_model_).setModel(RT, &rt_model_);
			
			return res;
	}
	
	void AveragineMatcher::dump_all_(ChargedIndexSet set, UInt sampling_size)
	{
			// dumping linear interpolation DS
			String filename = "dump_out_" + String(counter_) + String("_") + String(sampling_size);
			std::ofstream outfile(filename.c_str());
			
			for (UInt i=0;i<mz_lin_int_.getData().size();++i)
			{
			outfile << mz_lin_int_.index2key(i) << " ";
			outfile << (mz_lin_int_.getData()[i] ) << endl;
			}
			outfile.close();
		
			// feature region
			String filename2 = "dump_region_out_" + String(counter_) + String("_") + String(sampling_size);
			outfile.open(filename2.c_str());
			for (ChargedIndexSet::const_iterator citer = set.begin();
						citer != set.end();
						++citer)
			{
				outfile << traits_->getPeakRt(*citer) << " " << traits_->getPeakMz(*citer) << " " << traits_->getPeakIntensity(*citer) << std::endl;
			}
			outfile.close();
	
	}
	
	AveragineMatcher::QualityType AveragineMatcher::compute_mz_corr_(IntensityType& mz_data_sum, 
																																																				IsotopeModel& iso_model,
																																																				CoordinateType& mz_data_avg)
	{	
			// normalize m/z model
			IntensityType mz_model_sum = 0.0;
			DPeakArray< DPeak<1> > samples;
			iso_model.getSamples(samples);
			for(UInt i=0; i<samples.size();++i)
			{
				mz_model_sum += samples[i].getIntensity();	
			}
		
			// compute m/z model average (for the given set of data points)
			CoordinateType mz_model_avg = 0.0;		
			for (UInt i=0;i<mz_lin_int_.getData().size();++i)
			{
				mz_model_avg += ( iso_model.getIntensity(mz_lin_int_.index2key(i) ) / mz_model_sum);
			} 		
			mz_model_avg /= mz_lin_int_.getData().size();
			
			// compute Pearson correlation
			IntensityType cross_product_sum = 0;
			IntensityType data_square_sum   = 0;
			IntensityType model_square_sum  = 0; 
			
			for (UInt i=0;i<mz_lin_int_.getData().size();++i)
			{
					IntensityType m = (iso_model.getIntensity(mz_lin_int_.index2key(i) ) / mz_model_sum);
					IntensityType d = (mz_lin_int_.getData()[i] / mz_data_sum);
			
					cross_product_sum += ( m - mz_model_avg) * ( d - mz_data_avg);
					data_square_sum    += ( d - mz_data_avg)  * ( d - mz_data_avg);
					model_square_sum  += ( m - mz_model_avg)  * ( m - mz_model_avg);					
			}
		
			if (data_square_sum == 0 || model_square_sum == 0)
			{
				return -2.0;			
			}
			
			QualityType corr_mz = cross_product_sum / sqrt(data_square_sum * model_square_sum);
	
			return corr_mz;	
	}
	
	AveragineMatcher::QualityType AveragineMatcher::fit_mz_(ChargedIndexSet /*set*/, UInt /*samplingsize*/, MzFitting charge,Coordinate isotope_stdev)
	{			
		// new model
		IsotopeModel iso_model;
		Param iso_param = param_.copy("isotope_model:",true);
		iso_param.remove("stdev");
		iso_model.setParameters(iso_param);
		iso_model.setInterpolationStep(interpolation_step_mz_);
			
		QualityType max_corr         =  -std::numeric_limits<QualityType>::max();
		CoordinateType max_center =  -std::numeric_limits<QualityType>::max();
				
		// normalize data and compute mean position
		IntensityType mz_data_sum  = 0.0;
		CoordinateType mz_mean_pos = 0.0;
		for (UInt i=0;i<mz_lin_int_.getData().size();++i)
		{
			mz_mean_pos += (mz_lin_int_.index2key(i) * mz_lin_int_.getData()[i]);
			mz_data_sum += mz_lin_int_.getData()[i];
		}	
		
		if (mz_data_sum == 0.0)
		{
// 			dump_all_(set,(UInt) samplingsize);
			return -1.0;		
		}
		
		mz_mean_pos /= mz_data_sum;
			
// 		cout << "mz_mean_pos = " << mz_mean_pos << endl;
		
		// compute m/z data average
		CoordinateType mz_data_avg = 0.0;
		for (UInt i=0;i<mz_lin_int_.getData().size();++i)
		{
			mz_data_avg += ( mz_lin_int_.getData()[i] / mz_data_sum);	
		}	
		mz_data_avg /= mz_lin_int_.getData().size();
		
		for (CoordinateType pos = mz_mean_pos- 0.5;
		     pos <= mz_mean_pos+ 0.5;
				 pos += 0.1)
		{
		
			Param tmp;
			tmp.setValue("charge", static_cast<Int>(charge));
			tmp.setValue("isotope:stdev",isotope_stdev);
			tmp.setValue("statistics:mean", pos);
	
			iso_model.setParameters(tmp);
			iso_model.setSamples();
					
			// estimate goodness of m/z fit
			QualityType corr_mz = compute_mz_corr_(mz_data_sum, iso_model, mz_data_avg);			
			
// 			if (corr_mz == -2.0)
// 			{
// 					dump_all_(set,(UInt) samplingsize);
// 			}
			
			if (corr_mz > max_corr)
			{
				max_corr   = corr_mz;
				max_center = pos;
			}
				
		}		
		mz_model_ = IsotopeModel();
		Param p = param_.copy("isotope_model:",true);
		p.remove("stdev");
		mz_model_.setParameters(p);
		mz_model_.setInterpolationStep(interpolation_step_mz_);
			
		Param tmp;
		tmp.setValue("charge", static_cast<Int>(charge));
		tmp.setValue("isotope:stdev",isotope_stdev);
		tmp.setValue("statistics:mean", max_center);
		mz_model_.setParameters(tmp);
			
// 		cout << "Setting averagine model to center (finally) " << mz_model_.getCenter() << endl;
// 		cout << "Best center: " << max_center << " with correlation " << max_corr << endl;
// 		
// 		cout << "in fitOffset_mz_: " << 	mz_model_ << endl;
		
		return max_corr;
	}


	AveragineMatcher::QualityType AveragineMatcher::fitOffset_(	InterpolationModel* model,
																																													const IndexSet& set, const double stdev1,  const double stdev2,
																																													const Coordinate offset_step)
	{
		const Coordinate offset_min = model->getInterpolation().supportMin() - stdev1;
		const Coordinate offset_max = model->getInterpolation().supportMin() + stdev2;

		Coordinate offset;
		double correlation;

		//test model with default offset
		Coordinate max_offset = model->getInterpolation().getOffset();
		double max_correlation = quality_->evaluate(set, model2D_);

		//test different offsets
		for ( offset = offset_min; offset <= offset_max; offset += offset_step )
		{
			model->setOffset(offset);
			correlation = quality_->evaluate(set, model2D_);
			if ( correlation > max_correlation )
			{
				max_correlation = correlation;
				max_offset = offset;
			}
		}
		model->setOffset(max_offset);
		return max_correlation;
	}

	//create a vector with RT-values & Intensities and compute the parameters (intial values) for the EMG & Gauss function
	void AveragineMatcher::setData(const IndexSet& set)
	{
		// sum over all intensities
		double sum = 0.0;

		// iterate over all points of the signal
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); it++)
		{
			// store the current rt-position and signal
			float position = traits_->getPeakRt(*it);
			float signal = traits_->getPeakIntensity(*it);

			//float mz = traits_->getPeakMz(*it);
			sum+=signal;

			//orgFile << position << "  " << mz << " " << signal << "\n";

			// fill vectors with rt-postion and signal
			if (positionsDC2_.empty() || positionsDC2_.back()!=position) {
				positionsDC2_.push_back(position);
				signalDC2_.push_back(signal);
			}
			else {
				signal += signalDC2_.back();
				signalDC2_.pop_back();
				signalDC2_.push_back(signal);
			}
		}

		// calculate the median
		int median = 0;
		float count = 0.0;
		for (size_t current_point=0; current_point<positionsDC2_.size();current_point++)
		{
			count += signalDC2_[current_point];
			if (count <= sum/2)
				median = current_point;
		}

		double sumS = 0.0;
		for (size_t current_point=0; current_point<positionsDC2_.size();current_point++)
		{
			sumS += pow((positionsDC2_[current_point] - positionsDC2_[median]),2);
		}

		// calculate the stardard deviation
		standard_deviation_ = sqrt(sumS/(positionsDC2_.size() - 1));

		// set expected value
		expected_value_ = positionsDC2_[median];//rt_stat_.mean();

		// calculate the peak height
		height_ = signalDC2_[median];

		// calculate the width of the peak
		// rt-values with intensity zero are not allowed for calculation of the width
		width_ = abs(positionsDC2_[positionsDC2_.size()-1] - positionsDC2_[0]);

		// calculate retention time
		retention_ = positionsDC2_[median];

		// default is an asymmetric peak
		symmetric_ = false;

		// calculate the symmetry (fronted peak: s<1 , tailed peak: s>1)
		symmetry_ = abs(positionsDC2_.back() - positionsDC2_[median])/abs(positionsDC2_[median] - positionsDC2_.front());

		// check the symmetry
		if (isinf(symmetry_) || isnan(symmetry_))
		{
			symmetric_ = true;
			symmetry_ = 10;
		}

		// optimize the symmetry
		if (profile_=="LogNormal") {
			// The computations can lead to an overflow error at very low values of symmetry (s~0).
			
			if (symmetry_<=0.8)
				symmetry_=0.8;
			
			if (symmetry_==1)
				symmetry_=1.1;

			if (symmetry_>=1.5)
				symmetry_=1.4;

			// it is better to proceed from narrow peaks
			width_ /= 2;
		}
		else
		{
			// The computations can lead to an overflow error at very low values of symmetry (s~0). For s~5 the parameter can be aproximized by the Levenberg-Marquardt argorithms. (the other parameters are much greater than one)
				
			if (symmetry_<1)
				symmetry_+=5;

			// it is better for the emg function to proceed from narrow peaks
			width_ = symmetry_;
		}

		/* set the parameter r of the log normal function;
		  r is the ratio between h and the height at which w and s are computed;
		  r = 2, see "Mathematical functions for representation of chromatographic peaks", V.B. Di Marco(2001)
		*/
		r_ = 2;
	}

	//Evaluation of the target function for nonlinear optimization.
	int residualDC2(const gsl_vector* x, void* params , gsl_vector* f)
	{
		size_t n = ((struct ExpFitPolyData*)params)->n;
		String profile = ((struct ExpFitPolyData*)params)->profile;

		/// normal distribution (s = standard deviation, m = expected value)
		if (profile=="LmaGauss")
		{
			double normal_s = gsl_vector_get(x,0);
			double normal_m = gsl_vector_get(x,1);
			double normal_scale = gsl_vector_get(x,2);

			double Yi = 0.0;

			for (size_t i = 0; i < n; i++)
			{
				double t = positionsDC2_[i];

				Yi=(1/(sqrt(2*M_PI)*normal_s))*exp(-((t-normal_m)*(t-normal_m))/(2*normal_s*normal_s))*normal_scale;

				gsl_vector_set(f, i, (Yi - signalDC2_[i]));
			}
		}
		else
		{
			/// Simplified EMG
			if (profile=="EMG")
			{
				double h = gsl_vector_get(x,0);
				double w = gsl_vector_get(x,1);
				double s = gsl_vector_get(x,2);
				double z = gsl_vector_get(x,3);

				double Yi = 0.0;

				// iterate over all points of the signal
				for (size_t i = 0; i < n; i++)
				{
					double t = positionsDC2_[i];

					// Simplified EMG
					Yi=(h*w/s)*sqrt(2*M_PI)*exp(((w*w)/(2*s*s))-((t-z)/s))/(1+exp((-2.4055/sqrt(2))*(((t-z)/w)-w/s)));

					gsl_vector_set(f, i, (Yi - signalDC2_[i]));
				}
			}
			/// log normal
			else
			{
				double h = gsl_vector_get(x,0);
				double w = gsl_vector_get(x,1);
				double s = gsl_vector_get(x,2);
				double z = gsl_vector_get(x,3);
				double r = 2;//gsl_vector_get(x,4);

				double Yi = 0.0;

				for (size_t i = 0; i < n; i++)
				{
					double t = positionsDC2_[i];

					Yi = h*exp(-log(r)/(log(s)*log(s))*pow(log((t-z)*(s*s-1)/(w*s)+1),2));

					gsl_vector_set(f, i, (Yi - signalDC2_[i]));
				}
			}
		}

		return GSL_SUCCESS;
	}

	/** Compute the Jacobian of the residual, where each row of the matrix corresponds to a
	 *  point in the data. */
	int jacobianDC2(const gsl_vector* x, void* params, gsl_matrix* J)
	{

		size_t n = ((struct ExpFitPolyData*)params)->n;
		String profile = ((struct ExpFitPolyData*)params)->profile;

		// normal distribution (s = standard deviation, m = expected value)
		if (profile=="LmaGauss")
		{
			double normal_s = gsl_vector_get(x,0);
			double normal_m = gsl_vector_get(x,1);
			double normal_scale = gsl_vector_get(x,2);

			double derivative_normal_s, derivative_normal_m, derivative_normal_scale = 0.0;

			for (size_t i = 0; i < n; i++)
			{
				double t = positionsDC2_[i];

				// f'(normal_s)
				derivative_normal_s = -((1/sqrt(2*M_PI))/(normal_s*normal_s))*exp(-((t-normal_m)*(t-normal_m))/(2*normal_s*normal_s))*normal_scale+((1/sqrt(2*M_PI))/(normal_s*normal_s*normal_s*normal_s))*((t-normal_m)*(t-normal_m))*exp(-((t-normal_m)*(t-normal_m))/(2*normal_s*normal_s))*normal_scale;

				// f'(normal_m)
				derivative_normal_m = ((1/sqrt(2*M_PI))/(normal_s*normal_s*normal_s))*(t-normal_m)*exp(-((t-normal_m)*(t-normal_m))/(2*normal_s*normal_s))*normal_scale;

				// f'(normal_scale)
				derivative_normal_scale = ((1/sqrt(2*M_PI))/(normal_s))*exp(-((t-normal_m)*(t-normal_m))/(2*normal_s*normal_s));

				// set the jacobian matrix of the normal distribution
				gsl_matrix_set(J, i, 0, derivative_normal_s);
				gsl_matrix_set(J, i, 1, derivative_normal_m);
				gsl_matrix_set(J, i, 2, derivative_normal_scale);
			}
		}
		else
		{
			//Simplified EMG (sEMG)
			if (profile=="EMG")
			{
				double h = gsl_vector_get(x,0);
				double w = gsl_vector_get(x,1);
				double s = gsl_vector_get(x,2);
				double z = gsl_vector_get(x,3);

				const double emg_const = 2.4055;
				const double sqrt_2pi = sqrt(2*M_PI);
				const double sqrt_2   = sqrt(2);

				double exp1, exp2, exp3 = 0.0;
				double derivative_height, derivative_width, derivative_symmetry, derivative_retention = 0.0;

				// iterate over all points of the signal
				for (size_t i = 0; i < n; i++)
				{
					double t = positionsDC2_[i];

					exp1 = exp(((w*w)/(2*s*s))-((t-z)/s));
					exp2 = (1+exp((-emg_const/sqrt_2)*(((t-z)/w)-w/s)));
					exp3 = exp((-emg_const/sqrt_2)*(((t-z)/w)-w/s));

					// f'(h) - sEMG
					derivative_height = w/s*sqrt_2pi*exp1/exp2;

					// f'(h) - sEMG
					derivative_width = h/s*sqrt_2pi*exp1/exp2 + (h*w*w)/(s*s*s)*sqrt_2pi*exp1/exp2 + (emg_const*h*w)/s*sqrt_2pi*exp1*(-(t-z)/(w*w)-1/s)*exp3/((exp2*exp2)*sqrt_2);

					// f'(s) - sEMG
					derivative_symmetry = - h*w/(s*s)*sqrt_2pi*exp1/exp2 +  h*w/s*sqrt_2pi*(-(w*w)/(s*s*s)+(t-z)/(s*s))*exp1/exp2 + (emg_const*h*w*w)/(s*s*s)*sqrt_2pi*exp1*exp3/((exp2*exp2)*sqrt_2);

					// f'(z) - sEMG
					derivative_retention = h*w/(s*s)*sqrt_2pi*exp1/exp2 - (emg_const*h)/s*sqrt_2pi*exp1*exp3/((exp2*exp2)*sqrt_2);

					// set the jacobian matrix
					gsl_matrix_set(J, i, 0, derivative_height);
					gsl_matrix_set(J, i, 1, derivative_width);
					gsl_matrix_set(J, i, 2, derivative_symmetry);
					gsl_matrix_set(J, i, 3, derivative_retention);
				}
			}
			// log normal function
			else
			{
				double h = gsl_vector_get(x,0);
				double w = gsl_vector_get(x,1);
				double s = gsl_vector_get(x,2);
				double z = gsl_vector_get(x,3);
				double r = 2;//gsl_vector_get(x,4);

				double derivative_height, derivative_width, derivative_symmetry, derivative_retention,  derivative_r = 0.0;

				// iterate over all points of the signal
				for (size_t i = 0; i < n; i++)
				{
					double t = positionsDC2_[i];

					double exp1  = exp(-log(r)/(log(s)*log(s))*pow(log((t-z)*(s*s-1)/(w*s)+1),2));
					double term1 = (((t-z)*(s*s-1))/(w*s))+1;
					double log_s = log(s);
					double log_term1 = log(term1);
					double log_r = log(r);

					derivative_height = exp1;

					derivative_width = 2*h*log_r/(log_s*log_s)*log_term1*(t-z)*(s*s-1)/(w*w)/s/term1*exp1;

					derivative_symmetry = h*(2*log_r/(log_s*log_s*log_s)*(log_term1*log_term1)/s-2*log_r/(log_s*log_s)*log_term1*(2*(t-z)/w-(t-z)*(s*s-1)/(w*s*s))/term1)*exp1;

					derivative_retention = 2*h*log_r/(log_s*log_s)*log_term1*(s*s-1)/(w*s)/term1*exp1;

					derivative_r = 	-h/r/(log_s*log_s)*(log_term1*log_term1)*exp1;

					// set the jacobian matrix
					gsl_matrix_set(J, i, 0, derivative_height);
					gsl_matrix_set(J, i, 1, derivative_width);
					gsl_matrix_set(J, i, 2, derivative_symmetry);
					gsl_matrix_set(J, i, 3, derivative_retention);
					//gsl_matrix_set(J, i, 4, derivative_r);
				}
			}
		}

		return GSL_SUCCESS;
	}

	// Driver function for the evaluation of function and jacobian.
	int evaluateDC2(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
	{
		residualDC2(x, params, f);
		jacobianDC2(x, params, J);

		return GSL_SUCCESS;
	}

	// perform a nonlinear optimization
	void AveragineMatcher::optimize()
	{
		const gsl_multifit_fdfsolver_type *T;
		gsl_multifit_fdfsolver *s;

		int status;
		size_t iter = 0;
		const size_t n = positionsDC2_.size();

		// number of parameter to be optimize
		unsigned int p = 0;
		if (profile_=="LmaGauss") p = 3;
		else if (profile_=="LogNormal") p = 4;//5;
		else p = 4;

		// gsl always excepts N>=p or default gsl error handler invoked, cause Jacobian be rectangular M x N with M>=N	
		if (n<p)
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,"UnableToFit-FinalSet","Skipping feature, gsl always expects N>=p");
 
		gsl_matrix *covar = gsl_matrix_alloc(p,p);
		gsl_multifit_function_fdf f;

		double x_init_normal[3] = { standard_deviation_, expected_value_, height_ };
		double x_init_emg[4] = { height_, width_, symmetry_, retention_ };
		//double x_init_lognormal[5] = { height_, width_, symmetry_, retention_, r_ };
		double x_init_lognormal[4] = { height_, width_, symmetry_, retention_ };

		gsl_vector_view x;

		if (profile_=="LmaGauss")
			x = gsl_vector_view_array(x_init_normal, p);
		else
		{
			if (profile_=="LogNormal")
				x = gsl_vector_view_array(x_init_lognormal, p);
			else
				x = gsl_vector_view_array(x_init_emg, p);
		}

		const gsl_rng_type * type;
		gsl_rng * r;
		gsl_rng_env_setup();
		type = gsl_rng_default;
		r = gsl_rng_alloc (type);

		struct ExpFitPolyData d = {n, profile_};
		f.f = &residualDC2;
		f.df = &jacobianDC2;
		f.fdf = &evaluateDC2;
		f.n = n;
		f.p = p;
		f.params = &d;

		T = gsl_multifit_fdfsolver_lmsder;
		s = gsl_multifit_fdfsolver_alloc(T,n,p);
		gsl_multifit_fdfsolver_set(s, &f, &x.vector);

#ifdef DEBUG_FEATUREFINDER
		if (profile_=="LmaGauss") {
			printf ("before loop iter: %4u x = % 15.8f % 15.8f % 15.8f |f(x)| = %g\n", iter,
							gsl_vector_get(s->x,0),
							gsl_vector_get(s->x,1),
							gsl_vector_get(s->x,2),
							gsl_blas_dnrm2(s->f));
		}
		else {
			if (profile_=="EMG") {
				printf ("before loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
								gsl_vector_get(s->x,0),
								gsl_vector_get(s->x,1),
								gsl_vector_get(s->x,2),
								gsl_vector_get(s->x,3),
								gsl_blas_dnrm2(s->f));
			}
			else {
				printf ("before loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
								gsl_vector_get(s->x,0),
								gsl_vector_get(s->x,1),
								gsl_vector_get(s->x,2),
								gsl_vector_get(s->x,3),
								//gsl_vector_get(s->x,4),
								gsl_blas_dnrm2(s->f));

			}
		}
#endif

		// this is the loop for fitting
		do
		{
			iter++;
			status = gsl_multifit_fdfsolver_iterate (s);

#ifdef DEBUG_FEATUREFINDER
			// This is useful for debugging
			if (profile_=="LmaGauss") 	{
				printf ("in loop iter: %4u x = % 15.8f % 15.8f % 15.8f |f(x)| = %g\n", iter,
								gsl_vector_get(s->x,0),
								gsl_vector_get(s->x,1),
								gsl_vector_get(s->x,2),
								gsl_blas_dnrm2(s->f));
			}
			else {
				if (profile_=="EMG") {
					printf ("in loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
									gsl_vector_get(s->x,0),
									gsl_vector_get(s->x,1),
									gsl_vector_get(s->x,2),
									gsl_vector_get(s->x,3),
									gsl_blas_dnrm2(s->f));
				}
				else {
					printf ("in loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
									gsl_vector_get(s->x,0),
									gsl_vector_get(s->x,1),
									gsl_vector_get(s->x,2),
									gsl_vector_get(s->x,3),
									//gsl_vector_get(s->x,4),
									gsl_blas_dnrm2(s->f));
				}
			}
#endif

			// fit is done
			if (status) break;
			status = gsl_multifit_test_delta(s->dx, s->x, eps_abs_, eps_rel_);
		}
		while (status == GSL_CONTINUE && iter < max_iteration_);

		// This function uses Jacobian matrix J to compute the covariance matrix of the best-fit parameters, covar. The parameter epsrel (0.0) is used to remove linear-dependent columns when J is rank deficient.
		gsl_multifit_covar(s->J, 0.0, covar);

#ifdef DEBUG_FEATUREFINDER
		gsl_matrix_fprintf(stdout, covar, "%g");
#endif

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

		gsl_status_ = gsl_strerror(status);

#ifdef DEBUG_FEATUREFINDER
		cout << profile_ << " status: " << gsl_status_ << endl;
#endif

		if (profile_=="LmaGauss")
		{
#ifdef DEBUG_FEATUREFINDER
			printf("deviation          = %.5f +/- %.5f\n", FIT(0), ERR(0));
			printf("expected_value	   = %.5f +/- %.5f\n", FIT(1), ERR(1));
			printf("scale_factor       = %.5f +/- %.5f\n", FIT(2), ERR(2));
#endif
			standard_deviation_ = FIT(0);
			expected_value_     = FIT(1);
			scale_factor_       = FIT(2);
		}
		else
		{
			if (profile_=="EMG")
			{
#ifdef DEBUG_FEATUREFINDER
				printf("h = %.5f +/- %.5f\n", FIT(0), ERR(0));
				printf("w = %.5f +/- %.5f\n", FIT(1), ERR(1));
				printf("s = %.5f +/- %.5f\n", FIT(2), ERR(2));
				printf("z = %.5f +/- %.5f\n", FIT(3), ERR(3));
#endif
				height_     = FIT(0);
				width_      = FIT(1);
				symmetry_   = FIT(2);
				retention_  = FIT(3);
			}
			else
			{
#ifdef DEBUG_FEATUREFINDER
				printf("h = %.5f +/- %.5f\n", FIT(0), ERR(0));
				printf("w = %.5f +/- %.5f\n", FIT(1), ERR(1));
				printf("s = %.5f +/- %.5f\n", FIT(2), ERR(2));
				printf("z = %.5f +/- %.5f\n", FIT(3), ERR(3));
			//	printf("r = %.5f +/- %.5f\n", FIT(4), ERR(4));
#endif
				height_     = FIT(0);
				width_      = FIT(1);
				symmetry_   = FIT(2);
				retention_  = FIT(3);
				r_ 	    = r_;//FIT(4);
			}
		}

#ifdef DEBUG_FEATUREFINDER
		{
			// chi-squared value
			double chi = gsl_blas_dnrm2(s->f);
			printf("chisq/dof = %g\n", pow(chi, 2.0)/ (n-p));
		}
#endif

		// function free all memory associated with the solver s
		gsl_multifit_fdfsolver_free(s);

#ifdef DEBUG_FEATUREFINDER
		for (size_t current_point=0; current_point<positionsDC2_.size();current_point++)
			std::cout << positionsDC2_[current_point] << " " << signalDC2_[current_point] << std::endl;

		cout << "" << endl;
		cout << "*** parameter for optimization ***" << endl;
		cout << "       height:  " << height_ << endl;
		cout << "        width:  " << width_ << endl;
		cout << "     symmetry:  " << symmetry_ << endl;
		cout << "    retention:  " << retention_ << endl;
		cout << "std.deviation:  " << standard_deviation_ << endl;
		cout << "max_iteration:  " << max_iteration_ << endl;
		cout << "      eps_abs:  " << eps_abs_ << endl;
		cout << "      eps_rel:  " << eps_rel_ << endl;
		cout << "      profile:  " << profile_ << endl;
		cout << "" << endl;
#endif
		positionsDC2_.clear();
		signalDC2_.clear();
	}

 	AveragineMatcher::CoordinateType AveragineMatcher::getSymmetry() const 
	{ 
		return symmetry_; 
	}
	
	AveragineMatcher::CoordinateType AveragineMatcher::getHeight() const 
	{ 
		return height_; 
	}
	
	AveragineMatcher::CoordinateType AveragineMatcher::getWidth() const 
	{ 
		return width_; 
	}
	
	AveragineMatcher::CoordinateType AveragineMatcher::getRT() const 
	{ 
		return retention_; 
	}

	AveragineMatcher::CoordinateType AveragineMatcher::getStandardDeviation() const 
	{ 
		return standard_deviation_; 
	}

	AveragineMatcher::CoordinateType AveragineMatcher::getExpectedValue() const 
	{ 
		return expected_value_; 
	}

	AveragineMatcher::CoordinateType AveragineMatcher::getScaleFactor() const 
	{ 
		return scale_factor_; 
	}

	std::string AveragineMatcher::getGSLStatus() const 
	{ 
		return gsl_status_; 
	}

}
