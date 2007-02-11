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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseQuality.h>

#include <iostream>
#include <fstream>
#include <numeric>


using namespace std;

namespace OpenMS
{
	using namespace Internal;

	SimpleModelFitter::SimpleModelFitter()
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
		counter_(0)
	{
		setName(getProductName());
		
		defaults_.setValue("tolerance_stdev_bounding_box",3.0f);
		defaults_.setValue("feature_intensity_sum",1);
		defaults_.setValue("min_num_peaks:final",5);
		defaults_.setValue("min_num_peaks:extended",10);
		defaults_.setValue("intensity_cutoff_factor",0.05f);
		defaults_.setValue("rt:interpolation_step",0.2f);
		defaults_.setValue("mz:interpolation_step",0.03f);
		defaults_.setValue("mz:model_type:first",0);
		defaults_.setValue("mz:model_type:last",4);
		defaults_.setValue("quality:type","Correlation");
		defaults_.setValue("quality:minimum",0.65f);
		defaults_.setValue("isotope_model:stdev:first",0.04f);
		defaults_.setValue("isotope_model:stdev:last",0.12f);
		defaults_.setValue("isotope_model:stdev:step",0.04f);
		defaults_.setValue("isotope_model:averagines:C",0.0443f);
		defaults_.setValue("isotope_model:averagines:H",0.0f);
		defaults_.setValue("isotope_model:averagines:N",0.0037f);
		defaults_.setValue("isotope_model:averagines:O",0.022f);
		defaults_.setValue("isotope_model:averagines:S",0.0f);
		defaults_.setValue("isotope_model:isotope:trim_right_cutoff",0.001f);
		defaults_.setValue("isotope_model:isotope:maximum",1000000);
		defaults_.setValue("isotope_model:isotope:distance",1.000495f);
		
		defaultsToParam_();
	}

	SimpleModelFitter::~SimpleModelFitter()
	{
		delete quality_;
	}


  SimpleModelFitter::SimpleModelFitter(const SimpleModelFitter& rhs)
    : BaseModelFitter(rhs),
    	quality_(0)
  {
    updateMembers_();
  }
  
  SimpleModelFitter& SimpleModelFitter::operator= (const SimpleModelFitter& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseModelFitter::operator=(rhs);
    
    updateMembers_();
    
    return *this;
  }

	void SimpleModelFitter::updateMembers_()
	{
		if (quality_) delete quality_;
		quality_ = Factory<BaseQuality>::create(param_.getValue("quality:type"));
	}

  DFeature<2> SimpleModelFitter::fit(const IndexSet& set) throw (UnableToFit)
	{
		
		// not enough peaks to fit
		if (set.size() < static_cast<Size>(param_.getValue("min_num_peaks:extended")))
		{
			String mess = String("Skipping feature, IndexSet size too small: ") + set.size();
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__, "UnableToFit-IndexSet", mess.c_str());
		}
		
		quality_->setTraits(traits_);
		ModelDescription<2> model_desc;
		double quality = 0.0;
		double max_quality = -std::numeric_limits<double>::max();
		
		float stdev = param_.getValue("isotope_model:stdev:first");
		float last = param_.getValue("isotope_model:stdev:last");
		float step = param_.getValue("isotope_model:stdev:step");
		
		// Calculate statistics
		mz_stat_.update( 	IntensityIterator(set.begin(),traits_),
											IntensityIterator(set.end(),traits_),
											MzIterator(set.begin(),traits_) );

		rt_stat_.update ( IntensityIterator(set.begin(),traits_),
										 IntensityIterator(set.end(),traits_),
										 RtIterator(set.begin(),traits_) );

		// Calculate bounding box
		IndexSet::const_iterator it=set.begin();
		min_ = max_ = traits_->getPeakPos(*it);
		for ( ++it; it!=set.end(); ++it )
		{
			CoordinateType tmp = traits_->getPeakMz(*it);
			if (min_[MZ] > tmp) min_[MZ] = tmp;
			if (max_[MZ] < tmp) max_[MZ] = tmp;
			tmp = traits_->getPeakRt(*it);
			if (min_[RT] > tmp) min_[RT] = tmp;
			if (max_[RT] < tmp) max_[RT] = tmp;
		}
		
		double const tolerance_stdev_box = param_.getValue("tolerance_stdev_bounding_box");
		stdev_mz_ = sqrt ( mz_stat_.variance() ) * tolerance_stdev_box;
		min_[MZ] -= stdev_mz_;
		max_[MZ] += stdev_mz_;

		stdev_rt1_ = sqrt ( rt_stat_.variance1() ) * tolerance_stdev_box;
		stdev_rt2_ = sqrt ( rt_stat_.variance2() ) * tolerance_stdev_box;
		min_[RT] -= stdev_rt1_;
		max_[RT] += stdev_rt2_;

		// Test charges and stdevs
		const int first_model = param_.getValue("mz:model_type:first");
		const int last_model = param_.getValue("mz:model_type:last");
		for ( ; stdev <= last; stdev += step)
		{
			for (int mz_fit_type = first_model; mz_fit_type <= last_model; ++mz_fit_type)
			{
				quality = fit_(set, static_cast<MzFitting>(mz_fit_type), BIGAUSS, stdev);
				if (quality > max_quality)
				{
					max_quality = quality;
					model_desc = ModelDescription<2>(&model2D_);
				}
			}
		}
	
		// model with highest correlation
		ProductModel<2>* final = static_cast< ProductModel<2>* >(model_desc.createModel());
		
		// model_desc.createModel() returns 0 if class model_desc is not initialized
		// in this case something went wrong during the modelfitting and we stop.
		if (! final) 
		{
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,
										"UnableToFit-BadQuality","Zero quality after fitting. Skipping this feature");			
			delete final;   
		}

		// find peak with highest predicted intensity to use as cutoff
		float model_max = 0;
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			float model_int = final->getIntensity(traits_->getPeakPos(*it));
			if (model_int>model_max) model_max = model_int;
		}
		final->setCutOff( model_max * float(param_.getValue("intensity_cutoff_factor")));

		// Cutoff low intensities wrt to model maximum -> cutoff independent of scaling
		IndexSet model_set;
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it) 
		{
			if ( final->isContained( traits_->getPeakPos(*it) ))
			{
				model_set.insert(*it);
			}
			else		// free dismissed peak via setting the appropriate flag
			{
				traits_->getPeakFlag(*it) = FeaFiTraits::UNUSED;
			}
		}
		// Print number of selected peaks after cutoff
		std::cout << " Selected " << model_set.size() << " from " << set.size() << " peaks.\n";
	
		// not enough peaks left for feature
		if (model_set.size() < static_cast<Size>(param_.getValue("min_num_peaks:final")))
		{
			delete final;
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,
						"UnableToFit-FinalSet",
						String("Skipping feature, IndexSet size after cutoff too small: ") + model_set.size() );
		}
		max_quality = quality_->evaluate(model_set, *final); // recalculate quality after cutoff
		
		std::cout << "P-value : " << quality_->getPvalue() << std::endl;
			
		// fit has too low quality or fit was not possible i.e. because of zero stdev
		if (max_quality < static_cast<float>(param_.getValue("quality:minimum")))
		{
			delete final;
			String mess = String("Skipping feature, correlation too small: ") + max_quality;
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__, "UnableToFit-Correlation", mess.c_str());
		}

		// Calculate intensity scaling
		float model_sum = 0;
		float data_sum = 0;
		float data_max = 0;
		for (IndexSet::const_iterator it=model_set.begin(); it!=model_set.end(); ++it)
		{
			float model_int = final->getIntensity(traits_->getPeakPos(*it));
			model_sum += model_int;
			data_sum += traits_->getPeakIntensity(*it);
			if (traits_->getPeakIntensity(*it) > data_max) data_max = traits_->getPeakIntensity(*it);
		}

		// fit has too low quality or fit was not possible i.e. because of zero stdev
		if (model_sum == 0)
		{
			delete final;
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,
						"UnableToFit-ZeroSum","Skipping feature, model_sum zero.");
		}
		
		final->setScale(data_max/model_max);	// use max quotient instead of sum quotient
	
		// Build Feature
		// The feature coordinate in rt dimension is given
		// by the centroid of the rt model whereas the coordinate
		// in mz dimension is equal to the monoisotopic peak.
		DFeature<2> f;
		f.setModelDescription( ModelDescription<2>(final) );
		f.setOverallQuality(max_quality);
		f.getPosition()[RT] = static_cast<InterpolationModel<>*>(final->getModel(RT))->getCenter();
		f.getPosition()[MZ] = static_cast<InterpolationModel<>*>(final->getModel(MZ))->getCenter();
		
		// set feature charge												
		if (final->getModel(MZ)->getName() == "IsotopeModel")
		{
			f.setCharge(static_cast<IsotopeModel*>(final->getModel(MZ))->getCharge());
		}
		// if we used a simple Gaussian model to fit the feature, we can't say anything about
		// its charge state. The value 0 indicates that charge state is undetermined.
		else 
		{
			f.setCharge(0);		
		}
				
		int const intensity_choice = param_.getValue("feature_intensity_sum");
		double feature_intensity = 0.0;
		
		if (intensity_choice == 1)
		{
			// intensity of the feature is the sum of all included data points
			for (IndexSet::const_iterator it=model_set.begin(); it!=model_set.end(); ++it) 
			{
				feature_intensity += traits_->getPeakIntensity(*it);
			}
		}
		else
		{
			// feature intensity is the maximum intensity of all peaks
			for (IndexSet::const_iterator it=model_set.begin(); it!=model_set.end(); ++it) 
			{
				if (traits_->getPeakIntensity(*it) > feature_intensity)
					feature_intensity = traits_->getPeakIntensity(*it);
			}	
		}
		
		f.setIntensity(feature_intensity);
		traits_->addConvexHull(model_set, f);
		
		std::cout << Date::now() << " Feature " << counter_ << ": (" << f.getPosition()[RT]
							<< "," << f.getPosition()[MZ] << ") Qual.:" << max_quality << "\n";
		
		f.getQuality(RT) = quality_->evaluate(model_set, *final->getModel(RT), RT );
		f.getQuality(MZ) = quality_->evaluate(model_set, *(static_cast<InterpolationModel<>*>(final->getModel(MZ)) ),MZ );

		// save meta data in feature for TOPPView
		stringstream s;
		s <<  "Feature #" << counter_ << ", +" << f.getCharge() << ", " << set.size() << "->" << model_set.size() 
			<< ", Corr: (" << max_quality << "," << f.getQuality(RT) << "," << f.getQuality(MZ) << ")";
		f.setMetaValue(3,s.str());
		
		#ifdef DEBUG_FEATUREFINDER
		// write debug output
		CoordinateType rt = f.getPosition()[RT];
		CoordinateType mz = f.getPosition()[MZ];
		
		// write feature model 
		String fname = String("model") + counter_ + "_" + rt + "_" + mz;
		ofstream file(fname.c_str()); 
		for (IndexSet::const_iterator it=model_set.begin(); it!=model_set.end(); ++it) 
		{
			FeaFiTraits::PositionType2D p = traits_->getPeakPos(*it);
			if ( final->isContained(p) )
			{
				file << p[RT] << " " << p[MZ] << " " << final->getIntensity(p) << "\n";						
			}
		}
		
		// wrote peaks remaining after model fit
		fname = String("feature") + counter_ + "_" + rt + "_" + mz;
		ofstream file2(fname.c_str()); 
		for (IndexSet::const_iterator it=model_set.begin(); it!=model_set.end(); ++it) 
		{
			FeaFiTraits::PositionType2D p = traits_->getPeakPos(*it);
			if ( final->isContained(p) )
			{
				file2 << p[RT] << " " << p[MZ] << " " << traits_->getPeakIntensity(*it) << "\n";						
			}
		}
		file2.close();
		#endif
		++counter_;

		delete final;
	  return f;
	}

	double SimpleModelFitter::fit_(const IndexSet& set, MzFitting mz_fit, RtFitting rt_fit, Coordinate isotope_stdev)
	{
		const Coordinate interpolation_step_mz = param_.getValue("mz:interpolation_step");
		const Coordinate interpolation_step_rt = param_.getValue("rt:interpolation_step");

		// Build Models
		InterpolationModel<>* mz_model;
		if (mz_fit==MZGAUSS)
		{
			mz_model = new GaussModel();
			mz_model->setInterpolationStep(interpolation_step_mz);
			static_cast<GaussModel*>(mz_model)->setParam(mz_stat_, min_[MZ], max_[MZ]);
		}
		else
		{
			mz_model = new IsotopeModel();
			Param iso_param = param_.copy("isotope_model:",true);
			iso_param.remove("stdev");
			mz_model->setParameters(iso_param);
			mz_model->setInterpolationStep(interpolation_step_mz);
			static_cast<IsotopeModel*>(mz_model)->setParam(mz_stat_.mean(), mz_fit, isotope_stdev);
		}
		InterpolationModel<>* rt_model;
		if (rt_fit==RTGAUSS)
		{
			rt_model = new GaussModel();
			rt_model->setInterpolationStep(interpolation_step_rt);
			static_cast<GaussModel*>(rt_model)->setParam(rt_stat_, min_[RT], max_[RT]);
		}else
		{
			rt_model = new BiGaussModel();
			rt_model->setInterpolationStep(interpolation_step_rt);
			static_cast<BiGaussModel*>(rt_model)->setParam(rt_stat_.mean(), rt_stat_.variance1(),
					 rt_stat_.variance2(),	min_[RT], max_[RT]);
		}
		model2D_.setModel(MZ, mz_model).setModel(RT, rt_model);

		double res;
		res = fitOffset_(mz_model, set, stdev_mz_, stdev_mz_, interpolation_step_mz);
		res = fitOffset_(rt_model, set, stdev_rt1_, stdev_rt2_, interpolation_step_rt);
		return res;
	}


	double SimpleModelFitter::fitOffset_(	InterpolationModel<>* model, const IndexSet& set, const double stdev1,  const double stdev2, const Coordinate offset_step)
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

}



