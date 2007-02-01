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
// $Maintainer: Clemens Groepl, Marcel Grunert $
// --------------------------------------------------------------------------


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedModelFitter.h>
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
		// Helper struct for ExtendedModelFitter
		struct ExpFitPolyData
		{
			size_t n;
			string profile;
		};
	}

	using namespace Internal;
	using namespace std;

	//positions and signal values
	std::vector<double> positionsDC_;
	std::vector<double> signalDC_;

	ExtendedModelFitter::ExtendedModelFitter()
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
		defaults_.setValue("feature_intensity_sum",1);
		defaults_.setValue("rt:interpolation_step",0.2f);
		defaults_.setValue("rt:max_iteration",500);
		defaults_.setValue("rt:deltaAbsError",0.0001);
		defaults_.setValue("rt:deltaRelError",0.0001);
		defaults_.setValue("rt:profile","EMG");
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

	ExtendedModelFitter::~ExtendedModelFitter()
	{
		delete quality_;
	}

  ExtendedModelFitter::ExtendedModelFitter(const ExtendedModelFitter& rhs)
    : BaseModelFitter(rhs),
    	quality_(0)
  {
    updateMembers_();
  }
  
  ExtendedModelFitter& ExtendedModelFitter::operator= (const ExtendedModelFitter& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseModelFitter::operator=(rhs);
    
    updateMembers_();
    
    return *this;
  }

	void ExtendedModelFitter::updateMembers_()
	{
		if (quality_) delete quality_;
		quality_ = Factory<BaseQuality>::create(param_.getValue("quality:type"));
		
		max_iteration_ = param_.getValue("rt:max_iteration");
		eps_abs_ = param_.getValue("rt:deltaAbsError");
		eps_rel_ = param_.getValue("rt:deltaRelError");
		profile_ = (string)param_.getValue("rt:profile");
	}

	DFeature<2> ExtendedModelFitter::fit(const IndexSet& set) throw (UnableToFit)
	{
		// not enough peaks to fit
		if (set.size() < (UnsignedInt)(param_.getValue("min_num_peaks:extended")))
		{
			String mess = String("Skipping feature, IndexSet size too small: ") + set.size();
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,"UnableToFit-IndexSet", mess.c_str());
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
		IndexSetIter it=set.begin();
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

		/// create a vector with RT-values & Intensity
		/// compute the parameters (intial values) for the EMG & Gauss function and finally,
		/// optimize the parameters with Levenberg-Marquardt algorithms
		if (profile_=="LmaGauss" || profile_=="EMG" || profile_=="LogNormal")
		{
			setData(set);
			optimize();
		}

		/// Test different charges and stdevs
		const int first_model = param_.getValue("mz:model_type:first");
		const int last_model = param_.getValue("mz:model_type:last");
		for ( ; stdev <= last; stdev += step)
		{
			for (int mz_fit_type = first_model; mz_fit_type <= last_model; ++mz_fit_type)
			{
				if (profile_=="LmaGauss")
					quality = fit_(set, static_cast<MzFitting>(mz_fit_type), LMAGAUSS, stdev);
				else if (profile_=="EMG")
					quality = fit_(set, static_cast<MzFitting>(mz_fit_type), EMGAUSS, stdev);
				else if (profile_=="LogNormal")
					quality = fit_(set, static_cast<MzFitting>(mz_fit_type), LOGNORMAL, stdev);
				else
					quality = fit_(set, static_cast<MzFitting>(mz_fit_type), BIGAUSS, stdev);

				if (quality > max_quality)
				{
					max_quality = quality;
					model_desc = ModelDescription<2>(&model2D_);
				}
			}
		}

		// model with highest correlation
		ProductModel<2>* final = dynamic_cast< ProductModel<2>* >(model_desc.createModel());
		
		// model_desc.createModel() returns 0 if class model_desc is not initialized
		// in this case something went wrong during the modelfitting and we stop.
		if (! final)
		{
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__, "UnableToFit-BadQuality","Zero quality after fitting. Skipping this feature");
			delete final;
		}

		// find peak with highest predicted intensity to use as cutoff
		float model_max = 0;
		for (IndexSetIter it=set.begin(); it!=set.end(); ++it)
		{
			float model_int = final->getIntensity(traits_->getPeakPos(*it));
			if (model_int>model_max) model_max = model_int;
		}
		final->setCutOff( model_max * float(param_.getValue("intensity_cutoff_factor")));

		// Cutoff low intensities wrt to model maximum -> cutoff independent of scaling
		IndexSet model_set;
		for (IndexSetIter it=set.begin(); it!=set.end(); ++it) {
			if ( final->isContained(traits_->getPeakPos(*it)) )
			{
				model_set.insert(*it);
			}else		// free dismissed peak via setting the appropriate flag
			{
				traits_->getPeakFlag(*it) = FeaFiTraits::UNUSED;
			}
		}
		// Print number of selected peaks after cutoff
		std::cout << " Selected " << model_set.size() << " from " << set.size() << " peaks.\n";

		// not enough peaks left for feature
		if (model_set.size() < (UnsignedInt)(param_.getValue("min_num_peaks:final")))
		{
			delete final;
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,
												"UnableToFit-FinalSet",
												String("Skipping feature, IndexSet size after cutoff too small: ") + model_set.size() );
		}
 		max_quality = quality_->evaluate(model_set, *final); // recalculate quality after cutoff

		std::cout << "P-value : " << quality_->getPvalue() << std::endl;
		
		// fit has too low quality or fit was not possible i.e. because of zero stdev
		if (max_quality < (float)(param_.getValue("quality:minimum")))
		{
			delete final;
			String mess = String("Skipping feature, correlation too small: ") + max_quality;
			throw UnableToFit(__FILE__, __LINE__,__PRETTY_FUNCTION__,"UnableToFit-Correlation", mess.c_str());
		}

		// Calculate intensity scaling
		float model_sum = 0;
		float data_sum = 0;
		float data_max = 0;
		for (IndexSetIter it=model_set.begin(); it!=model_set.end(); ++it)
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
		f.getPosition()[RT] = dynamic_cast<InterpolationModel<>*>(final->getModel(RT))->getCenter();
		f.getPosition()[MZ] = dynamic_cast<InterpolationModel<>*>(final->getModel(MZ))->getCenter();
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

		int const intensity_choice = param_.getValue("feature_intensity_sum");
		double feature_intensity = 0.0;
		
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

		std::cout << Date::now() << " Feature " << counter_
							<< ": (" << f.getPosition()[RT]
							<< "," << f.getPosition()[MZ] << ") Qual.:"
							<< max_quality << "\n";


		f.getQuality(RT) = quality_->evaluate(model_set, *final->getModel(RT), RT );
		f.getQuality(MZ) = quality_->evaluate(model_set, *(dynamic_cast<InterpolationModel<>*>(final->getModel(MZ)) ), MZ );

		// save meta data in feature for TOPPView
		stringstream meta ;
		meta << "Feature #" << counter_ << ", +"	<< f.getCharge() << ", " << set.size() << "->" << model_set.size() 
				 << ", Corr: (" << max_quality << ","  << f.getQuality(RT) << "," << f.getQuality(MZ) << ")";
		f.setMetaValue(3,meta.str());

		#ifdef DEBUG_FEATUREFINDER
		// write debug output
		CoordinateType rt = f.getPosition()[RT];
		CoordinateType mz = f.getPosition()[MZ];
		
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
		file.close();
		#endif
		
		++counter_;

		delete final;

 		return f;
	}

	double ExtendedModelFitter::fit_(const IndexSet& set, MzFitting mz_fit, RtFitting rt_fit,
																	 Coordinate isotope_stdev)
	{

		const Coordinate interpolation_step_mz = param_.getValue("mz:interpolation_step");
		const Coordinate interpolation_step_rt = param_.getValue("rt:interpolation_step");

		// Build Models
		InterpolationModel<>* mz_model;
		if (mz_fit==MZGAUSS)
		{
			mz_model = new GaussModel();
			mz_model->setInterpolationStep(interpolation_step_mz);
			dynamic_cast<GaussModel*>(mz_model)->setParam(mz_stat_, min_[MZ], max_[MZ]);
		}
		else
		{
			// new model
			mz_model = new IsotopeModel();
			Param iso_param = param_.copy("isotope_model:",true);
			iso_param.remove("stdev");
			mz_model->setParameters(iso_param);
			mz_model->setInterpolationStep(interpolation_step_mz);
			dynamic_cast<IsotopeModel*>(mz_model)->setParam(mz_stat_.mean(), mz_fit, isotope_stdev);
		}
		InterpolationModel<>* rt_model;
		if (rt_fit==RTGAUSS)
		{
			rt_model = new GaussModel();
			rt_model->setInterpolationStep(interpolation_step_rt);
			dynamic_cast<GaussModel*>(rt_model)->setParam(rt_stat_, min_[RT], max_[RT]);

		}
		else if (rt_fit==LMAGAUSS)
		{
			rt_model = new LmaGaussModel();
			rt_model->setInterpolationStep(interpolation_step_rt);
			dynamic_cast<LmaGaussModel*>(rt_model)->setParam(rt_stat_, scale_factor_, standard_deviation_, expected_value_, min_[RT], max_[RT]);
		}
		else if (rt_fit==EMGAUSS)
		{
			rt_model = new EmgModel();
			rt_model->setInterpolationStep(interpolation_step_rt);
			dynamic_cast<EmgModel*>(rt_model)->setParam(rt_stat_, height_, width_, symmetry_, retention_, min_[RT], max_[RT]);
		}
		else if (rt_fit==LOGNORMAL)
		{
			rt_model = new LogNormalModel();
			rt_model->setInterpolationStep(interpolation_step_rt);
			dynamic_cast<LogNormalModel*>(rt_model)->setParam(rt_stat_, height_, width_, symmetry_, retention_, r_, min_[RT], max_[RT]);
		}
		else
		{
			rt_model = new BiGaussModel();
			rt_model->setInterpolationStep(interpolation_step_rt);
			dynamic_cast<BiGaussModel*>(rt_model)->setParam(rt_stat_.mean(), rt_stat_.variance1(), rt_stat_.variance2(), min_[RT], max_[RT]);
		}

		model2D_.setModel(MZ, mz_model).setModel(RT, rt_model);

		double res;
		res = fitOffset_(mz_model, set, stdev_mz_, stdev_mz_, interpolation_step_mz);
		if (profile_!="LmaGauss" && profile_!="EMG" && profile_!="LogNormal")
			res = fitOffset_(rt_model, set, stdev_rt1_, stdev_rt2_, interpolation_step_rt);

		return res;

	}


	double ExtendedModelFitter::fitOffset_(	InterpolationModel<>* model,
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
	void ExtendedModelFitter::setData(const IndexSet& set)
	{
		// start rt-value with intensity 0.0

		positionsDC_.push_back(min_[RT]);
		positionsDC_.push_back(min_[RT]+1);
		signalDC_.push_back(0.);
		signalDC_.push_back(0.);

		// sum over all intensities
		double sum = 0.0;

		//String fname = String("feature") + counter_ + "_orginal_" + profile_;
		//ofstream orgFile(fname.c_str());

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
			if (positionsDC_.empty() || positionsDC_.back()!=position) {
				positionsDC_.push_back(position);
				signalDC_.push_back(signal);
			}
			else {
				signal += signalDC_.back();
				signalDC_.pop_back();
				signalDC_.push_back(signal);
			}
		}

		// end rt-value with intensity 0.0

		positionsDC_.push_back(max_[RT]-1);
		positionsDC_.push_back(max_[RT]);
		signalDC_.push_back(0.);
		signalDC_.push_back(0.);

		// calculate the median
		int median = 0;
		float count = 0.0;
		for (size_t current_point=0; current_point<positionsDC_.size();current_point++)
		{
			count += signalDC_[current_point];
			if (count <= sum/2)
				median = current_point;
		}

		double sumS = 0.0;
		for (size_t current_point=0; current_point<positionsDC_.size();current_point++)
		{
			sumS += pow((positionsDC_[current_point] - positionsDC_[median]),2);
		}

		// calculate the stardard deviation
		deviation_ = sqrt(sumS/(positionsDC_.size() - 1));

		// calculate the heigth of the peak
		height_ = signalDC_[median];

		// calculate the width of the peak
		// rt-values with intensity zero are not allowed for calculation of the width
		width_ = abs(positionsDC_[positionsDC_.size()-3] - positionsDC_[2]);

		// calculate retention time
		retention_ = positionsDC_[median];

		// calculate the symmetry (fronted peak: s<1 , tailed peak: s>1)
		symmetry_ = abs(positionsDC_.back() - positionsDC_[median])/abs(positionsDC_[median] - positionsDC_.front());

		// optimize the symmetry
		if (profile_=="LogNormal") {
			// The computations can lead to an overflow error at very low values of symmetry (s~0).
			if (symmetry_<=0.8)
				symmetry_=0.8;

			if (symmetry_>=1.5)
				symmetry_=1.4;
		}
		else
		{
			// The computations can lead to an overflow error at very low values of symmetry (s~0). For s~5 the parameter can be aproximized by the Levenberg-Marquardt argorithms.
			if (symmetry_<1)
				symmetry_+=5;
			else
				symmetry_=5;
		}


		// calculate the parameter r of the log normal function
		// r is the ratio between h and the height at which w and s are computed
		r_ = 1.5;
	}

	//Evaluation of the target function for nonlinear optimization.
	int residualDC(const gsl_vector* x, void* params , gsl_vector* f)
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
				double t = positionsDC_[i];

				Yi=(1/(sqrt(2*M_PI)*normal_s))*exp(-((t-normal_m)*(t-normal_m))/(2*normal_s*normal_s))*normal_scale;

				gsl_vector_set(f, i, (Yi - signalDC_[i]));
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
					double t = positionsDC_[i];

					// Simplified EMG
					Yi=(h*w/s)*sqrt(2*M_PI)*exp(((w*w)/(2*s*s))-((t-z)/s))/(1+exp((-2.4055/sqrt(2))*(((t-z)/w)-w/s)));

					gsl_vector_set(f, i, (Yi - signalDC_[i]));
				}
			}
			/// log normal
			else
			{
				double h = gsl_vector_get(x,0);
				double w = gsl_vector_get(x,1);
				double s = gsl_vector_get(x,2);
				double z = gsl_vector_get(x,3);
				double r = gsl_vector_get(x,4);

				double Yi = 0.0;

				for (size_t i = 0; i < n; i++)
				{
					double t = positionsDC_[i];

					Yi = h*exp(-log(r)/(log(s)*log(s))*pow(log((t-z)*(s*s-1)/(w*s)+1),2));

					gsl_vector_set(f, i, (Yi - signalDC_[i]));
				}
			}
		}

		return GSL_SUCCESS;
	}

	/** Compute the Jacobian of the residual, where each row of the matrix corresponds to a
	 *  point in the data. */
	int jacobianDC(const gsl_vector* x, void* params, gsl_matrix* J)
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
				double t = positionsDC_[i];

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
					double t = positionsDC_[i];

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
			/// log normal function
			else
			{
				double h = gsl_vector_get(x,0);
				double w = gsl_vector_get(x,1);
				double s = gsl_vector_get(x,2);
				double z = gsl_vector_get(x,3);
				double r = gsl_vector_get(x,4);

				double derivative_height, derivative_width, derivative_symmetry, derivative_retention,  derivative_r = 0.0;

				// iterate over all points of the signal
				for (size_t i = 0; i < n; i++)
				{
					double t = positionsDC_[i];

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
					gsl_matrix_set(J, i, 4, derivative_r);
				}
			}
		}

		return GSL_SUCCESS;
	}

	// Driver function for the evaluation of function and jacobian.
	int evaluateDC(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
	{
		residualDC(x, params, f);
		jacobianDC(x, params, J);

		return GSL_SUCCESS;
	}

	// perform a nonlinear optimization
	void ExtendedModelFitter::optimize()
	{
		const gsl_multifit_fdfsolver_type *T;
		gsl_multifit_fdfsolver *s;

		int status;
		size_t iter = 0;
		const size_t n = positionsDC_.size();

		// number of parameter to be optimize
		unsigned int p = 0;
		if (profile_=="LmaGauss") p = 3;
		else if (profile_=="LogNormal") p = 5;
		else p = 4;

		gsl_matrix *covar = gsl_matrix_alloc(p,p);
		gsl_multifit_function_fdf f;

		double x_init_normal[3] = { deviation_, rt_stat_.mean(), height_ };
		double x_init_emg[4] = { height_, width_, symmetry_, retention_ };
		double x_init_lognormal[5] = { height_, width_, symmetry_, retention_, r_ };

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
		f.f = &residualDC;
		f.df = &jacobianDC;
		f.fdf = &evaluateDC;
		f.n = n;
		f.p = p;
		f.params = &d;

		T = gsl_multifit_fdfsolver_lmsder;
		s = gsl_multifit_fdfsolver_alloc(T,n,p);
		gsl_multifit_fdfsolver_set(s, &f, &x.vector);

#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
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
				printf ("before loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f % 15.8f |f(x)| = %g\n", iter,
								gsl_vector_get(s->x,0),
								gsl_vector_get(s->x,1),
								gsl_vector_get(s->x,2),
								gsl_vector_get(s->x,3),
								gsl_vector_get(s->x,4),
								gsl_blas_dnrm2(s->f));
			}
		}
#endif

		// this is the loop for fitting
		do
		{
			iter++;
			status = gsl_multifit_fdfsolver_iterate (s);

#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
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
					printf ("in loop iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f % 15.8f |f(x)| = %g\n", iter,
									gsl_vector_get(s->x,0),
									gsl_vector_get(s->x,1),
									gsl_vector_get(s->x,2),
									gsl_vector_get(s->x,3),
									gsl_vector_get(s->x,4),
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

#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
		gsl_matrix_fprintf(stdout, covar, "%g");
#endif

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
		cout << "status: " << gsl_strerror(status) << endl;
		cout << "symmetry: " << symmetry_ << endl;
#endif
		if (profile_=="LmaGauss")
		{
#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
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
#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
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
#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
				printf("h = %.5f +/- %.5f\n", FIT(0), ERR(0));
				printf("w = %.5f +/- %.5f\n", FIT(1), ERR(1));
				printf("s = %.5f +/- %.5f\n", FIT(2), ERR(2));
				printf("z = %.5f +/- %.5f\n", FIT(3), ERR(3));
				printf("r = %.5f +/- %.5f\n", FIT(4), ERR(4));
#endif
				height_     = FIT(0);
				width_      = FIT(1);
				symmetry_   = FIT(2);
				retention_  = FIT(3);
				r_ 	    = FIT(4);
			}
		}

#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
		{
			// chi-squared value
			double chi = gsl_blas_dnrm2(s->f);
			printf("chisq/dof = %g\n", pow(chi, 2.0)/ (n-p));
		}
#endif

		// function free all memory associated with the solver s
		gsl_multifit_fdfsolver_free(s);

#ifdef DEBUG_EXTENDEDMODELFITTER_LMA
		for (size_t current_point=0; current_point<positionsDC_.size();current_point++)
			std::cout << positionsDC_[current_point] << " " << signalDC_[current_point] << std::endl;

		cout << "" << endl;
		cout << "*** parameter for optimization ***" << endl;
		cout << "       height:  " << height_ << endl;
		cout << "        width:  " << width_ << endl;
		cout << "     symmetry:  " << symmetry_ << endl;
		cout << "    retention:  " << retention_ << endl;
		cout << "std.deviation:  " << deviation_ << endl;
		cout << "max_iteration:  " << max_iteration_ << endl;
		cout << "      eps_abs:  " << eps_abs_ << endl;
		cout << "      eps_rel:  " << eps_rel_ << endl;
		cout << "      profile:  " << profile_ << endl;
		cout << "" << endl;
#endif

		positionsDC_.clear();
		signalDC_.clear();
	}

}
