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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>

#include <OpenMS/SIMULATION/IsotopeModelGeneral.h>
#include <OpenMS/SIMULATION/MixtureModel.h>
#include <OpenMS/SIMULATION/ElutionModel.h>

#include <vector>
using std::vector;

namespace OpenMS {

  /**
   TODO:
    * add removeDuplicate Points
    * review Ole's methods for improvment/changes
    * BIG TODO: initialize members!!!!
  */
  RawMSSignalSimulation::RawMSSignalSimulation(const gsl_rng * random_generator)
  : DefaultParamHandler("RawSignalSimulation"), rnd_gen_(random_generator)
  {
    setDefaultParams_();
    updateMembers_();
  }


  RawMSSignalSimulation::RawMSSignalSimulation()
    : DefaultParamHandler("RawSignalSimulation")
  {
    setDefaultParams_();
    updateMembers_();
  }

  RawMSSignalSimulation::RawMSSignalSimulation(const RawMSSignalSimulation& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
  }

  RawMSSignalSimulation& RawMSSignalSimulation::operator = (const RawMSSignalSimulation& source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
    return *this;
  }

  RawMSSignalSimulation::~RawMSSignalSimulation()
  {}

  void RawMSSignalSimulation::setDefaultParams_()
  {
    // noise params
    // m/z error
    defaults_.setValue("mz:error_mean",0.0,"Average systematic m/z error (Da)");
    defaults_.setValue("mz:error_stddev",0.0,"Standard deviation for m/z errors. Set to 0 to disable simulation of m/z errors.");

    // maximal size of map in mz dimension
    defaults_.setValue("mz:upper_measurement_limit",2500.0,"Upper m/z detecter limit.");
    defaults_.setValue("mz:lower_measurement_limit",200.0,"Lower m/z detecter limit.");
    // mz sampling rate
    defaults_.setValue("mz:sampling_rate",0.12,"MS hardware resolution (e.g. bin size in m/z).");

    // intensity error
    defaults_.setValue("int:error_mean",0,"Average systematic intensity error.");
    defaults_.setValue("int:error_stddev",0.0,"Standard deviation for peak intensities (relative to peak height). Set to 0 to disable intensity errors.");

    // peak and instrument parameter
    defaults_.setValue("peak_fwhm",0.5,"FWHM (full width at half maximum) of simulated peaks (Da).");

    // shot noise
    defaults_.setValue("noise_rate",0,"Poisson rate of shot noise. Set to 0 to disable simulation of shot noise.");
    defaults_.setValue("noise_int_mean",50.0,"Shot noise intensity mean.");

    // baseline
    defaults_.setValue("baseline_scaling",0.0,"Scale of baseline (zero disables baseline simulation)");

    defaultsToParam_();
  }

  void RawMSSignalSimulation::updateMembers_()
  {
    double tmp    =  param_.getValue("peak_fwhm");
    peak_std_     = (tmp / 2.355);			// Approximation for Gaussian-shaped signals
    mz_sampling_rate_ = param_.getValue("mz:sampling_rate");

		// TODO: this needs to be known to IonizationsSim as well...
		// or every signal outside the range must be discarded here (favoring the first - use param-mirroring in MSSim)
    maximal_mz_measurement_limit_ = param_.getValue("mz:upper_measurement_limit");
    minimal_mz_measurement_limit_ = param_.getValue("mz:lower_measurement_limit");

    mz_error_mean_    = param_.getValue("mz:error_mean");
    mz_error_stddev_  = param_.getValue("mz:error_stddev");

    intensity_error_mean_   = param_.getValue("int:error_mean");
    intensity_error_stddev_ = param_.getValue("int:error_stddev");

  }

  void RawMSSignalSimulation::generateRawSignals(FeatureMapSim & features, MSSimExperiment & experiment)
  {
    if(experiment.size() == 1)
    {
      for(FeatureMap< >::iterator feature_it = features.begin();
          feature_it != features.end();
          ++feature_it)
      {
        add1DSignal_(*feature_it,experiment);
      }
    }
    else
    {
      for(FeatureMap< >::iterator feature_it = features.begin();
          feature_it != features.end();
          ++feature_it)
      {
        add2DSignal_(*feature_it, experiment);
      }
      addShotNoise_(experiment);
      addBaseLine_(experiment);
    }
    compressSignals_(experiment);
  }

  void RawMSSignalSimulation::add1DSignal_(Feature & active_feature, MSSimExperiment & experiment)
  {
    Param p1;

    // was: 3000 TODO: ???? why 1500
    SimIntensityType scale = active_feature.getIntensity() * 1500;
    // TODO: ahhh. Does that mean every time I use this funtion (luckily we only do it once, mean_scaling_ is increased?)
    mean_scaling_ += scale;
    ++ion_count_;

    SimChargeType charge = active_feature.getCharge();
    EmpiricalFormula feature_ef = active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula();
    SimCoordinateType mz = active_feature.getMZ();

    // don't show ions with m/z higher than the MS detection limit
    if (mz > maximal_mz_measurement_limit_ || mz < minimal_mz_measurement_limit_)
    {
      return; // ignore the current feature
    }

    p1.setValue("statistics:mean", active_feature.getMZ() );
    p1.setValue("interpolation_step", 0.001);
    p1.setValue("isotope:stdev", peak_std_);
    p1.setValue("intensity_scaling", scale);
    p1.setValue("charge", charge);

    IsotopeModelGeneral isomodel;
    isomodel.setSamples(feature_ef);
    isomodel.setParameters(p1);

    SimCoordinateType mz_start = isomodel.getInterpolation().supportMin();
    SimCoordinateType mz_end = isomodel.getInterpolation().supportMax();
        
    samplePeptideModel1D_(isomodel, mz_start, mz_end, experiment, active_feature);
  }

  void RawMSSignalSimulation::add2DSignal_(Feature & active_feature, MSSimExperiment & experiment)
  {
    // was: 3000 TODO: ???? why 1500
    SimIntensityType scale = active_feature.getIntensity() * 1500;
    mean_scaling_ += scale;
    ++ion_count_;

    SimCoordinateType mz = active_feature.getMZ();

    // don't show ions with m/z higher than the MS detection limit
    if (mz > maximal_mz_measurement_limit_ || mz < minimal_mz_measurement_limit_)
    {
      return; // ignore the current feature
    }
    
    Param p1;
    p1.setValue("statistics:mean", active_feature.getMZ() );
    p1.setValue("interpolation_step", 0.001);
    p1.setValue("isotope:stdev", peak_std_);
    p1.setValue("charge", active_feature.getCharge());

    IsotopeModelGeneral* isomodel = new IsotopeModelGeneral();
    isomodel->setSamples(active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula());
    isomodel->setParameters(p1);

		if (experiment.size()<2)
		{
			throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, experiment.size());
		}
		DoubleReal rt_sampling_rate = experiment[1].getRT() - experiment[0].getRT();
		ElutionModel* elutionmodel = new ElutionModel();
    chooseElutionProfile_(elutionmodel, active_feature, scale, rt_sampling_rate, experiment);
    ProductModel<2> pm;
    pm.setModel(0, elutionmodel); // new'ed models will be deleted by the pm! no need to delete them manually
		pm.setModel(1, isomodel);			// new'ed models will be deleted by the pm! no need to delete them manually
    pm.setScale(scale);

    // start and end points of the sampling are entirely arbitrary
    // and should be modified at some point
    SimCoordinateType rt_start = elutionmodel->getInterpolation().supportMin();
    SimCoordinateType rt_end = elutionmodel->getInterpolation().supportMax();
    SimCoordinateType mz_start = isomodel->getInterpolation().supportMin();
    SimCoordinateType mz_end = isomodel->getInterpolation().supportMax();
    std::cout << ((DoubleReal)elutionmodel->getParameters().getValue("bounding_box:min")-rt_start) << " - " << ((DoubleReal)elutionmodel->getParameters().getValue("bounding_box:max")-rt_end) << "\n";
		std::cout <<  mz_start << " - " << mz_end << " [mz] \n";
		    
    // add peptide to global MS map
    samplePeptideModel2D_(pm, mz_start, mz_end, rt_start, rt_end, experiment, active_feature);
  }


  void RawMSSignalSimulation::samplePeptideModel1D_(const IsotopeModelGeneral & pm,
																										const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
																										MSSimExperiment & experiment, Feature & active_feature)
  {
    // start and end points of the sampling are entirely arbitrary
    // and should be modified at some point

    // (cg) commented this out since it cuts off fronted elution profiles!!
    // (ost) Why should this happen ?

    SimIntensityType intensity_sum = 0.0;

    std::cout << "Sampling at [mz] " << mz_start << ":" << mz_end << std::endl;

    /// TODO: think of better error checking

    SimPointType point;

    for (SimCoordinateType mz = mz_start; mz < mz_end; mz += mz_sampling_rate_)
    {
      //++it;
      point.setMZ(mz);
      point.setIntensity( pm.getIntensity( DPosition<1>( mz ) ) );

      if ( point.getIntensity() > 10.0)
      {
        // add m/z and itensity error (both Gaussian distributed)
        double it_err  = gsl_ran_gaussian(rnd_gen_, (point.getIntensity() * intensity_error_stddev_ ) ) + intensity_error_mean_ ;

        // this is a quick fix only, should be improved to prevent simulation of negative intensities
        if (it_err < 0.0) it_err = fabs(it_err);

        point.setIntensity( point.getIntensity( ) + it_err );

        double mz_err = gsl_ran_gaussian(rnd_gen_, mz_error_stddev_) + mz_error_mean_;
        point.setMZ( point.getMZ() + mz_err );

        intensity_sum += point.getIntensity();
        experiment[0].push_back(point);
      }
    }
    active_feature.setIntensity(intensity_sum);
  }


  void RawMSSignalSimulation::samplePeptideModel2D_(const ProductModel<2> & pm,
                                    const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
                                    SimCoordinateType rt_start, SimCoordinateType rt_end,
                                    MSSimExperiment & experiment, Feature & active_feature)
  {
    if (rt_start <=0) rt_start = 0;

    MSExperiment<Peak1D>::iterator exp_iter = experiment.RTBegin(rt_start);
    if(exp_iter == experiment.end() )
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 0);
    }


    SimIntensityType intensity_sum = 0.0;
    std::cout << "Sampling at [RT] " << rt_start << ":" << rt_end << " [mz] " << mz_start << ":" << mz_end << std::endl;

    SimPointType point;
    vector< DPosition<2> > points;
    
    Int start_scan = exp_iter - experiment.begin();
    Int end_scan  = -5;

		SimCoordinateType rt = rt_start;

    for (; rt < rt_end && exp_iter != experiment.end(); ++exp_iter)
    {
			rt = exp_iter->getRT();
			
      for (SimCoordinateType mz = mz_start; mz < mz_end; mz += mz_sampling_rate_)
      {
        //++it;

        point.setMZ(mz);
        point.setIntensity( pm.getIntensity( DPosition<2>( rt, mz) ) );

        if ( point.getIntensity() > 10.0)
        {
          // add m/z and intensity error (both Gaussian distributed)
          double it_err  = gsl_ran_gaussian(rnd_gen_, (point.getIntensity() * intensity_error_stddev_ ) ) + intensity_error_mean_ ;

          // this is a quick fix only, should be improved to prevent simulation of negative intensities
          if (it_err < 0.0) it_err = fabs(it_err);

          point.setIntensity( point.getIntensity( ) + it_err );

          double mz_err = gsl_ran_gaussian(rnd_gen_, mz_error_stddev_) + mz_error_mean_;
          point.setMZ( point.getMZ() + mz_err );

          intensity_sum += point.getIntensity();
          points.push_back( DPosition<2>( rt, mz) );		// store position
          exp_iter->push_back(point);
          //std::cout << "Sampling intensity: " << point.getIntensity() << std::endl;

          //update last scan affected
          end_scan = exp_iter - experiment.begin();
        }
      }
    }

    // This is a clear misuse of the Feature data structure
    // but at this point, we couldn't care less ;-)
    // TODO: this was currentFeature -> reimplement
    active_feature.setQuality(0,start_scan);
    active_feature.setQuality(1,end_scan);

    active_feature.setIntensity(intensity_sum);
    // store convex hull
    active_feature.getConvexHulls().clear();
    // adding ALL points of the feature, the ConvexHull2D assignment operator will take care of the rest
    active_feature.getConvexHulls().push_back(points);

  }

  void RawMSSignalSimulation::chooseElutionProfile_(ElutionModel*& elutionmodel, const Feature& feature, const double scale, const DoubleReal rt_sampling_rate, const MSSimExperiment & experiment)
  {
    SimCoordinateType f_rt = feature.getRT();
    DoubleReal f_symmetry = (DoubleReal) feature.getMetaValue("rt_symmetry");
    DoubleReal f_width = (DoubleReal) feature.getMetaValue("rt_width");

    // Exponentially modified Gaussian
    const DoubleReal decay_stretch = 5.0;
    // it doesn't matter what bounding box I set here, it always cuts of the fronting !!! :-(
    Param p;
    p.setValue("bounding_box:min", f_rt - decay_stretch*(f_width+fabs(f_symmetry)) );
    p.setValue("bounding_box:max", f_rt + decay_stretch*(f_width+fabs(f_symmetry)) );
    // WARNING: step used to be 'rt_sampling_rate / 3.0', but distortion is not part of RT sim, and thus only
    //          modelled 1:1
    p.setValue("interpolation_step", rt_sampling_rate );
    p.setValue("statistics:variance", 1.0);
    p.setValue("statistics:mean", f_rt);
    p.setValue("emg:height", scale);
    p.setValue("emg:width", f_width);
    p.setValue("emg:symmetry", f_symmetry);
    p.setValue("emg:retention", f_rt);
    elutionmodel->setParameters(p); // does the calculation

    //----------------------------------------------------------------------

    // Hack away constness :-P   We know what we want.
    ElutionModel::ContainerType &data = const_cast<ElutionModel::ContainerType &>(elutionmodel->getInterpolation().getData());

		SimCoordinateType rt_em_start = elutionmodel->getInterpolation().supportMin();
	
		// find scan in experiment at which our elution starts
		MSSimExperiment::ConstIterator exp_it = experiment.RTBegin(rt_em_start);
    for ( Size i = 1; (i < data.size() - 1) && exp_it!=experiment.end(); ++i )
    { // .. and disturb values by (an already smoothed) distortion diced in RTSimulation
			data[i] *= (DoubleReal) exp_it->getMetaValue("distortion");
    }
  
      
  }

  void RawMSSignalSimulation::addShotNoise_(MSSimExperiment & experiment)
  {
    // we model the amount of (background) noise as Poisson process
    // i.e. the number of noise data points per unit m/z interval follows a Poisson
    // distribution. Noise intensity is assumed to be Gaussian-distributed.

    DoubleReal rate    = param_.getValue("noise_rate");
    DoubleReal it_mean = param_.getValue("noise_int_mean");

    const UInt num_intervals = 100;
    SimCoordinateType interval_size = ( maximal_mz_measurement_limit_ -  minimal_mz_measurement_limit_) / num_intervals;
    SimPointType point;

    std::cout << "Adding shot noise to spectra...." << std::endl;
    std::cout << "Interval size: "  << interval_size << " poisson rate: " << rate << std::endl;

    // TODO: switch to iterator ??
    for (Size i=0 ;i < experiment.size() ; ++i)
    {

      for (Size j=0;j<num_intervals;++j)
      {
        UInt counts = gsl_ran_poisson ( rnd_gen_, rate);
        SimCoordinateType mz_lw = j * interval_size + minimal_mz_measurement_limit_;
        SimCoordinateType mz_up = (j+1) * interval_size + minimal_mz_measurement_limit_;

        for (UInt c=0; c<counts;++c)
        {
          SimCoordinateType mz  = gsl_ran_flat(rnd_gen_, mz_lw, mz_up );
          SimCoordinateType it = gsl_ran_exponential(rnd_gen_,it_mean);
          point.setIntensity(it);
          point.setMZ(mz);
          experiment[i].push_back(point);
        }

      }
    } // end of each scan

    experiment.updateRanges();

  }

  void RawMSSignalSimulation::addBaseLine_(MSSimExperiment & experiment)
  {
    DoubleReal scale = param_.getValue("baseline_scaling");
    if (scale == 0.0) return;

    for ( Size i = 0; i < experiment.size() ; ++i )
    {
      for ( Size j = 0 ; j < experiment[i].size() ; ++j )
      {
        SimCoordinateType x = (experiment[i][j].getMZ() - minimal_mz_measurement_limit_);
        // TODO: what is this?! Baseline up to 1000Th and then what?
        // We should consider using an exp() for MALDI and nothing(?) for ESI?
        if (x >= 1000.0) continue; // speed-up

        DoubleReal b = gsl_ran_exponential_pdf(x,125.0);
        b *= scale;
        experiment[i][j].setIntensity( experiment[i][j].getIntensity() + b );
      }
    }
  }

  void RawMSSignalSimulation::compressSignals_(MSSimExperiment & experiment)
  {
		// assume we changed every scan a priori
    //changed_scans_.resize(experiment.size(), true);
  
    Size count = 0;
    do
    { // TODO here we can improve on speed a lot by better compression! check if running times are significant
      // TODO: either use in compressSignalsRun_ or delete the member
	    //changed_scans_[ i ]  = false;
      count = compressSignalsRun_(experiment);
    } while (count != 0);
  }

  Size RawMSSignalSimulation::compressSignalsRun_(MSSimExperiment & experiment)
  {
    SimCoordinateType diff_mz = 0.0;
    SimPointType p;

    Size count = 0;
    bool change = false;

    for( Size i = 0 ; i < experiment.size() ; ++i )
    {
      experiment[i].sortByPosition();

      // copy Spectrum and remove Peaks ..
      MSSimExperiment::SpectrumType cont = experiment[i];
      cont.clear();

      for ( Size j = 0 ; j < (experiment[i].size() - 1) ; ++j )
      {
        diff_mz = fabs(experiment[i][ (j+1) ].getMZ() - experiment[i][j].getMZ());

        if (diff_mz < mz_sampling_rate_)
        {
          change = true;
          // sum intensities
          SimCoordinateType it1 = experiment[i][ (j+1) ].getIntensity();
          SimCoordinateType it2 = experiment[i][ (j) ].getIntensity();
          SimCoordinateType it =  it1 + it2;
          p.setIntensity( it );

          // keep m/z of point with higher intensity
          SimCoordinateType mz1 = experiment[i][ (j+1) ].getMZ();
          SimCoordinateType mz2 = experiment[i][ (j) ].getMZ();
          SimCoordinateType mz =  it1 > it2 ? mz1 : mz2;
          p.setMZ( mz );
          cont.push_back(p);

          ++j;
          ++count;
        }
        else
        {
          cont.push_back( experiment[i][j] );
        }
      }
      // don't forget the last one
      if (!change) cont.push_back( experiment[i][ (experiment[i].size() - 1) ] );
      experiment[i] = cont;
    }

#ifdef DEBUG_SIM
    cout << "Done " << endl;
    cout << "Count: " << count << endl;
#endif

    return count;

  }
}
