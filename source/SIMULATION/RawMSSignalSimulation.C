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

#include <OpenMS/SIMULATION/IsotopeModelGeneral.h>

#include <vector>
using std::vector;

namespace OpenMS {

  /**
   * TODO: review signal compression code
   * TODO: review baseline and noise code
   */
  RawMSSignalSimulation::RawMSSignalSimulation(const gsl_rng * random_generator)
  : DefaultParamHandler("RawSignalSimulation"), mz_sampling_rate_(), mz_error_mean_(), mz_error_stddev_(),
  intensity_error_mean_(), intensity_error_stddev_(), peak_std_(), rnd_gen_(random_generator)
  {
    setDefaultParams_();
    updateMembers_();
  }


  RawMSSignalSimulation::RawMSSignalSimulation()
    : DefaultParamHandler("RawSignalSimulation"), mz_sampling_rate_(), mz_error_mean_(), mz_error_stddev_(),
    intensity_error_mean_(), intensity_error_stddev_(), peak_std_()
  {
    setDefaultParams_();
    updateMembers_();
  }

  RawMSSignalSimulation::RawMSSignalSimulation(const RawMSSignalSimulation& source)
    : DefaultParamHandler(source), mz_sampling_rate_(source.mz_sampling_rate_), mz_error_mean_(source.mz_error_mean_), mz_error_stddev_(source.mz_error_stddev_),
    intensity_error_mean_(source.intensity_error_mean_), intensity_error_stddev_(source.intensity_error_stddev_), peak_std_(source.peak_std_)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
  }

  RawMSSignalSimulation& RawMSSignalSimulation::operator = (const RawMSSignalSimulation& source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;

    mz_error_mean_ = source.mz_error_mean_;
    mz_error_stddev_ = source.mz_error_stddev_;
    mz_sampling_rate_ = source.mz_sampling_rate_;

    intensity_error_mean_ = source.intensity_error_mean_;
    intensity_error_stddev_ = source.intensity_error_stddev_;

    peak_std_ = source.peak_std_;

    updateMembers_();
    return *this;
  }

  RawMSSignalSimulation::~RawMSSignalSimulation()
  {}

  void RawMSSignalSimulation::setDefaultParams_()
  {
		defaults_.setValue("enabled","true","Enable RAW signal simulation?");
		defaults_.setValidStrings("enabled", StringList::create("true,false"));
		
    // noise params
    // m/z error
    defaults_.setValue("mz:error_mean",0.0,"Average systematic m/z error (Da)");
    defaults_.setValue("mz:error_stddev",0.0,"Standard deviation for m/z errors. Set to 0 to disable simulation of m/z errors.");

    // mz sampling rate
    defaults_.setValue("mz:sampling_rate",0.12,"MS hardware resolution (e.g. bin size in m/z).");

    // intensity error
    defaults_.setValue("int:error_mean",0,"Average systematic intensity error.");
    defaults_.setValue("int:error_stddev",0.0,"Standard deviation for peak intensities (relative to peak height). Set to 0 to disable intensity errors.");

    // peak and instrument parameter
    defaults_.setValue("peak_fwhm",0.5,"FWHM (full width at half maximum) of simulated peaks (Da).");

    // shot noise
    defaults_.setValue("noise:rate",0.0,"Poisson rate of shot noise. Set to 0 to disable simulation of shot noise.");
    defaults_.setMinFloat("noise:rate",0.0);
    defaults_.setValue("noise:int-mean",50.0,"Shot noise intensity mean.");

    // baseline
    defaults_.setValue("baseline:scaling",0.0,"Scale of baseline. Set to 0 to disable simulation of baseline.");
    defaults_.setMinFloat("baseline:scaling",0.0);
    defaults_.setValue("baseline:shape",0.5, "The baseline is modeled by an exponential probability density function (pdf) with f(x) = shape*e^(- shape*x)");
    defaults_.setMinFloat("baseline:shape",0.0);

    defaultsToParam_();
  }

  void RawMSSignalSimulation::updateMembers_()
  {
    double tmp    =  param_.getValue("peak_fwhm");
    peak_std_     = (tmp / 2.355);			// Approximation for Gaussian-shaped signals
    mz_sampling_rate_ = param_.getValue("mz:sampling_rate");

    mz_error_mean_    = param_.getValue("mz:error_mean");
    mz_error_stddev_  = param_.getValue("mz:error_stddev");

    intensity_error_mean_   = param_.getValue("int:error_mean");
    intensity_error_stddev_ = param_.getValue("int:error_stddev");

  }

  void RawMSSignalSimulation::generateRawSignals(FeatureMapSim & features, MSSimExperiment & experiment)
  {
    if (param_.getValue("enabled") == "false")
    {
			return;
		}
		// retrieve mz boundary parameters from experiment:
		SimCoordinateType minimal_mz_measurement_limit = experiment[0].getInstrumentSettings().getScanWindows()[0].begin;
		SimCoordinateType maximal_mz_measurement_limit = experiment[0].getInstrumentSettings().getScanWindows()[0].end;
		
    if(experiment.size() == 1)
    {
      for(FeatureMap< >::iterator feature_it = features.begin();
          feature_it != features.end();
          ++feature_it)
      {
        add1DSignal_(*feature_it,experiment);
      }
      addShotNoise_(experiment, minimal_mz_measurement_limit, maximal_mz_measurement_limit);
    }
    else
    {
      for(Size idx=0; idx<features.size(); ++idx)
      {
        add2DSignal_(features[idx], experiment);
        if (idx % (features.size()/10+1) == 0) std::cout << idx << " of " << features.size() << " MS1 features generated...\n";
      }
      addShotNoise_(experiment, minimal_mz_measurement_limit, maximal_mz_measurement_limit);
    }
    // we should trigger this based on ionization type
    addBaseLine_(experiment, minimal_mz_measurement_limit);
    compressSignals_(experiment);
  }

  void RawMSSignalSimulation::add1DSignal_(Feature & active_feature, MSSimExperiment & experiment)
  {
    Param p1;

    // was: 3000 TODO: ???? why 1500
    SimIntensityType scale = active_feature.getIntensity() * 150;

    SimChargeType q = active_feature.getCharge();
    EmpiricalFormula ef = active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula();
    ef += active_feature.getMetaValue("charge_adducts"); // adducts
    ef -= String("H")+String(q);ef.setCharge(q);				 // effectively substract q electrons
    p1.setValue("statistics:mean", ef.getAverageWeight() / q);		

    p1.setValue("statistics:mean", active_feature.getMZ() );
    p1.setValue("interpolation_step", 0.001);
    p1.setValue("isotope:stdev", peak_std_);
    p1.setValue("intensity_scaling", scale);
    p1.setValue("charge", q);

    IsotopeModelGeneral isomodel;
    isomodel.setSamples(ef);
    isomodel.setParameters(p1);

    SimCoordinateType mz_start = isomodel.getInterpolation().supportMin();
    SimCoordinateType mz_end = isomodel.getInterpolation().supportMax();
        
    samplePeptideModel1D_(isomodel, mz_start, mz_end, experiment, active_feature);
  }

  void RawMSSignalSimulation::add2DSignal_(Feature & active_feature, MSSimExperiment & experiment)
  {
    // was: 3000 TODO: ???? why 1500
    // TODO: we need to improve this
    SimIntensityType scale = active_feature.getIntensity() * 1500;

    Param p1;
    SimChargeType q = active_feature.getCharge();
    EmpiricalFormula ef = active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula();
    //std::cout << "current feature: " << active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().toString() << " with scale " << scale << std::endl;
    ef += active_feature.getMetaValue("charge_adducts"); // adducts
    ef -= String("H")+String(q);ef.setCharge(q);				 // effectively substract q electrons
    p1.setValue("statistics:mean", ef.getAverageWeight() / q);
    p1.setValue("interpolation_step", 0.001);
    p1.setValue("isotope:stdev", peak_std_);
    p1.setValue("charge", q);

    IsotopeModelGeneral* isomodel = new IsotopeModelGeneral();
    isomodel->setSamples(ef); // this already includes adducts
    isomodel->setParameters(p1);

		if (experiment.size()<2)
		{
			throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, experiment.size());
		}
		DoubleReal rt_sampling_rate = experiment[1].getRT() - experiment[0].getRT();
		EGHModel* elutionmodel = new EGHModel();
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

    // add peptide to global MS map
    samplePeptideModel2D_(pm, mz_start, mz_end, rt_start, rt_end, experiment, active_feature);
  }


  void RawMSSignalSimulation::samplePeptideModel1D_(const IsotopeModelGeneral & pm,
																										const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
																										MSSimExperiment & experiment, Feature & active_feature)
  {
    SimIntensityType intensity_sum = 0.0;

    //std::cout << "Sampling at [mz] " << mz_start << ":" << mz_end << std::endl;

    /// TODO: think of better error checking

    SimPointType point;

    for (SimCoordinateType mz = mz_start; mz < mz_end; mz += mz_sampling_rate_)
    {
      point.setMZ(mz);
      point.setIntensity( pm.getIntensity( DPosition<1>( mz ) ) );

      if ( point.getIntensity() > 10.0)
      {
        // add m/z and intensity error (both Gaussian distributed)
        double it_err  = gsl_ran_gaussian(rnd_gen_, (point.getIntensity() * intensity_error_stddev_ ) ) + intensity_error_mean_ ;
        point.setIntensity( std::max(0., point.getIntensity( ) + it_err) );
				
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

    MSSimExperiment::iterator exp_iter = experiment.RTBegin(rt_start);
    if(exp_iter == experiment.end() )
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 0);
    }

    //std::cout << "Sampling at [RT] " << rt_start << ":" << rt_end << " [mz] " << mz_start << ":" << mz_end << std::endl;

    SimIntensityType intensity_sum = 0.0;
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

        point.setMZ(mz);
        point.setIntensity( pm.getIntensity( DPosition<2>( rt, mz) ) );

        if ( point.getIntensity() > 10.0)
        {
          // add m/z and intensity error (both Gaussian distributed)
          double it_err  = gsl_ran_gaussian(rnd_gen_, (point.getIntensity() * intensity_error_stddev_ ) ) + intensity_error_mean_ ;
          point.setIntensity( std::max(0., point.getIntensity( ) + it_err) );

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

    active_feature.setQuality(0,start_scan);
    active_feature.setQuality(1,end_scan);

    active_feature.setIntensity(intensity_sum);
    // store convex hull
    active_feature.getConvexHulls().clear();
    // adding ALL points of the feature, the ConvexHull2D assignment operator will take care of the rest
    active_feature.getConvexHulls().push_back(points);

  }

  void RawMSSignalSimulation::chooseElutionProfile_(EGHModel*& elutionmodel, const Feature& feature, const double scale, const DoubleReal rt_sampling_rate, const MSSimExperiment & experiment)
    {
      SimCoordinateType f_rt = feature.getRT();

      Param p;
      // WARNING: step used to be 'rt_sampling_rate / 3.0', but distortion is not part of RT sim, and thus only
      //          modeled 1:1
      p.setValue("interpolation_step", rt_sampling_rate );
      p.setValue("statistics:variance", 1.0);
      p.setValue("statistics:mean", f_rt);

      p.setValue("egh:height", scale);
      p.setValue("egh:A", 50.0);
      p.setValue("egh:B", 90.0);

      elutionmodel->setParameters(p); // does the calculation

      //----------------------------------------------------------------------

      // Hack away constness :-P   We know what we want.
      EmgModel::ContainerType &data = const_cast<EGHModel::ContainerType &>(elutionmodel->getInterpolation().getData());

      SimCoordinateType rt_em_start = elutionmodel->getInterpolation().supportMin();

      // find scan in experiment at which our elution starts
      MSSimExperiment::ConstIterator exp_it = experiment.RTBegin(rt_em_start);
      for ( Size i = 1; (i < data.size() - 1) && exp_it!=experiment.end(); ++i, ++exp_it )
      { // .. and disturb values by (an already smoothed) distortion diced in RTSimulation
        data[i] *= (DoubleReal) exp_it->getMetaValue("distortion");
      }


    }


  void RawMSSignalSimulation::chooseElutionProfile_(EmgModel*& elutionmodel, const Feature& feature, const double scale, const DoubleReal rt_sampling_rate, const MSSimExperiment & experiment)
  {
    SimCoordinateType f_rt = feature.getRT();
    DoubleReal f_symmetry = (DoubleReal) feature.getMetaValue("rt_symmetry");
    DoubleReal f_width = (DoubleReal) feature.getMetaValue("rt_width");

    DoubleReal f_bb_min = (DoubleReal) feature.getMetaValue("rt_bb_min");
    DoubleReal f_bb_max = (DoubleReal) feature.getMetaValue("rt_bb_max");

    // it doesn't matter what bounding box I set here, it always cuts of the fronting !!! :-(
    Param p;
    p.setValue("bounding_box:min", f_bb_min );
    p.setValue("bounding_box:max", f_bb_max );
    // WARNING: step used to be 'rt_sampling_rate / 3.0', but distortion is not part of RT sim, and thus only
    //          modeled 1:1
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
    EmgModel::ContainerType &data = const_cast<EmgModel::ContainerType &>(elutionmodel->getInterpolation().getData());

		SimCoordinateType rt_em_start = elutionmodel->getInterpolation().supportMin();
	
		// find scan in experiment at which our elution starts
		MSSimExperiment::ConstIterator exp_it = experiment.RTBegin(rt_em_start);
    for ( Size i = 1; (i < data.size() - 1) && exp_it!=experiment.end(); ++i, ++exp_it )
    { // .. and disturb values by (an already smoothed) distortion diced in RTSimulation
			data[i] *= (DoubleReal) exp_it->getMetaValue("distortion");
    }
  
      
  }

  void RawMSSignalSimulation::addShotNoise_(MSSimExperiment & experiment, SimCoordinateType minimal_mz_measurement_limit, SimCoordinateType maximal_mz_measurement_limit)
  {
    // we model the amount of (background) noise as Poisson process
    // i.e. the number of noise data points per unit m/z interval follows a Poisson
    // distribution. Noise intensity is assumed to be Gaussian-distributed.

    DoubleReal rate    = param_.getValue("noise:rate");
    if (rate == 0.0) return;
    DoubleReal it_mean = param_.getValue("noise:int-mean");

    const UInt num_intervals = 100;
    SimCoordinateType interval_size = ( maximal_mz_measurement_limit -  minimal_mz_measurement_limit) / num_intervals;
    SimPointType point;

    std::cout << "Adding shot noise to spectra...." << std::endl;
    std::cout << "Interval size: "  << interval_size << " poisson rate: " << rate << std::endl;

    // TODO: switch to iterator ??
    for (Size i=0 ;i < experiment.size() ; ++i)
    {

      for (Size j=0;j<num_intervals;++j)
      {
        UInt counts = gsl_ran_poisson ( rnd_gen_, rate);
        SimCoordinateType mz_lw = j * interval_size + minimal_mz_measurement_limit;
        SimCoordinateType mz_up = (j+1) * interval_size + minimal_mz_measurement_limit;

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

  void RawMSSignalSimulation::addBaseLine_(MSSimExperiment & experiment, SimCoordinateType minimal_mz_measurement_limit)
  {
    DoubleReal scale = param_.getValue("baseline:scaling");
    DoubleReal shape = param_.getValue("baseline:shape");

    if (scale == 0.0) return;

    for ( Size i = 0; i < experiment.size() ; ++i )
    {
      for ( Size j = 0 ; j < experiment[i].size() ; ++j )
      {
        SimCoordinateType x = (experiment[i][j].getMZ() - minimal_mz_measurement_limit);        

        //if (x >= 1000.0) continue; // speed-up TODO: revise this ..

        DoubleReal b = gsl_ran_exponential_pdf(x, shape);
        b *= scale;
        experiment[i][j].setIntensity( experiment[i][j].getIntensity() + b );
      }
    }
  }

  void RawMSSignalSimulation::compressSignals_(MSSimExperiment & experiment)
  {
		// assume we changed every scan a priori

    Size count = 0;
    do
    { // TODO here we can improve on speed a lot by better compression! check if running times are significant
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

			if (experiment[i].size() >=2)
			{
				// copy Spectrum and remove Peaks ..
				MSSimExperiment::SpectrumType cont = experiment[i];
				cont.clear(false);
				
				for ( Size j = 0 ; j < experiment[i].size() -1 ; ++j )
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
						change = false;
						cont.push_back( experiment[i][j] );
					}
				}
				// don't forget the last one
				if (!change) cont.push_back( experiment[i][ (experiment[i].size() - 1) ] );
	      experiment[i] = cont;
			}
    }

#ifdef DEBUG_SIM
    cout << "Done " << endl;
    cout << "Count: " << count << endl;
#endif

    return count;

  }
}
