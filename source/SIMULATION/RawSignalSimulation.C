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

#include <OpenMS/SIMULATION/RawSignalSimulation.h>
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
    * add shot noise
    * add base line
  */
  RawSignalSimulation::RawSignalSimulation(const gsl_rng * random_generator)
  : DefaultParamHandler("RawSignalSimulation"), rnd_gen_(random_generator)
  {
    setDefaultParams_();
    updateMembers_();
  }


  RawSignalSimulation::RawSignalSimulation()
    : DefaultParamHandler("RawSignalSimulation")
  {
    setDefaultParams_();
    updateMembers_();
  }

  RawSignalSimulation::RawSignalSimulation(const RawSignalSimulation& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
  }

  RawSignalSimulation& RawSignalSimulation::operator = (const RawSignalSimulation& source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
    return *this;
  }

  RawSignalSimulation::~RawSignalSimulation()
  {}

  void RawSignalSimulation::setDefaultParams_()
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

    // rt parameters
    defaults_.setValue("rt:sampling_rate",2.0,"Time interval (in seconds) between consecutive scans");

    // peak and instrument parameter
    defaults_.setValue("peak_fwhm",0.5,"FWHM (full width at half maximum) of simulated peaks (Da).");

    // shot noise
    defaults_.setValue("noise_rate",0,"Poisson rate of shot noise. Set to 0 to disable simulation of shot noise.");
    defaults_.setValue("noise_int_mean",50.0,"Shot noise intensity mean.");
    
    // baseline
    defaults_.setValue("baseline_scaling",0.0,"Scale of baseline (zero disables baseline simulation)");
    
    // column conditions
    defaults_.setValue("column_condition:preset","medium","LC condition (none|good|medium|poor) if set to none the explicit values will be used.");
    StringList valid_presets = StringList::create("none,good,medium,poor");
    defaults_.setValidStrings("column_condition:preset", valid_presets);
    
    defaults_.setValue("column_condition:distortion", 1.0, "LC distortion (ignored if preset is not set to none)");
    defaults_.setValue("column_condition:symetry_up", -60.0, "LC symetry up (ignored if preset is not set to none)");
    defaults_.setValue("column_condition:symetry_down", +60.0, "LC symetry down (ignored if preset is not set to none)");
    
    defaultsToParam_();
  }

  void RawSignalSimulation::updateMembers_()
  {
    double tmp    =  param_.getValue("peak_fwhm");
    peak_std_     = (tmp / 2.355);			// Approximation for Gaussian-shaped signals
    mz_sampling_rate_ = param_.getValue("mz:sampling_rate");

    maximal_mz_measurement_limit_ = param_.getValue("mz:upper_measurement_limit");
    minimal_mz_measurement_limit_ = param_.getValue("mz:lower_measurement_limit");

    mz_error_mean_    = param_.getValue("mz:error_mean");
    mz_error_stddev_  = param_.getValue("mz:error_stddev");

    intensity_error_mean_   = param_.getValue("int:error_mean");
    intensity_error_stddev_ = param_.getValue("int:error_stddev");

    rt_sampling_rate_ = param_.getValue("rt:sampling_rate");
    
    String column_preset = param_.getValue("column_condition:preset");
    if (column_preset == "poor")
    {
      distortion_    = 2.0;
      symmetry_down_ = -100;
      symmetry_up_   = +100;
    }
    else if (column_preset == "medium")
    {
      distortion_    = 1.0;
      symmetry_down_ = -60;
      symmetry_up_   = +60;
    }
    else if (column_preset == "good")
    {
      distortion_    = 0.0;
      symmetry_down_ = -15;
      symmetry_up_   = +15;
    }
    else 	// default is "none" so get user set parameters
    {
      distortion_    = param_.getValue("column_condition:distortion");
      symmetry_up_   = param_.getValue("column_condition:symetry_up");
      symmetry_down_ = param_.getValue("column_condition:symetry_down");
    }
    
  }

  void RawSignalSimulation::generateRawSignals(FeatureMapSim & features, MSSimExperiment & experiment)
  {
    changed_scans_.resize(experiment.size());
    for(FeatureMap< >::iterator feature_it = features.begin();
        feature_it != features.end();
        ++feature_it)
    {
      addMSSignal(*feature_it, experiment);
    }
    addShotNoise_(experiment);
    compressSignals_(experiment);
    addBaseLine_(experiment);
  }
  
  void RawSignalSimulation::addMSSignal(Feature & activeFeature, MSSimExperiment & expirement)
  {
    ProductModel<2> pm;
    Param p1;

    // was: 3000 TODO: ???? why 1500
    SimIntensityType scale = activeFeature.getIntensity() * 1500;
    mean_scaling_ += scale;
    ++ion_count_;

    SimChargeType charge = activeFeature.getCharge();

    //TODO use H+ weight instead of 1*c ?
    EmpiricalFormula feature_ef = activeFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula();

    SimCoordinateType mz = ( (feature_ef.getMonoWeight() + charge) / charge) ;

    // TODO: this should not be necessary, check ...
    /*
    allow_overlaps_ = (unsigned int) param_.getValue("allow_overlaps");
    if ( (allow_overlaps_ == 0) && checkForOverlaps_( DPosition<2>(mz,rt) ) )
    {
      return;
    }
    */

    // don't show ions with m/z higher than the MS detection limit
    // TODO: integrate map ranges
    if (mz > maximal_mz_measurement_limit_ || mz < minimal_mz_measurement_limit_)
    {
      // ignore the current feature
      return;
    }
    activeFeature.setMZ(mz);

    //TODO use H+ weight instead of 1*c ?
    p1.setValue("statistics:mean",((feature_ef.getAverageWeight()+charge)/charge));
    p1.setValue("interpolation_step",0.001);
    p1.setValue("isotope:stdev",peak_std_);
    p1.setValue("charge",(int) charge);

    IsotopeModelGeneral* isomodel = new IsotopeModelGeneral();
    isomodel->setSamples(feature_ef);
    isomodel->setParameters(p1);

    chooseElutionProfile_(pm,activeFeature.getRT(),scale);
    pm.setModel(1,isomodel);
    pm.setScale(scale);

    // add peptide to global MS map
    // TODO: use flexible boundaries for rt/mz depending on abundance/charge?
    // TODO: remove this current_feature_ dependency
    // TODO: store all simulated features in a separate FeatureMap ..

    samplePeptideModel_(pm, (mz - 2.5),(mz + 5.0), (activeFeature.getRT() - 160.0),(activeFeature.getRT() + 280.0), expirement, activeFeature);
    /*
    if (current_feature_.getConvexHulls().begin()->getPoints().size() > 0)
    {
      features_.push_back(current_feature_);
    }
    */
  }

  void RawSignalSimulation::samplePeptideModel_(const ProductModel<2> & pm,
                                    const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
                                    SimCoordinateType rt_start, SimCoordinateType rt_end,
                                    MSSimExperiment & expirement, Feature & activeFeature)
  {
    // start and end points of the sampling are entirely arbitrary
    // and should be modified at some point

    // (cg) commented this out since it cuts off fronted elution profiles!!
    // (ost) Why should this happen ?
    if (rt_start <=0) rt_start = 0;

    MSExperiment<Peak1D>::iterator exp_iter = expirement.RTBegin(rt_start);

    SimIntensityType intensity_sum = 0.0;
    vector< DPosition<2> > points;

//#ifdef DEBUG_SIM
    std::cout << "Sampling at: " << mz_start << " " << mz_end << " ";
    std::cout << rt_start << " " << rt_end << std::endl;
//#endif

    /// TODO: think of better error checking
    if(exp_iter == expirement.end() )
    {
      std::cout << "error ! " << std::endl; // ;-) should not happen
      return;
    }

    SimPointType point;

    Int start_scan = -5;
    Int end_scan  = -5;

    //UInt it = 0;
    //UInt pit = 0;

    for (SimCoordinateType rt = rt_start; rt < rt_end && exp_iter != expirement.end(); rt += rt_sampling_rate_, ++exp_iter)
    {
      for (SimCoordinateType mz = mz_start; mz < mz_end; mz += mz_sampling_rate_)
      {
        //++it;

        point.setMZ(mz);
        point.setIntensity( pm.getIntensity( DPosition<2>( rt, mz) ) );

        if ( point.getIntensity() > 10.0)
        {
          if (start_scan == -5)
          {
            start_scan = exp_iter - expirement.begin();
            //std::cout << "start_scan: " << start_scan << std::endl;
          }

          if (! changed_scans_.at( exp_iter - expirement.begin() ) )
          {
            changed_scans_.at( exp_iter - expirement.begin() )  = true;
          }
          //++pit;

          // add m/z and itensity error (both Gaussian distributed)
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

          //update last scan
          end_scan = exp_iter - expirement.begin();
        }
      }
    }

    //cout << "End of sampling: " << it << " vs " << pit << endl;

    // do not set this here, because it might include low intensity points and is inconsistent with the convex hull
    //end_scan = exp_iter - exp_.begin();


    //cout << "end_scan: " << end_scan << endl;

    // This is a clear misuse of the Feature data structure
    // but at this point, we couldn't care less ;-)
    // TODO: this was currentFeature -> reimplement
    activeFeature.setQuality(0,start_scan);
    activeFeature.setQuality(1,end_scan);

    activeFeature.setIntensity(intensity_sum);
    // store convex hull
    activeFeature.getConvexHulls().clear();
    activeFeature.getConvexHulls().resize( activeFeature.getConvexHulls().size()+1);
    activeFeature.getConvexHulls()[ activeFeature.getConvexHulls().size()-1] = points;

//#ifdef DEBUG_SIM
//    current_feature_.setModelDescription( ModelDescription<2>( &pm ) );
//#endif
  }

  void RawSignalSimulation::chooseElutionProfile_(ProductModel<2>& pm, const SimCoordinateType rt,const double scale)
  {
      Param p;
      double symmetry = gsl_ran_flat (rnd_gen_, symmetry_down_, symmetry_up_);

      double width = gsl_ran_flat (rnd_gen_, 5, 15);

      // Exponentially modified Gaussian
      const DoubleReal decay_stretch = 5.0;
      // it doesn't matter what bounding box I set here, it always cuts of the fronting !!! :-(
      p.setValue("bounding_box:min",rt - decay_stretch*(width+abs(symmetry)) );
      p.setValue("bounding_box:max",rt + decay_stretch*(width+abs(symmetry)) );
      p.setValue("interpolation_step", rt_sampling_rate_ / 3.0);
      p.setValue("statistics:variance",1.0);
      p.setValue("statistics:mean",rt );
      p.setValue("emg:height",scale);
      p.setValue("emg:width",width);
      p.setValue("emg:symmetry",symmetry);
      p.setValue("emg:retention",rt);
      ElutionModel* elutionmodel = new ElutionModel();
      elutionmodel->setParameters(p);
      elutionmodel->setScalingFactor(scale);

      //----------------------------------------------------------------------

      // So far we had randomness in the parameters of the EMG.  Now let us add
      // some random ups and downs within the elution profile.

      // Hack away constness :-P   We know what we want.
      ElutionModel::ContainerType &data = const_cast<ElutionModel::ContainerType &>(elutionmodel->getInterpolation().getData());

      // Member distortion_ is the logarithm of the maximal factor which can be applied to
      // each point. But note that elution profile is smoothed again afterward.

      // first and last entry shall not be changed
      for ( Size i = 1; i < data.size() - 1; ++i )
      {
        data[i] *= exp(gsl_ran_flat (rnd_gen_, -distortion_, +distortion_));
      }
      // moving average filter (width 3), implemented very inefficiently (guess why!)
      if ( distortion_ != 0.0 ) // otherwise we want perfect EMG shape!
      {
        ElutionModel::ContainerType tmp;
        tmp.resize(data.size());
        for ( Size i = 1; i < data.size() - 1; ++i )
        {
          tmp[i] = ( data[i-1] + data[i] + data[i+1] ) / 3.0;
        }
        for ( Size i = 1; i < data.size() - 1; ++i )
        {
          data[i] = tmp[i];
        }

        const int num_rounds = 10;
        for ( int rounds = 0; rounds < num_rounds; ++rounds)
        {
          data.swap(tmp);
          for ( Size i = 1; i < data.size() - 1; ++i )
          {
            tmp[i] = ( data[i-1] + data[i] + data[i+1] ) / 3.0;
          }
          for ( Size i = 1; i < data.size() - 1; ++i )
          {
            data[i] = tmp[i];
          }
        }

      }
      pm.setModel(0,elutionmodel);
  }
  
  SimCoordinateType RawSignalSimulation::getRTSamplingRate() const
  {
    return rt_sampling_rate_;
  }
  
  
  void RawSignalSimulation::addShotNoise_(MSSimExperiment & experiment)
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
  
  void RawSignalSimulation::addBaseLine_(MSSimExperiment & experiment)
  {
    DoubleReal scale = param_.getValue("baseline_scaling");
    if (scale == 0.0) return;
    
    for ( Size i = 0; i < experiment.size() ; ++i )
    {
      for ( Size j = 0 ; j < experiment[i].size() ; ++j )
      {
        SimCoordinateType x = (experiment[i][j].getMZ() - minimal_mz_measurement_limit_);
        if (x >= 1000.0) continue; // speed-up
        
        DoubleReal b = gsl_ran_exponential_pdf(x,125.0);
        b *= scale;
        experiment[i][j].setIntensity( experiment[i][j].getIntensity() + b );
      }
    }
  }
  
  void RawSignalSimulation::compressSignals_(MSSimExperiment & experiment)
  {
    Size count = 0;
    do
    {
      count = compressSignalsRun_(experiment);
    } while (count != 0);
  }
  
  Size RawSignalSimulation::compressSignalsRun_(MSSimExperiment & experiment)
  {
    SimCoordinateType diff_mz = 0.0;
    SimPointType p;
    
    Size count = 0;
    bool change = false;
    
    for( Size i = 0 ; i < experiment.size() ; ++i )
    {
      experiment[i].sortByPosition();
      if (experiment[i].size() <= 2) continue;
      
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
