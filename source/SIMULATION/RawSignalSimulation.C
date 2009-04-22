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
// $Authors: Stephan Aiche Chris Bielow$
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

  void RawSignalSimulation::generateRawSignals(FeatureMap< > &, MSExperiment< Peak1D > &)
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
  }

  void RawSignalSimulation::addMSSignal(Feature & activeFeature, MSExperiment< Peak1D > & expirement)
  {
    ProductModel<2> pm;
    Param p1;

    SimIntensityType scale = activeFeature.getIntensity() * 1500; // was: 3000 TODO: ???? why 1500
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
    // TODO: store all simulated features ..

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
                                    MSExperiment<Peak1D> & expirement, Feature & activeFeature)
  {
    // start and end points of the sampling are entirely arbitrary
    // and should be modified at some point

    // (cg) commented this out since it cuts off fronted elution profiles!!
    // (ost) Why should this happen ?
    if (rt_start <=0) rt_start = 0;

    MSExperiment<Peak1D>::iterator exp_iter = expirement.RTBegin(rt_start);

    SimIntensityType intensity_sum = 0.0;
    vector< DPosition<2> > points;

#ifdef DEBUG_SIM
    std::cout << "Sampling at: " << mz_start << " " << mz_end << " ";
    std::cout << rt_start << " " << rt_end << std::endl;
#endif

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
            //cout << "start_scan: " << start_scan << endl;
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
          //cout << "Sampling intensity: " << point.getIntensity() << endl;

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
      double symmetry = gsl_ran_flat (rnd_gen_, -100.0, 100.0);

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
      for ( UInt i = 1; i < data.size() - 1; ++i )
      {
        data[i] *= exp(gsl_ran_flat (rnd_gen_, -distortion_, +distortion_));
      }
      // moving average filter (width 3), implemented very inefficiently (guess why!)
      if ( distortion_ != 0.0 ) // otherwise we want perfect EMG shape!
      {
        ElutionModel::ContainerType tmp;
        tmp.resize(data.size());
        for ( UInt i = 1; i < data.size() - 1; ++i )
        {
          tmp[i] = ( data[i-1] + data[i] + data[i+1] ) / 3.0;
        }
        for ( UInt i = 1; i < data.size() - 1; ++i )
        {
          data[i] = tmp[i];
        }

        const int num_rounds = 10;
        for ( int rounds = 0; rounds < num_rounds; ++rounds)
        {
          data.swap(tmp);
          for ( UInt i = 1; i < data.size() - 1; ++i )
          {
            tmp[i] = ( data[i-1] + data[i] + data[i+1] ) / 3.0;
          }
          for ( UInt i = 1; i < data.size() - 1; ++i )
          {
            data[i] = tmp[i];
          }
        }

      }
      pm.setModel(0,elutionmodel);
  }
}
