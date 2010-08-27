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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <OpenMS/SIMULATION/IsotopeModelGeneral.h>

#include <vector>
using std::vector;

namespace OpenMS {

  /**
   * TODO: review baseline and noise code
   */
  RawMSSignalSimulation::RawMSSignalSimulation(const SimRandomNumberGenerator& rng)
  : DefaultParamHandler("RawSignalSimulation"), mz_sampling_rate_(), mz_error_mean_(), mz_error_stddev_(),
  intensity_scale_(), intensity_scale_stddev_(), peak_std_(), rnd_gen_(&rng)
  {
    setDefaultParams_();
    updateMembers_();
  }


  RawMSSignalSimulation::RawMSSignalSimulation()
    : DefaultParamHandler("RawSignalSimulation"), mz_sampling_rate_(), mz_error_mean_(), mz_error_stddev_(),
    intensity_scale_(), intensity_scale_stddev_(), peak_std_()
  {
    setDefaultParams_();
    updateMembers_();
  }

  RawMSSignalSimulation::RawMSSignalSimulation(const RawMSSignalSimulation& source)
    : DefaultParamHandler(source), mz_sampling_rate_(source.mz_sampling_rate_), mz_error_mean_(source.mz_error_mean_), mz_error_stddev_(source.mz_error_stddev_),
    intensity_scale_(source.intensity_scale_), intensity_scale_stddev_(source.intensity_scale_stddev_), peak_std_(source.peak_std_)
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

    intensity_scale_ = source.intensity_scale_;
    intensity_scale_stddev_ = source.intensity_scale_stddev_;

    peak_std_ = source.peak_std_;

    updateMembers_();
    return *this;
  }

  RawMSSignalSimulation::~RawMSSignalSimulation()
  {
  }

  void RawMSSignalSimulation::setDefaultParams_()
  {
		defaults_.setValue("enabled","true","Enable RAW signal simulation? (select 'false' if you only need feature-maps)");
		defaults_.setValidStrings("enabled", StringList::create("true,false"));
		
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    defaults_.setValidStrings("ionization_type", StringList::create("MALDI,ESI"));

    // peak and instrument parameter
    defaults_.setValue("resolution",50000,"Instrument resolution at 400Th.");

    // baseline
    defaults_.setValue("baseline:scaling",0.0,"Scale of baseline. Set to 0 to disable simulation of baseline.");
    defaults_.setMinFloat("baseline:scaling",0.0);
    defaults_.setValue("baseline:shape",0.5, "The baseline is modeled by an exponential probability density function (pdf) with f(x) = shape*e^(- shape*x)");
    defaults_.setMinFloat("baseline:shape",0.0);
    defaults_.setSectionDescription("baseline","Baseline modeling for MALDI ionization");

    // mz sampling rate
    // TODO: investigate if this can be hidden from the user by estimating it from "resolution"
    defaults_.setValue("mz:sampling_rate",0.12,"detector interval(e.g. bin size in m/z).");

    // contaminants:
    defaults_.setValue("contaminants:file","","Contaminants file with sum formula and absolute RT interval.");

    // noise params

    // VARIATION

    // m/z error
    // todo: also plan for affine trafo (as in RT shift?)
    defaults_.setValue("variation:mz:error_stddev",0.0,"Standard deviation for m/z errors. Set to 0 to disable simulation of m/z errors.");
    defaults_.setValue("variation:mz:error_mean",0.0,"Average systematic m/z error (Da)");

    defaults_.setValue("variation:intensity:scale", 1.0 , "Constant scale factor of the feature intensity. Set to 1.0 to get the real intensity values.");
    defaults_.setMinFloat("variation:intensity:scale", 0.0);
    defaults_.setValue("variation:intensity:scale_stddev", 0.0 ,"Standard deviation of peak intensity (relative to the scaled peak height). Set to 0 to get simple rescaled intensities.");
    defaults_.setMinFloat("variation:intensity:scale_stddev", 0.0);

    defaults_.setSectionDescription("variation:mz", "Shifts in mass to charge dimension of the simulated signals.");
    defaults_.setSectionDescription("variation:intensity", "Variations in intensity to model randomness in feature intensity.");
    defaults_.setSectionDescription("variation", "Random components that simulate biological and technical variations of the simulated data.");


    // shot noise
    defaults_.setValue("noise:shot:rate",0.0,"Poisson rate of shot noise per unit m/z. Set to 0 to disable simulation of shot noise.");
    defaults_.setMinFloat("noise:shot:rate",0.0);
    defaults_.setValue("noise:shot:int-mean",50.0,"Shot noise intensity mean (gaussian distributed).");

    // white noise
    defaults_.setValue("noise:white:mean", 0.0, "Mean value of the white noise that is added to each measured signal.");
    defaults_.setValue("noise:white:stddev", 0.0, "Mean value of the white noise that is added to each measured signal.");

    defaults_.setSectionDescription("noise", "Parameters modelling noise in mass spectrometry measurements.");
    defaults_.setSectionDescription("noise:shot", "Parameters regarding shot noise modelling.");
    defaults_.setSectionDescription("noise:white", "Parameters regarding white noise modelling.");

    defaultsToParam_();
  }

  void RawMSSignalSimulation::updateMembers_()
  {
    // convert from resolution @ 400th --> FWHM
    DoubleReal tmp = 400.00 / (double) param_.getValue("resolution");
    peak_std_     = (tmp / 2.355);			// Approximation for Gaussian-shaped signals
    mz_sampling_rate_ = param_.getValue("mz:sampling_rate");

    mz_error_mean_    = param_.getValue("variation:mz:error_mean");
    mz_error_stddev_  = param_.getValue("variation:mz:error_stddev");

    intensity_scale_ = param_.getValue("variation:intensity:scale");
    intensity_scale_stddev_ = param_.getValue("variation:intensity:scale_stddev");

    // contaminants:
    String file = param_.getValue("contaminants:file");
    if (file.trim().size()!=0)
    {
      if (!File::readable(file)) throw Exception::FileNotReadable(__FILE__,__LINE__,__PRETTY_FUNCTION__,file);
      // read & parse file:
      TextFile tf(file,true);
      contaminants_.clear();
      const UInt COLS_EXPECTED = 8;
      for (Size i=0;i<tf.size();++i)
      {
        if (tf[i].size()==0 || tf[i].hasPrefix("#")) continue; // skip comments
        StringList cols;
        tf[i].removeWhitespaces().split(',', cols);
        if (cols.size()!=COLS_EXPECTED) throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,tf[i],"Expected " + String(COLS_EXPECTED) + " components, got" + String(cols.size()));
        ContaminantInfo c;
        c.name = cols[0];
        try
        {
          c.sf = EmpiricalFormula(cols[1]);
        }
        catch (...)
        {
          throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,cols[1],"Could not parse line " + String(i+1) + " in " + file + ".");
        }
        try
        {
          c.rt_start = cols[2].toDouble();
          c.rt_end = cols[3].toDouble();
          c.intensity = cols[4].toDouble();
          c.q = cols[5].toInt();
        }
        catch (...)
        {
          throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,tf[i],"Could not parse line " + String(i+1) + " in " + file + ".");
        }
        if (cols[6]=="rec") c.shape = RT_RECTANGULAR;
        else if (cols[6]=="gauss") c.shape = RT_GAUSSIAN;
        else throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,tf[i], "Unknown shape type: " + cols[6] + " in line " + String(i) + " of '" + file + "'");

        if (cols[7]=="ESI") c.im = IM_ESI;
        else if (cols[7]=="MALDI") c.im = IM_MALDI;
        else if (cols[7]=="ALL") c.im = IM_ALL;
        else throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,tf[i], "Unknown ionization type: " + cols[7] + " in line " + String(i) + " of '" + file + "'");

        contaminants_.push_back(c);
      }
    }

  }

  void RawMSSignalSimulation::generateRawSignals(FeatureMapSim & features, MSSimExperiment & experiment, FeatureMapSim & c_map)
  {

    // TODO: check if signal intensities scale linear with actual abundance, e.g. DOI: 10.1021/ac0202280 for NanoFlow-ESI


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
    }
    else
    {
      for(Size idx=0; idx<features.size(); ++idx)
      {
        add2DSignal_(features[idx], experiment);
      }
      // build contaminant feature map & add raw signal
      createContaminants_(c_map, experiment);
    }

    if ((String)param_.getValue("ionization_type")=="MALDI")
    {
      addBaseLine_(experiment, minimal_mz_measurement_limit);
    }
    addShotNoise_(experiment, minimal_mz_measurement_limit, maximal_mz_measurement_limit);
    compressSignals_(experiment);

    // finally add a white noise to the simulated data
    addWhiteNoise_(experiment);
  }


  void RawMSSignalSimulation::add1DSignal_(Feature & active_feature, MSSimExperiment & experiment)
  {
    Param p1;

    SimIntensityType scale = getFeatureScaledIntensity_(active_feature.getIntensity(), 150.0);

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
    SimIntensityType scale = getFeatureScaledIntensity_(active_feature.getIntensity(), 1.0);

    Param p1;
    SimChargeType q = active_feature.getCharge();
    EmpiricalFormula ef;
    if (active_feature.metaValueExists("sum_formula"))
    {
      ef = active_feature.getMetaValue("sum_formula");
    }
    else
    {
      ef = active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula();
      std::cout << "current feature: " << active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().toString() << " with scale " << scale << std::endl;
    }
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
    // add CH and new intensity to feature
    samplePeptideModel2D_(pm, mz_start, mz_end, rt_start, rt_end, experiment, active_feature);
  }


  void RawMSSignalSimulation::samplePeptideModel1D_(const IsotopeModelGeneral & pm,
																										const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
																										MSSimExperiment & experiment, Feature & active_feature)
  {
    SimIntensityType intensity_sum = 0.0;

    LOG_DEBUG << "Sampling at [mz] " << mz_start << ":" << mz_end << std::endl;

    SimPointType point;

    for (SimCoordinateType mz = mz_start; mz < mz_end; mz += mz_sampling_rate_)
    {
      point.setMZ(mz);
      point.setIntensity( pm.getIntensity( DPosition<1>( mz ) ) );

      if ( point.getIntensity() > 10.0)
      {
        // add gaussian distributed m/z error
        double mz_err = gsl_ran_gaussian(rnd_gen_->technical_rng, mz_error_stddev_) + mz_error_mean_;
        point.setMZ( point.getMZ() + mz_err );

        intensity_sum += point.getIntensity();
        experiment[0].push_back(point);
      }
    }
    active_feature.setIntensity(intensity_sum);
  }


  void RawMSSignalSimulation::samplePeptideModel2D_(const ProductModel<2> & pm,
                                                    const SimCoordinateType mz_start,
                                                    const SimCoordinateType mz_end,
                                                    SimCoordinateType rt_start,
                                                    SimCoordinateType rt_end,
                                                    MSSimExperiment & experiment,
                                                    Feature & active_feature)
  {
    if (rt_start <=0) rt_start = 0;

    MSSimExperiment::iterator exp_iter = experiment.RTBegin(rt_start);
    if(exp_iter == experiment.end() )
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 0);
    }

    LOG_INFO << "Sampling at [RT] " << rt_start << ":" << rt_end << " [mz] " << mz_start << ":" << mz_end << std::endl;

    SimIntensityType intensity_sum = 0.0;
    vector< DPosition<2> > points;
    
    Int start_scan = exp_iter - experiment.begin();
    Int end_scan  = -5;

		SimCoordinateType rt = rt_start;

    for (; rt < rt_end && exp_iter != experiment.end(); ++exp_iter)
    {
			rt = exp_iter->getRT();
			
      for (SimCoordinateType mz = mz_start; mz < mz_end; mz += mz_sampling_rate_)
      {
        SimPointType point;
        point.setMZ(mz);
        point.setIntensity( pm.getIntensity( DPosition<2>( rt, mz) ) );

        // add gaussian distributed m/z error
        double mz_err = gsl_ran_gaussian(rnd_gen_->technical_rng, mz_error_stddev_) + mz_error_mean_;
        point.setMZ( point.getMZ() + mz_err );

        intensity_sum += point.getIntensity();
        points.push_back( DPosition<2>( rt, mz) );		// store position
        exp_iter->push_back(point);
      }
      //update last scan affected
      end_scan = exp_iter - experiment.begin();
    }

    OPENMS_POSTCONDITION(end_scan  != -5, "RawMSSignalSimulation::samplePeptideModel2D_(): setting RT bounds failed!");

    active_feature.setQuality(0,start_scan);
    active_feature.setQuality(1,end_scan);

    active_feature.setIntensity(intensity_sum);
    // store convex hull
    active_feature.getConvexHulls().clear();
    // adding ALL points of the feature
		ConvexHull2D hull;
		hull.addPoints(points);
    active_feature.getConvexHulls().push_back(hull);

  }

  void RawMSSignalSimulation::chooseElutionProfile_(EGHModel*& elutionmodel, Feature& feature, const double scale, const DoubleReal rt_sampling_rate, const MSSimExperiment & experiment)
  {
      SimCoordinateType f_rt = feature.getRT();

      Param p;
      // WARNING: step used to be 'rt_sampling_rate / 3.0', but distortion is not part of RT sim, and thus only
      //          modeled 1:1
      p.setValue("interpolation_step", rt_sampling_rate );
      p.setValue("statistics:variance", 1.0);
      p.setValue("statistics:mean", f_rt);

      p.setValue("egh:height", scale);
      p.setValue("egh:retention", f_rt);

      if (feature.metaValueExists("RT_width_gaussian"))
      { // this is for contaminants only (we want the gaussian distribution width A+B at 5% of maximal height)
        p.setValue("egh:alpha", 0.05);
        p.setValue("egh:A", double(feature.getMetaValue("RT_width_gaussian"))/2.0 * 0.9); // make width a little smaller as this is only the 5% height cutoff
        p.setValue("egh:B", double(feature.getMetaValue("RT_width_gaussian"))/2.0 * 0.9);
      }
      else
      {      // TODO remove fixed values .. come up with a meaningfull model for elutionprofile shapes
        p.setValue("egh:A", 50.0);
        p.setValue("egh:B", 60.0);
      }

      elutionmodel->setParameters(p); // does the calculation

      //----------------------------------------------------------------------

      // Hack away constness :-P   We know what we want.
      EmgModel::ContainerType &data = const_cast<EGHModel::ContainerType &>(elutionmodel->getInterpolation().getData());

      SimCoordinateType rt_em_start = elutionmodel->getInterpolation().supportMin();

      // find scan in experiment at which our elution starts
      MSSimExperiment::ConstIterator exp_it = experiment.RTBegin(rt_em_start);
      DoubleList elution_intensities;
      DoubleList elution_bounds;
      elution_bounds.resize(4); // store min and max RT (in seconds and index terms)
      elution_bounds[0] = std::distance(experiment.begin(), exp_it);
      elution_bounds[1] = exp_it->getRT();
      elution_bounds[2] = elution_bounds[0];
      elution_bounds[3] = elution_bounds[1];

      for ( Size i = 0; (i < data.size()) && exp_it!=experiment.end(); ++i, ++exp_it )
      { // .. and disturb values by (an already smoothed) distortion diced in RTSimulation
        data[i] *= (DoubleReal) exp_it->getMetaValue("distortion");
        // store elution profile in feature MetaValue
        elution_intensities.push_back(data[i]);
        elution_bounds[2] = std::distance(experiment.begin(), exp_it);
        elution_bounds[3] = exp_it->getRT();
      }
      // set elution profile details in feature -> used for MS^E precursor selection in tandemMS later
      feature.setMetaValue("elution_profile_intensities", elution_intensities);
      feature.setMetaValue("elution_profile_bounds", elution_bounds);
  }

  void RawMSSignalSimulation::createContaminants_(FeatureMapSim & c_map, MSSimExperiment & exp)
  {
    IONIZATIONMETHOD this_im = (String)param_.getValue("ionization_type")=="ESI" ? IM_ESI : IM_MALDI;
    c_map.clear(true);

    std::cerr << "\n\n CREATING " << contaminants_.size() << " CONTAMINANTS\n\n";

    for (Size i=0;i<contaminants_.size();++i)
    {
      std::cerr << "\n\n CREATING " << i << " th CONTAMINANTS with IM: " << Int(contaminants_[i].im) << " " << Int(this_im) << "\n\n";

      if (contaminants_[i].im!=IM_ALL && contaminants_[i].im!=this_im) continue;

      
      // ... create contaminants...
      FeatureMapSim::FeatureType feature;
      feature.setRT( (contaminants_[i].rt_end+contaminants_[i].rt_start)/2 );
      feature.setMZ( contaminants_[i].sf.getMonoWeight() );
      feature.setIntensity(contaminants_[i].intensity);
      if (contaminants_[i].shape == RT_RECTANGULAR)
      {
        feature.setMetaValue("RT_width_gaussian", 1e6);
      }
      else
      {
        feature.setMetaValue("RT_width_gaussian", contaminants_[i].rt_end-contaminants_[i].rt_start);
      }
      feature.setMetaValue("sum_formula", contaminants_[i].sf.getString());
      feature.setCharge(contaminants_[i].q);
      feature.setMetaValue("charge_adducts","H"+String(contaminants_[i].q));
      add2DSignal_(feature, exp);
      c_map.push_back(feature);
    }
  }


  void RawMSSignalSimulation::addShotNoise_(MSSimExperiment & experiment, SimCoordinateType minimal_mz_measurement_limit, SimCoordinateType maximal_mz_measurement_limit)
  {
    // we model the amount of (background) noise as Poisson process
    // i.e. the number of noise data points per unit m/z interval follows a Poisson
    // distribution. Noise intensity is assumed to be Gaussian-distributed.

    DoubleReal rate    = param_.getValue("noise:shot:rate");
    if (rate == 0.0) return;
    DoubleReal it_mean = param_.getValue("noise:shot:int-mean");

    const UInt num_intervals = 100;
    SimCoordinateType interval_size = ( maximal_mz_measurement_limit -  minimal_mz_measurement_limit) / num_intervals;
    SimPointType point;

    LOG_INFO << "Adding shot noise to spectra ..." << std::endl;
    LOG_INFO << "Interval size: "  << interval_size << " poisson rate: " << rate << std::endl;

    // TODO: switch to iterator ??
    for (Size i=0 ;i < experiment.size() ; ++i)
    {

      for (Size j=0;j<num_intervals;++j)
      {
        UInt counts = gsl_ran_poisson ( rnd_gen_->technical_rng, rate);
        SimCoordinateType mz_lw = j * interval_size + minimal_mz_measurement_limit;
        SimCoordinateType mz_up = (j+1) * interval_size + minimal_mz_measurement_limit;

        for (UInt c=0; c<counts;++c)
        {
          SimCoordinateType mz  = gsl_ran_flat(rnd_gen_->technical_rng, mz_lw, mz_up );
          SimCoordinateType it = gsl_ran_exponential(rnd_gen_->technical_rng,it_mean);
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

    // TODO: switch to iterator
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

  void RawMSSignalSimulation::addWhiteNoise_(MSSimExperiment &experiment)
  {
    // get white noise parameters
    DoubleReal white_noise_mean = param_.getValue("noise:white:mean");
    DoubleReal white_noise_stddev = param_.getValue("noise:white:stddev");

    for(MSSimExperiment::iterator spectrum_it = experiment.begin() ; spectrum_it != experiment.end() ; ++spectrum_it)
    {
      for(MSSimExperiment::SpectrumType::iterator peak_it = (*spectrum_it).begin() ; peak_it != (*spectrum_it).end() ; ++peak_it)
      {
        SimIntensityType intensity = peak_it->getIntensity() + gsl_ran_gaussian(rnd_gen_->technical_rng, white_noise_stddev * peak_it->getIntensity()) + white_noise_mean;
        peak_it->setIntensity( (intensity > 0.0 ? intensity : 0.0) );
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

    // this is necessary to avoid
    // 0.1 < 0.1 = true
    // due to numerical instability
    const SimCoordinateType numerical_correction = 0.0001 * mz_sampling_rate_;

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

          if ((diff_mz +  numerical_correction) < mz_sampling_rate_)
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
            //SimCoordinateType mz =  (mz1 + mz2) / 2.0;
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


  SimIntensityType RawMSSignalSimulation::getFeatureScaledIntensity_(const SimIntensityType feature_intensity, const SimIntensityType natural_scaling_factor)
  {
    SimIntensityType intensity = feature_intensity * natural_scaling_factor * intensity_scale_;
    
    // add some noise
    intensity += gsl_ran_gaussian(rnd_gen_->technical_rng, intensity_scale_stddev_ * intensity);

    return intensity;
  }
}
