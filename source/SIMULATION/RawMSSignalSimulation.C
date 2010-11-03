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
#include <OpenMS/CONCEPT/Constants.h>

#include <vector>
using std::vector;

#ifdef _OPENMP
#include <omp.h>
#endif

#include <OpenMS/FORMAT/MzMLFile.h>

namespace OpenMS {

  /**
   * TODO: review baseline and noise code
   */
  RawMSSignalSimulation::RawMSSignalSimulation(const SimRandomNumberGenerator& rng)
  : DefaultParamHandler("RawSignalSimulation"), 
    mz_sampling_rate_(),
    mz_error_mean_(),
    mz_error_stddev_(),
    intensity_scale_(),
    intensity_scale_stddev_(),
    res_model_(RES_CONSTANT),
    res_base_(0),
    rnd_gen_(&rng),
    contaminants_(),
    contaminants_loaded_(false)
  {
    setDefaultParams_();
    updateMembers_();
  }


  RawMSSignalSimulation::RawMSSignalSimulation()
    : DefaultParamHandler("RawSignalSimulation"),
      mz_sampling_rate_(),
      mz_error_mean_(),
      mz_error_stddev_(),
      intensity_scale_(),
      intensity_scale_stddev_(),
      res_model_(RES_CONSTANT),
      res_base_(0),
      contaminants_(),
      contaminants_loaded_(false)
  {
    setDefaultParams_();
    updateMembers_();
  }

  RawMSSignalSimulation::RawMSSignalSimulation(const RawMSSignalSimulation& source)
    : DefaultParamHandler(source),
    mz_sampling_rate_(source.mz_sampling_rate_),
    mz_error_mean_(source.mz_error_mean_),
    mz_error_stddev_(source.mz_error_stddev_),
    intensity_scale_(source.intensity_scale_),
    intensity_scale_stddev_(source.intensity_scale_stddev_),
    res_model_(source.res_model_),
    res_base_(source.res_base_),
    contaminants_(),
    contaminants_loaded_(false)
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

    res_model_ = source.res_model_;
    res_base_ = source.res_base_;

    contaminants_ = source.contaminants_;
    contaminants_loaded_ = source.contaminants_loaded_;

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
    defaults_.setValue("resolution:value",50000,"Instrument resolution at 400Th.");
    defaults_.setValue("resolution:type","linear","How does resolution change with increasing m/z?! QTOFs usually show 'linear' behaviour, FTs have linear degradation, and on Orbitraps the resolution decreases with square root of mass.");
    defaults_.setValidStrings("resolution:type", StringList::create("constant,linear,sqrt"));
    

    // baseline
    defaults_.setValue("baseline:scaling",0.0,"Scale of baseline. Set to 0 to disable simulation of baseline.");
    defaults_.setMinFloat("baseline:scaling",0.0);
    defaults_.setValue("baseline:shape",0.5, "The baseline is modeled by an exponential probability density function (pdf) with f(x) = shape*e^(- shape*x)");
    defaults_.setMinFloat("baseline:shape",0.0);
    defaults_.setSectionDescription("baseline","Baseline modeling for MALDI ionization");

    // mz sampling rate
    // TODO: investigate if this can be hidden from the user by estimating it from "resolution"
    //       e.g. http://www.adronsystems.com/faqs.htm#rate states 8 points per peak on low-res instruments --> ~4 points at FWHM
    defaults_.setValue("mz:sampling_rate",0.01,"detector interval(e.g. bin size in m/z).");

    // contaminants:
    defaults_.setValue("contaminants:file","examples/simulation/contaminants.csv","Contaminants file with sum formula and absolute RT interval.");

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

  DoubleReal RawMSSignalSimulation::getResolution_(const DoubleReal query_mz, const DoubleReal resolution, const RESOLUTIONMODEL model) const
  {
    switch ( model )
    {
      case RES_CONSTANT:
        return resolution;
      case RES_LINEAR:
        return resolution * (400/query_mz);
      case RES_SQRT:
        return resolution * (std::sqrt(400.0)/sqrt(query_mz));
      default:
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unknown RESOLUTIONMODEL encountered!");
    }
  }

  void RawMSSignalSimulation::updateMembers_()
  {
    res_base_ = (double) param_.getValue("resolution:value");
    String model = param_.getValue("resolution:type");
    if (model=="constant") res_model_ = RES_CONSTANT;
    else if (model=="linear") res_model_ = RES_LINEAR;
    else if (model=="sqrt") res_model_ = RES_SQRT;
    else throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Resolution:type given in parameters is unknown");

    mz_sampling_rate_ = param_.getValue("mz:sampling_rate");

    mz_error_mean_    = param_.getValue("variation:mz:error_mean");
    mz_error_stddev_  = param_.getValue("variation:mz:error_stddev");

    intensity_scale_ = param_.getValue("variation:intensity:scale");
    intensity_scale_stddev_ = param_.getValue("variation:intensity:scale_stddev");
    
    contaminants_loaded_=false;
  }

  void RawMSSignalSimulation::loadContaminants()
  {
    // contaminants:
    String contaminants_file = param_.getValue("contaminants:file");

    if (contaminants_file.trim().size()!=0)
    {
	    if (! File::readable( contaminants_file ) )
      { // look in OPENMS_DATA_PATH
        contaminants_file = File::find( contaminants_file );
      }
      if (!File::readable(contaminants_file)) throw Exception::FileNotReadable(__FILE__,__LINE__,__PRETTY_FUNCTION__,contaminants_file);
      // read & parse file:
      TextFile tf(contaminants_file,true);
      contaminants_.clear();
      const UInt COLS_EXPECTED = 8;
      for (Size i=0;i<tf.size();++i)
      {
        if (tf[i].size()==0 || tf[i].hasPrefix("#")) continue; // skip comments
        StringList cols;
        tf[i].removeWhitespaces().split(',', cols, true);
        if (cols.size()!=COLS_EXPECTED) throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,tf[i],"Expected " + String(COLS_EXPECTED) + " components, got " + String(cols.size()));
        ContaminantInfo c;
        c.name = cols[0];
        try
        {
          c.sf = EmpiricalFormula(cols[1]);
        }
        catch (...)
        {
          throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,cols[1],"Could not parse line " + String(i+1) + " in " + contaminants_file + ".");
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
          throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,tf[i],"Could not parse line " + String(i+1) + " in " + contaminants_file + ".");
        }
        if (cols[6].toUpper()=="REC") c.shape = RT_RECTANGULAR;
        else if (cols[6].toUpper()=="GAUSS") c.shape = RT_GAUSSIAN;
        else throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,tf[i], "Unknown shape type: " + cols[6] + " in line " + String(i) + " of '" + contaminants_file + "'");

        if (cols[7].toUpper()=="ESI") c.im = IM_ESI;
        else if (cols[7].toUpper()=="MALDI") c.im = IM_MALDI;
        else if (cols[7].toUpper()=="ALL") c.im = IM_ALL;
        else throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,tf[i], "Unknown ionization type: " + cols[7] + " in line " + String(i) + " of '" + contaminants_file + "'");

        contaminants_.push_back(c);
      }
    }
    contaminants_loaded_=true;
  }

  void RawMSSignalSimulation::generateRawSignals(FeatureMapSim & features, MSSimExperiment & experiment, FeatureMapSim & c_map)
  {
		LOG_INFO << "Raw MS1 Simulation ... ";
    // TODO: check if signal intensities scale linear with actual abundance, e.g. DOI: 10.1021/ac0202280 for NanoFlow-ESI


    if (param_.getValue("enabled") == "false")
    {
      LOG_INFO << "disabled" << std::endl;
			return;
		}
    else
    {
      LOG_INFO << "started" << std::endl;
    }

		// retrieve mz boundary parameters from experiment:
		SimCoordinateType minimal_mz_measurement_limit = experiment[0].getInstrumentSettings().getScanWindows()[0].begin;
		SimCoordinateType maximal_mz_measurement_limit = experiment[0].getInstrumentSettings().getScanWindows()[0].end;
		
    LOG_INFO << "  Simulating signal for " << features.size() << " features ..." << std::endl;

    this->startProgress(0,features.size(),"RawMSSignal");

    Size progress(0);
    // we have a bit of code duplication here but this eases the parallelization
    // step
    if(experiment.size() == 1) // MS only
    {
      for(FeatureMap< >::iterator feature_it = features.begin();
          feature_it != features.end();
          ++feature_it,++progress)
      {
        add1DSignal_(*feature_it,experiment);
        this->setProgress(progress);
      }
    }
    else // LC/MS
    {
      std::vector < MSSimExperiment*> experiments; // pointer to experiment(s)
      experiments.push_back(&experiment); // the master thread gets the original (just a reference, no copying here)

#ifdef _OPENMP
      // prepare random numbers for the different threads
      // each possible thread gets his own set of random
      // numbers
      Size thread_count = omp_get_max_threads();

      threaded_random_numbers_.resize(thread_count);
      threaded_random_numbers_index_.resize(thread_count);
      experiments.reserve(thread_count); // !reserve! 
      std::vector < MSSimExperiment> experiments_tmp(thread_count-1); // holds MSExperiments for slave threads

      for(Size i=0 ; i<thread_count; ++i)
      {
        threaded_random_numbers_[i].resize(THREADED_RANDOM_NUMBER_POOL_SIZE);
        threaded_random_numbers_index_[i] = THREADED_RANDOM_NUMBER_POOL_SIZE;
      }

      if (thread_count>1)
      {
        // prepare a temporary experiment to store the results
        MSSimExperiment e_tmp = experiment;
        // remove actual data
        for(MSSimExperiment::Iterator temp_it = e_tmp.begin(); temp_it != e_tmp.end(); ++temp_it)
        {
          temp_it->clear(false);
        }
        // each slave thread gets a copy
        for (Size i=1; i<thread_count; ++i)
        {
          experiments_tmp[i-1] = e_tmp;
          // assign it to the list of experiments (this is no real copy, but a reference!)
          experiments.push_back(&(experiments_tmp[i-1]));
        }
      }
#else
      Size thread_count=1;
#endif

      Size compress_size_intermediate = 20000 / thread_count; // compress map every X features, (10.000 feature are ~ 2 GB at 0.002 sampling rate)
      Size compress_count = 0; // feature count (for each thread)
      int current_thread(0);

#pragma omp parallel for private(current_thread) firstprivate(compress_count)
      for (SignedSize f = 0 ; f < (SignedSize)features.size() ; ++f)
      {
        #ifdef _OPENMP // update experiment index if necessary
          current_thread = omp_get_thread_num();
        #endif
        add2DSignal_(features[f], *(experiments[current_thread]));

        // progresslogger, only master thread sets progress (no barrier here)
        #pragma omp atomic
        ++progress;
        if (current_thread == 0) this->setProgress(progress);

        // intermediate compress to avoid memory problems
        ++compress_count;
        if (compress_count > compress_size_intermediate)
        {
          compress_count = 0;
          compressSignals_(*(experiments[current_thread]));
        }
      } // ! raw signal sim

#ifdef _OPENMP // merge back other experiments
      for (Size i=1;i<experiments.size();++i)
      {
        MSSimExperiment::Iterator org_it = experiment.begin();
        MSSimExperiment::Iterator temp_it = experiments[i]->begin();

        // copy peak data from temporal experiment
        for(; org_it != experiment.end() ; ++org_it, ++temp_it)
        {
          if(temp_it->empty()) continue; // we do not care if the spectrum wasn't touched at all
          // append all points from temp to org
          org_it->insert(org_it->end(), temp_it->begin(), temp_it->end());
          // delete from child experiment to save memory (otherwise the merge would double it!)
          temp_it->clear(false);
        }
      }
#endif
    
    } // ! 1D or 2D

    this->endProgress();

    // finally sort generated data
    experiment.sortSpectra(true);
    experiment.updateRanges();

    // build contaminant feature map & add raw signal
    if (experiment.size() > 1) // LC/MS only currently
    {
      createContaminants_(c_map, experiment);
    }

    if ((String)param_.getValue("ionization_type")=="MALDI")
    {
      addBaseLine_(experiment, minimal_mz_measurement_limit);
    }
    addShotNoise_(experiment, minimal_mz_measurement_limit, maximal_mz_measurement_limit);
    compressSignals_(experiment);

    // finally add white noise to the simulated data
    addWhiteNoise_(experiment);
  }

  DoubleReal RawMSSignalSimulation::getPeakSD_(const DoubleReal mz) const
  {
    // convert from resolution @ current m/z --> FWHM
    // Approximation for Gaussian-shaped signals,
    // i.e. sqrt(2*ln(2))*2 = 2.35482
    // , relating FWHM to gaussian width
    DoubleReal tmp = mz / getResolution_(mz, res_base_, res_model_);
    DoubleReal peak_std = (tmp / 2.35482);
    return peak_std;
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
    p1.setValue("isotope:stdev", getPeakSD_(active_feature.getMZ()));
    p1.setValue("intensity_scaling", scale);
    p1.setValue("charge", q);

    IsotopeModel isomodel;
    isomodel.setParameters(p1);
    isomodel.setSamples(ef);

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
      //std::cout << "current feature: " << active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().toString() << " with scale " << scale << std::endl;
    }
    ef += active_feature.getMetaValue("charge_adducts"); // adducts
    ef -= String("H")+String(q);ef.setCharge(q);				 // effectively substract q electrons
    p1.setValue("statistics:mean", ef.getAverageWeight() / q);
    p1.setValue("interpolation_step", 0.001);
    p1.setValue("isotope:stdev", getPeakSD_(active_feature.getMZ()));
    p1.setValue("charge", q);

    IsotopeModel* isomodel = new IsotopeModel();
    isomodel->setParameters(p1); // this needs to come BEFORE setSamples() - otherwise the default setSamples() is called here!
    isomodel->setSamples(ef); // this already includes adducts

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

    // start and end points of the sampling
    SimCoordinateType rt_start ( elutionmodel->getInterpolation().supportMin() );
    SimCoordinateType rt_end ( elutionmodel->getInterpolation().supportMax() );
    if (active_feature.metaValueExists("RT_width_start") && active_feature.metaValueExists("RT_width_end"))
    { // this is a contaminant with sampling restrictions
      rt_start=active_feature.getMetaValue("RT_width_start");
      rt_end=active_feature.getMetaValue("RT_width_end");
    }
    SimCoordinateType mz_start ( isomodel->getInterpolation().supportMin() );
    SimCoordinateType mz_end ( isomodel->getInterpolation().supportMax() );

    // add peptide to GLOBAL MS map
    // add CH and new intensity to feature
    samplePeptideModel2D_(pm, mz_start, mz_end, rt_start, rt_end, experiment, active_feature);
  }


  void RawMSSignalSimulation::samplePeptideModel1D_(const IsotopeModel & pm,
																										const SimCoordinateType mz_start,  const SimCoordinateType mz_end,
																										MSSimExperiment & experiment, Feature & active_feature)
  {
    SimIntensityType intensity_sum = 0.0;

    //LOG_DEBUG << "Sampling at [mz] " << mz_start << ":" << mz_end << std::endl;

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

    SimIntensityType intensity_sum = 0.0;
    vector< DPosition<2> > points;
    
    Int start_scan = exp_iter - experiment.begin();
    Int end_scan  = -5;

		SimCoordinateType rt = rt_start;
    // Sample the model ...
    for (; rt < rt_end && exp_iter != experiment.end(); ++exp_iter)
    {
			rt = exp_iter->getRT();
			
      for (SimCoordinateType mz = mz_start; mz < mz_end; mz += mz_sampling_rate_)
      {
        ProductModel<2>::IntensityType intensity = pm.getIntensity( DPosition<2>( rt, mz) );

        if(intensity < 1) // intensity cutoff (below that we don't want to see a signal)
        {
          continue;
        }

        SimPointType point;
        point.setMZ(mz);
        point.setIntensity( intensity );

        //LOG_ERROR << "Sampling " << rt << " , " << mz << " -> " << point.getIntensity() << std::endl;

        // add gaussian distributed m/z error
        double mz_err = 0.0;
#ifdef _OPENMP
        int CURRENT_THREAD = omp_get_thread_num();
        // check if we need to refill the random number pool for this thread
        if(threaded_random_numbers_index_[ CURRENT_THREAD ] == THREADED_RANDOM_NUMBER_POOL_SIZE)
        {
          if(mz_error_stddev_ != 0.0)
          {
            #pragma omp critical(generate_random_number_for_thread)
            {
              for(Size i = 0 ; i < THREADED_RANDOM_NUMBER_POOL_SIZE ; ++i)
              {
                threaded_random_numbers_[CURRENT_THREAD][i] = gsl_ran_gaussian(rnd_gen_->technical_rng, mz_error_stddev_) + mz_error_mean_;
              }
            }
          }
          else
          {
            // we do not need to care about concurrency here
            fill(threaded_random_numbers_[CURRENT_THREAD].begin(), threaded_random_numbers_[CURRENT_THREAD].end() , mz_error_mean_);
          }
          // reset index for this thread to first position
          threaded_random_numbers_index_[CURRENT_THREAD] = 0;
        }

        mz_err = threaded_random_numbers_[CURRENT_THREAD][threaded_random_numbers_index_[CURRENT_THREAD]++];
#else
        // we can use the normal gaussian ran-gen if we do not use OPENMP
        mz_err = gsl_ran_gaussian(rnd_gen_->technical_rng, mz_error_stddev_) + mz_error_mean_;
#endif
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

  void RawMSSignalSimulation::chooseElutionProfile_(EGHModel* const elutionmodel, Feature& feature, const double scale, const DoubleReal rt_sampling_rate, const MSSimExperiment & experiment)
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
        p.setValue("egh:A", 10.0);
        p.setValue("egh:B", 10.0);
      }

      elutionmodel->setParameters(p); // does the calculation
      //----------------------------------------------------------------------

      // Hack away constness :-P   We know what we want.
      EGHModel::ContainerType &data = const_cast<EGHModel::ContainerType &>(elutionmodel->getInterpolation().getData());

      SimCoordinateType rt_em_start = elutionmodel->getInterpolation().supportMin();

      // find scan in experiment at which our elution starts
      MSSimExperiment::ConstIterator exp_it = experiment.RTBegin(rt_em_start);
      if (exp_it==experiment.end()) --exp_it; // we need the last valid RT below, so .end() is not useful

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
    if(exp.size() == 1)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__); // not implemented for 1D yet
    }

    if (!contaminants_loaded_) loadContaminants();

    IONIZATIONMETHOD this_im = (String)param_.getValue("ionization_type")=="ESI" ? IM_ESI : IM_MALDI;
    c_map.clear(true);

    Size out_of_range_RT(0),out_of_range_MZ(0);
		SimCoordinateType minimal_mz_measurement_limit = exp[0].getInstrumentSettings().getScanWindows()[0].begin;
		SimCoordinateType maximal_mz_measurement_limit = exp[0].getInstrumentSettings().getScanWindows()[0].end;

    for (Size i=0;i<contaminants_.size();++i)
    {
      if (contaminants_[i].im!=IM_ALL && contaminants_[i].im!=this_im) continue;

      if (exp.getMinRT() > contaminants_[i].rt_end || contaminants_[i].rt_start > exp.getMaxRT())
      {
        ++out_of_range_RT;
        continue;
      }
      // ... create contaminants...
      FeatureMapSim::FeatureType feature;
      feature.setRT( (contaminants_[i].rt_end+contaminants_[i].rt_start) / 2 );
      feature.setMZ( (contaminants_[i].sf.getMonoWeight() / contaminants_[i].q) + Constants::PROTON_MASS_U ); // m/z (incl. protons)
      if (!(minimal_mz_measurement_limit < feature.getMZ() && feature.getMZ() < maximal_mz_measurement_limit))
      {
        ++out_of_range_MZ;
        continue;
      }
      feature.setIntensity(contaminants_[i].intensity);
      if (contaminants_[i].shape == RT_RECTANGULAR)
      {
        feature.setMetaValue("RT_width_gaussian", 1e6);
        feature.setMetaValue("RT_width_start", contaminants_[i].rt_start);
        feature.setMetaValue("RT_width_end", contaminants_[i].rt_end);
      }
      else
      {
        feature.setMetaValue("RT_width_gaussian", contaminants_[i].rt_end-contaminants_[i].rt_start);
      }
      feature.setMetaValue("sum_formula", contaminants_[i].sf.getString()); // formula without adducts
      feature.setCharge(contaminants_[i].q);
      feature.setMetaValue("charge_adducts","H"+String(contaminants_[i].q));  // adducts separately
      add2DSignal_(feature, exp);
      c_map.push_back(feature);
    }

    LOG_INFO << "Contaminants out-of-RT-range: " << out_of_range_RT << " / " << contaminants_.size() << "\n";
    LOG_INFO << "Contaminants out-of-MZ-range: " << out_of_range_MZ << " / " << contaminants_.size() << "\n";

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
    LOG_INFO << "Interval size: "  << interval_size << ", poisson rate: " << rate << std::endl;

    for (MSSimExperiment::Iterator it_exp=experiment.begin(); it_exp != experiment.end() ; ++it_exp)
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
          if(it > 0.0)
          {
            point.setIntensity(it);
            point.setMZ(mz);
            it_exp->push_back(point);
          }
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
        SimIntensityType intensity = peak_it->getIntensity() + gsl_ran_gaussian(rnd_gen_->technical_rng, white_noise_stddev) + white_noise_mean;
        peak_it->setIntensity( (intensity > 0.0 ? intensity : 1.0) );
      }
    }
  }

  // TODO: add instrument specific sampling technique
  void RawMSSignalSimulation::compressSignals_(MSSimExperiment & experiment)
  {
		// assume we changed every scan a priori

    Size count = 0;
    do
    { // this loop is not required any longer (one run will always suffice) -- but for now we leave it in
      count = compressSignalsRun_(experiment);
      std::cerr << "compress called with " << count << "\n";
    } while (count != 0);
  }

  Size RawMSSignalSimulation::compressSignalsRun_(MSSimExperiment & experiment)
  {
    //SimCoordinateType diff_mz = 0.0;
    SimPointType p;

    Size count = 0;
    bool change = false;

    // this is necessary to avoid
    // 0.1 < 0.1 = true
    // due to numerical instability
    const SimCoordinateType mz_sampling_rate_num = 0.99 * mz_sampling_rate_; // be generous here (99%) - everything closer to 1 might be unstable

    Size point_count_before = 0, point_count_after = 0;

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
          Size point_count = 1;
          DoubleReal mz_start = experiment[i][j].getMZ();
          DoubleReal mz_sum = mz_start;
          DoubleReal int_sum = experiment[i][j].getIntensity();
          change = false;
          while ((j < (experiment[i].size()-1)) // join points within sampling rate
              && (experiment[i][j+1].getMZ() < (mz_start + mz_sampling_rate_num)))
          {
            ++j;
            ++point_count;
            change = true;
            mz_sum += experiment[i][j].getMZ();
            int_sum += experiment[i][j].getIntensity();
          }
          if (change)
          {
            p.setIntensity(int_sum);
            p.setMZ( mz_sum/point_count ); // use average m/z
						cont.push_back(p);
          }
          else cont.push_back( experiment[i][j] );
/*
// OLD Version
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
					}*/
				}
				// don't forget the last one
				if (!change) cont.push_back( experiment[i][ (experiment[i].size() - 1) ] );

        point_count_before += experiment[i].size();
	      experiment[i] = cont;
        point_count_after += experiment[i].size();

			}
    }

    std::cerr << " point count: " <<  point_count_before << " --> " << point_count_after << " (" << (point_count_after*100/point_count_before) << "%)\n";

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
    // TODO: variables model für den intensitäts-einfluss
    // e.g. sqrt(intensity) || ln(intensity)
    intensity += gsl_ran_gaussian(rnd_gen_->technical_rng, intensity_scale_stddev_ * intensity);

    return intensity;
  }
}
