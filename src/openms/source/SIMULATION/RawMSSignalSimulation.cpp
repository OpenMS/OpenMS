// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------


#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/random/uniform_real.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions.hpp>

#ifdef _OPENMP
#include <omp.h>

#endif




namespace OpenMS
{

  /**
   * TODO: review baseline and noise code
   */
  RawMSSignalSimulation::RawMSSignalSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr rng) :
    DefaultParamHandler("RawSignalSimulation"),
    ProgressLogger(),
    mz_error_mean_(),
    mz_error_stddev_(),
    intensity_scale_(),
    intensity_scale_stddev_(),
    res_model_(RES_CONSTANT),
    res_base_(0),
    rnd_gen_(rng),
    contaminants_(),
    contaminants_loaded_(false)
  {
    setDefaultParams_();
    updateMembers_();
  }

  RawMSSignalSimulation::RawMSSignalSimulation() :
    DefaultParamHandler("RawSignalSimulation"),
    mz_error_mean_(),
    mz_error_stddev_(),
    intensity_scale_(),
    intensity_scale_stddev_(),
    res_model_(RES_CONSTANT),
    res_base_(0),
    rnd_gen_(),
    contaminants_(),
    contaminants_loaded_(false)
  {
    setDefaultParams_();
    updateMembers_();
  }

  RawMSSignalSimulation::RawMSSignalSimulation(const RawMSSignalSimulation& source) :
    DefaultParamHandler(source),
    ProgressLogger(source),
    mz_error_mean_(source.mz_error_mean_),
    mz_error_stddev_(source.mz_error_stddev_),
    intensity_scale_(source.intensity_scale_),
    intensity_scale_stddev_(source.intensity_scale_stddev_),
    res_model_(source.res_model_),
    res_base_(source.res_base_),
    contaminants_(),
    contaminants_loaded_(false)
  {
    setParameters(source.getParameters());
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
  }

  RawMSSignalSimulation& RawMSSignalSimulation::operator=(const RawMSSignalSimulation& source)
  {
    setParameters(source.getParameters());
    rnd_gen_ = source.rnd_gen_;

    mz_error_mean_ = source.mz_error_mean_;
    mz_error_stddev_ = source.mz_error_stddev_;

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
    defaults_.setValue("enabled", "true", "Enable RAW signal simulation? (select 'false' if you only need feature-maps)");
    defaults_.setValidStrings("enabled", ListUtils::create<String>("true,false"));

    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    defaults_.setValidStrings("ionization_type", ListUtils::create<String>("MALDI,ESI"));

    // peak and instrument parameter
    defaults_.setValue("resolution:value", 50000, "Instrument resolution at 400 Th");
    defaults_.setValue("resolution:type", "linear", "How does resolution change with increasing m/z?! QTOFs usually show 'constant' behavior, FTs have linear degradation, and on Orbitraps the resolution decreases with square root of mass");
    defaults_.setValidStrings("resolution:type", ListUtils::create<String>("constant,linear,sqrt"));

    defaults_.setValue("peak_shape", "Gaussian", "Peak Shape used around each isotope peak (be aware that the area under the curve is constant for both types, but the maximal height will differ (~ 2:3 = Lorentz:Gaussian) due to the wider base of the Lorentzian");
    defaults_.setValidStrings("peak_shape", ListUtils::create<String>("Gaussian,Lorentzian"));


    // baseline
    defaults_.setValue("baseline:scaling", 0.0, "Scale of baseline. Set to 0 to disable simulation of baseline");
    defaults_.setMinFloat("baseline:scaling", 0.0);
    defaults_.setValue("baseline:shape", 0.5, "The baseline is modeled by an exponential probability density function (pdf) with f(x) = shape*e^(- shape*x)");
    defaults_.setMinFloat("baseline:shape", 0.0);
    defaults_.setSectionDescription("baseline", "Baseline modeling for MALDI ionization");

    // mz sampling rate
    //       e.g. http://www.adronsystems.com/faqs.htm#rate states 8 points per peak on low-res instruments --> ~4 points at FWHM
    defaults_.setValue("mz:sampling_points", 3, "Number of raw data points per FWHM of the peak");
    defaults_.setMinInt("mz:sampling_points", 2);

    // contaminants:
    defaults_.setValue("contaminants:file", "examples/simulation/contaminants.csv", "Contaminants file with sum formula and absolute RT interval. See 'OpenMS/examples/simulation/contaminants.txt' for details");

    // VARIATION

    // m/z error
    // todo: also plan for affine trafo (as in RT shift?)
    defaults_.setValue("variation:mz:error_mean", 0.0, "Average systematic m/z error (in Da)");
    defaults_.setValue("variation:mz:error_stddev", 0.0, "Standard deviation for m/z errors. Set to 0 to disable simulation of m/z errors");
    defaults_.setSectionDescription("variation:mz", "Shifts in mass to charge dimension of the simulated signals");

    defaults_.setValue("variation:intensity:scale", 100.0, "Constant scale factor of the feature intensity. Set to 1.0 to get the real intensity values provided in the FASTA file");
    defaults_.setMinFloat("variation:intensity:scale", 0.0);
    defaults_.setValue("variation:intensity:scale_stddev", 0.0, "Standard deviation of peak intensity (relative to the scaled peak height). Set to 0 to get simple rescaled intensities");
    defaults_.setMinFloat("variation:intensity:scale_stddev", 0.0);
    defaults_.setSectionDescription("variation:intensity", "Variations in intensity to model randomness in feature intensity");

    defaults_.setSectionDescription("variation", "Random components that simulate biological and technical variations of the simulated data");

    // NOISE

    // shot noise
    // we model the amount of (background) noise as Poisson process
    // i.e. the number of noise data points per unit m/z interval follows a Poisson
    // distribution. Noise intensity is assumed to be exponentially-distributed.
    defaults_.setValue("noise:shot:rate", 0.0, "Poisson rate of shot noise per unit m/z (random peaks in m/z, where the number of peaks per unit m/z follows a Poisson distribution). Set this to 0 to disable simulation of shot noise");
    defaults_.setMinFloat("noise:shot:rate", 0.0);
    defaults_.setValue("noise:shot:intensity-mean", 1.0, "Shot noise intensity mean (exponentially distributed with given mean)");
    defaults_.setSectionDescription("noise:shot", "Parameters of Poisson and Exponential for shot noise modeling (set :rate OR :mean = 0 to disable)");

    // white noise
    defaults_.setValue("noise:white:mean", 0.0, "Mean value of white noise (Gaussian) being added to each *measured* signal intensity");
    defaults_.setValue("noise:white:stddev", 0.0, "Standard deviation of white noise being added to each *measured* signal intensity");
    defaults_.setSectionDescription("noise:white", "Parameters of Gaussian distribution for white noise modeling (set :mean AND :stddev = 0 to disable). No new peaks are generated; only intensity of existing ones is changed");

    // detector noise
    defaults_.setValue("noise:detector:mean", 0.0, "Mean intensity value of the detector noise (Gaussian distribution)");
    defaults_.setValue("noise:detector:stddev", 0.0, "Standard deviation of the detector noise (Gaussian distribution)");
    defaults_.setSectionDescription("noise:detector", "Parameters of Gaussian distribution for detector noise modeling (set :mean AND :stddev = 0 to disable). If enabled, ALL possible m/z positions (up to sampling frequency of detector) will receive an intensity increase/decrease according to the specified Gaussian intensity distribution (similar to a noisy baseline)");

    defaults_.setSectionDescription("noise", "Parameters modeling noise in mass spectrometry measurements");

    defaultsToParam_();
  }

  double RawMSSignalSimulation::getResolution_(const double query_mz, const double resolution, const RESOLUTIONMODEL model) const
  {
    switch (model)
    {
    case RES_CONSTANT:
      return resolution;

    case RES_LINEAR:
      return resolution * (400 / query_mz);

    case RES_SQRT:
      return resolution * (std::sqrt(400.0) / sqrt(query_mz));

    default:
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown RESOLUTIONMODEL encountered!");
    }
  }

  void RawMSSignalSimulation::updateMembers_()
  {
    res_base_ = (double) param_.getValue("resolution:value");
    String model = param_.getValue("resolution:type");
    if (model == "constant")
      res_model_ = RES_CONSTANT;
    else if (model == "linear")
      res_model_ = RES_LINEAR;
    else if (model == "sqrt")
      res_model_ = RES_SQRT;
    else
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Resolution:type given in parameters is unknown");

    sampling_points_per_FWHM_ = (Int) param_.getValue("mz:sampling_points") - 1;

    mz_error_mean_    = param_.getValue("variation:mz:error_mean");
    mz_error_stddev_  = param_.getValue("variation:mz:error_stddev");

    intensity_scale_ = param_.getValue("variation:intensity:scale");
    intensity_scale_stddev_ = param_.getValue("variation:intensity:scale_stddev");

    contaminants_loaded_ = false;
  }

  void RawMSSignalSimulation::loadContaminants()
  {
    // contaminants:
    String contaminants_file = param_.getValue("contaminants:file");

    if (contaminants_file.trim().size() != 0)
    {
      if (!File::readable(contaminants_file)) // look in OPENMS_DATA_PATH
      {
        contaminants_file = File::find(contaminants_file);
      }
      if (!File::readable(contaminants_file))
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, contaminants_file);
      // read & parse file:
      TextFile tf(contaminants_file, true);
      contaminants_.clear();
      const UInt COLS_EXPECTED = 8;
      Size line_number = 1;
      for (TextFile::ConstIterator tf_it = tf.begin(); tf_it != tf.end(); ++tf_it, ++line_number)
      {
        if (tf_it->empty() || tf_it->hasPrefix("#")) continue; // skip comments

        String line = *tf_it;
        StringList cols;
        line.removeWhitespaces().split(',', cols, true);
        if (cols.size() != COLS_EXPECTED)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "Expected " + String(COLS_EXPECTED) + " components, got " + String(cols.size()));
        }
        ContaminantInfo c;
        c.name = cols[0];
        try
        {
          c.sf = EmpiricalFormula(cols[1]);
          if (c.sf.getCharge() != 0)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, cols[1], "Line " + String(line_number) + " in " + contaminants_file + " contains forbidden charged sum formulas. Charges must be specified in another column. Remove all '+' or '-'!");
          }
        }
        catch (...)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, cols[1], "Could not parse line " + String(line_number) + " in " + contaminants_file + ".");
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
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "Could not parse line " + String(line_number) + " in " + contaminants_file + ".");
        }
        if (cols[6].toUpper() == "REC")
        {
          c.shape = RT_RECTANGULAR;
        }
        else if (cols[6].toUpper() == "GAUSS")
        {
          c.shape = RT_GAUSSIAN;
        }
        else
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "Unknown shape type: " + cols[6] + " in line " + String(line_number) + " of '" + contaminants_file + "'");
        }

        if (cols[7].toUpper() == "ESI")
        {
          c.im = IM_ESI;
        }
        else if (cols[7].toUpper() == "MALDI")
        {
          c.im = IM_MALDI;
        }
        else if (cols[7].toUpper() == "ALL")
        {
          c.im = IM_ALL;
        }
        else
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "Unknown ionization type: " + cols[7] + " in line " + String(line_number) + " of '" + contaminants_file + "'");
        }

        contaminants_.push_back(c);
      }
    }
    contaminants_loaded_ = true;
  }

  void RawMSSignalSimulation::generateRawSignals(SimTypes::FeatureMapSim& features, SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& experiment_ct, SimTypes::FeatureMapSim& c_map)
  {
    OPENMS_LOG_INFO << "Raw MS1 Simulation ... ";
    // TODO: check if signal intensities scale linear with actual abundance, e.g. DOI: 10.1021/ac0202280 for NanoFlow-ESI

    // we rely on the same size of Raw and Peak Map
    if (experiment.size() != experiment_ct.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, experiment_ct.size());
    }

    if (param_.getValue("enabled") == "false")
    {
      OPENMS_LOG_INFO << "disabled" << std::endl;
      return;
    }
    else
    {
      OPENMS_LOG_INFO << "started" << std::endl;
    }

    // retrieve mz boundary parameters from experiment:
    SimTypes::SimCoordinateType minimal_mz_measurement_limit = experiment[0].getInstrumentSettings().getScanWindows()[0].begin;
    SimTypes::SimCoordinateType maximal_mz_measurement_limit = experiment[0].getInstrumentSettings().getScanWindows()[0].end;

    // grid is constant over scans, so we compute it only once
    getSamplingGrid_(grid_, minimal_mz_measurement_limit, maximal_mz_measurement_limit, 5); // every 5 Da we adjust the sampling width by local FWHM

    OPENMS_LOG_INFO << "  Simulating signal for " << features.size() << " features ..." << std::endl;

    this->startProgress(0, features.size(), "RawMSSignal");

    Size progress(0);
    // we have a bit of code duplication here but this eases the
    // parallelization step
    if (experiment.size() == 1) // MS only
    {
      for (FeatureMap::iterator feature_it = features.begin();
           feature_it != features.end();
           ++feature_it, ++progress)
      {
        add1DSignal_(*feature_it, experiment, experiment_ct);
        this->setProgress(progress);
      }
    }
    else // LC/MS
    {
      std::vector<SimTypes::MSSimExperiment*> experiments; // pointer to experiment(s)
      experiments.push_back(&experiment); // the master thread gets the original (just a reference, no copying here)

      std::vector<SimTypes::MSSimExperiment*> experiments_ct; // pointer to experiment(s)
      experiments_ct.push_back(&experiment_ct); // the master thread gets the original (just a reference, no copying here)


#ifdef _OPENMP
      // prepare random numbers for the different threads
      // each possible thread gets his own set of random
      // numbers
      Size thread_count = omp_get_max_threads();

      threaded_random_numbers_.resize(thread_count);
      threaded_random_numbers_index_.resize(thread_count);
      experiments.reserve(thread_count); // !reserve!
      experiments_ct.reserve(thread_count); // !reserve!
      std::vector<SimTypes::MSSimExperiment> experiments_tmp(thread_count - 1); // holds MSExperiments for slave threads
      std::vector<SimTypes::MSSimExperiment> experiments_ct_tmp(thread_count - 1); // holds MSExperiments (centroided) for slave threads

      for (Size i = 0; i < thread_count; ++i)
      {
        threaded_random_numbers_[i].resize(THREADED_RANDOM_NUMBER_POOL_SIZE_);
        threaded_random_numbers_index_[i] = THREADED_RANDOM_NUMBER_POOL_SIZE_;
      }

      if (thread_count > 1)
      {
        // prepare a temporary experiment to store the results
        SimTypes::MSSimExperiment e_tmp = experiment;
        SimTypes::MSSimExperiment e_ct_tmp = experiment_ct;
        // remove actual data
        for (Size i = 0; i < e_tmp.size(); ++i)
        {
          e_tmp[i].clear(false);
          e_ct_tmp[i].clear(false);
        }
        // each slave thread gets a copy
        for (Size i = 1; i < thread_count; ++i)
        {
          experiments_tmp[i - 1] = e_tmp;
          experiments_ct_tmp[i - 1] = e_ct_tmp;
          // assign it to the list of experiments (this is no real copy, but a reference!)
          experiments.push_back(&(experiments_tmp[i - 1]));
          experiments_ct.push_back(&(experiments_ct_tmp[i - 1]));
        }
      }
#else
      Size thread_count = 1;
#endif

      Size compress_size_intermediate = 20000 / thread_count; // compress map every X features, (10.000 feature are ~ 2 GB at 0.002 sampling rate)
      Size compress_count = 0; // feature count (for each thread)

#ifdef _OPENMP
#pragma omp parallel for firstprivate(compress_count)
#endif
      for (SignedSize f = 0; f < (SignedSize)features.size(); ++f)
      {
#ifdef _OPENMP // update experiment index if necessary
        const int current_thread = omp_get_thread_num();
#else
        const int current_thread(0);
#endif
        add2DSignal_(features[f], *(experiments[current_thread]), *(experiments_ct[current_thread]));

        // progresslogger, only master thread sets progress (no barrier here)
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++progress;
        if (current_thread == 0)
        {
          this->setProgress(progress);
        }

        // intermediate compress to avoid memory problems
        ++compress_count;
        if (compress_count > compress_size_intermediate)
        {
          compress_count = 0;
          compressSignals_(*(experiments[current_thread]));
        }
      } // ! raw signal sim

#ifdef _OPENMP // merge back other experiments
      for (Size i = 1; i < experiments.size(); ++i)
      {
        // copy peak data from temporal experiment
        for (Size scan = 0; scan < experiment.size(); ++scan)
        {
          if ((*experiments[i])[scan].empty())
            continue; // we do not care if the spectrum wasn't touched at all
          // append all points from temp to org
          experiment[scan].insert(experiment[scan].end(), (*experiments[i])[scan].begin(), (*experiments[i])[scan].end());
          // delete from child experiment to save memory (otherwise the merge would double it!)
          (*experiments[i])[scan].clear(false);

          // peak GT ( small, so no need to compress)
          experiment_ct[scan].insert(experiment_ct[scan].end(), (*experiments_ct[i])[scan].begin(), (*experiments_ct[i])[scan].end());

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
      createContaminants_(c_map, experiment, experiment_ct);
    }

    if ((String)param_.getValue("ionization_type") == "MALDI")
    {
      addBaseLine_(experiment, minimal_mz_measurement_limit);
    }
    addShotNoise_(experiment, minimal_mz_measurement_limit, maximal_mz_measurement_limit);
    compressSignals_(experiment);

    // add white noise to the simulated data
    addWhiteNoise_(experiment);

    // add detector noise the simulated data
    addDetectorNoise_(experiment);
  }

  double RawMSSignalSimulation::getPeakWidth_(const double mz, const bool is_gaussian) const
  {
    double mz_local = std::max(mz, 400.0); // at least assume m/z=400, as otherwise FWHM might get ridiculously small
    // convert from resolution @ current m/z --> FWHM
    double fwhm = mz_local / getResolution_(mz_local, res_base_, res_model_);
    // Approximation for Gaussian-shaped signals,
    // i.e. sqrt(2*ln(2))*2 = 2.35482
    // , relating FWHM to Gaussian width
    if (is_gaussian)
      fwhm /= 2.35482;
    else
    {
    } // for Lorentzian, we do nothing as the scale parameter is exactly the FWHM
    return fwhm;
  }

  void RawMSSignalSimulation::add1DSignal_(Feature& active_feature, SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& experiment_ct)
  {
    SimTypes::SimIntensityType scale = getFeatureScaledIntensity_(active_feature.getIntensity(), 100.0);

    SimTypes::SimChargeType q = active_feature.getCharge();
    EmpiricalFormula ef = active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula();
    ef += EmpiricalFormula(active_feature.getMetaValue("charge_adducts")); // adducts
    ef -= EmpiricalFormula(String("H") + String(q));
    ef.setCharge(q); // effectively subtract q electrons

    Param p1;
    p1.setValue("statistics:mean", ef.getAverageWeight() / q);
    p1.setValue("interpolation_step", 0.001);
    p1.setValue("isotope:mode:mode", param_.getValue("peak_shape"));
    p1.setValue("intensity_scaling", 0.001 * scale); // this removes the problem of to big isotope-model values
    p1.setValue("charge", q);
    double fwhm;
    if (param_.getValue("peak_shape") == "Gaussian")
    {
      fwhm = getPeakWidth_(active_feature.getMZ(), true);
      p1.setValue("isotope:mode:GaussianSD", fwhm);
    }
    else
    {
      fwhm = getPeakWidth_(active_feature.getMZ(), false);
      p1.setValue("isotope:mode:LorentzFWHM", fwhm);
    }

    IsotopeModel isomodel;
    isomodel.setParameters(p1);
    isomodel.setSamples(ef);

    SimTypes::SimCoordinateType mz_start = isomodel.getInterpolation().supportMin();
    SimTypes::SimCoordinateType mz_end = isomodel.getInterpolation().supportMax();

    samplePeptideModel1D_(isomodel, mz_start, mz_end, experiment, experiment_ct, active_feature);
  }

  void RawMSSignalSimulation::add2DSignal_(Feature& active_feature, SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& experiment_ct)
  {
    SimTypes::SimIntensityType scale = getFeatureScaledIntensity_(active_feature.getIntensity(), 1.0);

    SimTypes::SimChargeType q = active_feature.getCharge();
    EmpiricalFormula ef;
    if (active_feature.metaValueExists("sum_formula"))
    {
      ef = EmpiricalFormula(active_feature.getMetaValue("sum_formula"));
    }
    else
    {
      ef = EmpiricalFormula(active_feature.getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula());
    }
    ef += EmpiricalFormula(active_feature.getMetaValue("charge_adducts")); // adducts
    ef -= EmpiricalFormula(String("H") + String(q));
    ef.setCharge(q); // effectively subtract q electrons

    Param p1;
    p1.setValue("statistics:mean", ef.getAverageWeight() / q);
    p1.setValue("interpolation_step", 0.001);
    p1.setValue("isotope:mode:mode", param_.getValue("peak_shape"));
    p1.setValue("intensity_scaling", 0.001); // this removes the problem of to big isotope-model values
    p1.setValue("charge", q);
    double fwhm;
    if (param_.getValue("peak_shape") == "Gaussian")
    {
      fwhm = getPeakWidth_(active_feature.getMZ(), true);
      p1.setValue("isotope:mode:GaussianSD", fwhm);
    }
    else
    {
      fwhm = getPeakWidth_(active_feature.getMZ(), false);
      p1.setValue("isotope:mode:LorentzFWHM", fwhm);
    }

    IsotopeModel* isomodel = new IsotopeModel();
    isomodel->setParameters(p1); // this needs to come BEFORE setSamples() - otherwise the default setSamples() is called here!
    isomodel->setSamples(ef); // this already includes adducts

    if (experiment.size() < 2)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, experiment.size());
    }
    double rt_sampling_rate = experiment[1].getRT() - experiment[0].getRT();
    EGHModel* elutionmodel = new EGHModel();
    chooseElutionProfile_(elutionmodel, active_feature, 1.0, rt_sampling_rate, experiment);
    ProductModel<2> pm;
    pm.setModel(0, elutionmodel); // new'ed models will be deleted by the pm! no need to delete them manually
    pm.setModel(1, isomodel); // new'ed models will be deleted by the pm! no need to delete them manually
    pm.setScale(scale); // scale

    // start and end points of the sampling
    SimTypes::SimCoordinateType rt_start(elutionmodel->getInterpolation().supportMin());
    SimTypes::SimCoordinateType rt_end(elutionmodel->getInterpolation().supportMax());
    if (active_feature.metaValueExists("RT_width_start") && active_feature.metaValueExists("RT_width_end")) // this is a contaminant with sampling restrictions
    {
      rt_start = active_feature.getMetaValue("RT_width_start");
      rt_end = active_feature.getMetaValue("RT_width_end");
    }
    SimTypes::SimCoordinateType mz_start(isomodel->getInterpolation().supportMin());
    SimTypes::SimCoordinateType mz_end(isomodel->getInterpolation().supportMax());

    // add peptide to GLOBAL MS map
    // add CH and new intensity to feature
    samplePeptideModel2D_(pm, mz_start, mz_end, rt_start, rt_end, experiment, experiment_ct, active_feature);
  }

  void RawMSSignalSimulation::samplePeptideModel1D_(const IsotopeModel& pm,
                                                    const SimTypes::SimCoordinateType mz_start,
                                                    const SimTypes::SimCoordinateType mz_end,
                                                    SimTypes::MSSimExperiment& experiment,
                                                    SimTypes::MSSimExperiment& experiment_ct,
                                                    Feature& active_feature)
  {
    SimTypes::SimIntensityType intensity_sum = 0.0;

    //OPENMS_LOG_DEBUG << "Sampling at [mz] " << mz_start << ":" << mz_end << std::endl;

    SimTypes::SimPointType point;

    // centroided GT
    for (IsotopeDistribution::const_iterator iter = pm.getIsotopeDistribution().begin();
         iter != pm.getIsotopeDistribution().end(); ++iter)
    {
      point.setMZ(iter->getMZ());
      point.setIntensity(iter->getIntensity());

      if (point.getIntensity() <= 0.0)
        continue;

      experiment_ct[0].push_back(point);
    }

    std::vector<SimTypes::SimCoordinateType>::const_iterator it_grid = lower_bound(grid_.begin(), grid_.end(), mz_start);
    boost::normal_distribution<double> ndist(mz_error_mean_, mz_error_stddev_);
    for (; it_grid != grid_.end() && (*it_grid) < mz_end; ++it_grid)
    {
      point.setMZ(*it_grid);
      point.setIntensity(pm.getIntensity(DPosition<1>(*it_grid)));

      if (point.getIntensity() <= 0.0)
        continue;

      // add Gaussian distributed m/z error
      double mz_err = ndist(rnd_gen_->getTechnicalRng());
      point.setMZ(fabs(point.getMZ() + mz_err));

      intensity_sum += point.getIntensity();
      experiment[0].push_back(point);
    }
    active_feature.setIntensity(intensity_sum);
  }

  void RawMSSignalSimulation::samplePeptideModel2D_(const ProductModel<2>& pm,
                                                    const SimTypes::SimCoordinateType mz_start,
                                                    const SimTypes::SimCoordinateType mz_end,
                                                    SimTypes::SimCoordinateType rt_start,
                                                    SimTypes::SimCoordinateType rt_end,
                                                    SimTypes::MSSimExperiment& experiment,
                                                    SimTypes::MSSimExperiment& experiment_ct,
                                                    Feature& active_feature)
  {
    if (rt_start <= 0)
      rt_start = 0;

    SimTypes::MSSimExperiment::iterator exp_start = experiment.RTBegin(rt_start);
    SimTypes::MSSimExperiment::iterator exp_ct_start = experiment_ct.RTBegin(rt_start);

    if (exp_start == experiment.end())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0);
    }

    SimTypes::SimIntensityType intensity_sum(0.0);

#ifdef OPENMS_ASSERTIONS
    Int end_scan = std::numeric_limits<Int>::min(); // only used in Debug build
#endif

    IsotopeModel* isomodel = static_cast<IsotopeModel*>(pm.getModel(1));
    IsotopeDistribution iso_dist = isomodel->getIsotopeDistribution();
    SimTypes::SimCoordinateType mz_mono = active_feature.getMZ();
    SimTypes::SimCoordinateType iso_peakdist = isomodel->getParameters().getValue("isotope:distance");
    Int q = active_feature.getCharge();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sample the model ...
    SimTypes::SimCoordinateType rt(0);
    SimTypes::MSSimExperiment::iterator exp_iter = exp_start;
    SimTypes::MSSimExperiment::iterator exp_ct_iter = exp_ct_start;
    for (; rt < rt_end && exp_iter != experiment.end(); ++exp_iter, ++exp_ct_iter)
    {
      rt = exp_iter->getRT();
      double distortion = double(exp_iter->getMetaValue("distortion"));
      double rt_intensity = ((EGHModel*)pm.getModel(0))->getIntensity(rt);

      // centroided GT
      Size iso_pos(0);
      SimTypes::SimPointType point;
      for (IsotopeDistribution::const_iterator iter = iso_dist.begin(); iter != iso_dist.end(); ++iter, ++iso_pos)
      {
        point.setMZ(mz_mono + (iso_pos * iso_peakdist / q));
        point.setIntensity(iter->getIntensity() * rt_intensity * distortion);

        if (point.getIntensity() <= 0.0)
          continue;

        exp_ct_iter->push_back(point);
      }

      // RAW signal (sample it on the grid)
      std::vector<SimTypes::SimCoordinateType>::const_iterator it_grid = lower_bound(grid_.begin(), grid_.end(), mz_start);
      for (; it_grid != grid_.end() && (*it_grid) < mz_end; ++it_grid)
      {
        ProductModel<2>::IntensityType intensity = pm.getIntensity(DPosition<2>(rt, *it_grid)) * distortion;
        if (intensity <= 0.0)
          continue; // intensity cutoff (below that we don't want to see a signal)

        point.setMZ(*it_grid);
        point.setIntensity(intensity);

        //OPENMS_LOG_ERROR << "Sampling " << rt << " , " << mz << " -> " << point.getIntensity() << std::endl;

        // add Gaussian distributed m/z error
#ifdef _OPENMP
        int CURRENT_THREAD = omp_get_thread_num();
        // check if we need to refill the random number pool for this thread
        if (threaded_random_numbers_index_[CURRENT_THREAD] == THREADED_RANDOM_NUMBER_POOL_SIZE_)
        {
          if (mz_error_stddev_ != 0.0)
          {
#pragma omp critical(generate_random_number_for_thread)
            {
              boost::normal_distribution<double> ndist(mz_error_mean_, mz_error_stddev_);
              for (Size i = 0; i < THREADED_RANDOM_NUMBER_POOL_SIZE_; ++i)
              {
                threaded_random_numbers_[CURRENT_THREAD][i] = ndist(rnd_gen_->getTechnicalRng());
              }
            }
          }
          else
          {
            // we do not need to care about concurrency here
            fill(threaded_random_numbers_[CURRENT_THREAD].begin(), threaded_random_numbers_[CURRENT_THREAD].end(), mz_error_mean_);
          }
          // reset index for this thread to first position
          threaded_random_numbers_index_[CURRENT_THREAD] = 0;
        }

        const double mz_err = threaded_random_numbers_[CURRENT_THREAD][threaded_random_numbers_index_[CURRENT_THREAD]++];
#else
        // we can use the normal Gaussian ran-gen if we do not use OPENMP
        boost::normal_distribution<double> ndist(mz_error_mean_, mz_error_stddev_);
        const double mz_err = ndist(rnd_gen_->getTechnicalRng());
#endif
        point.setMZ(std::fabs(point.getMZ() + mz_err));
        exp_iter->push_back(point);

        intensity_sum += point.getIntensity();
      }
      //update last scan affected
#ifdef OPENMS_ASSERTIONS
      end_scan = exp_iter - experiment.begin();
#endif
    }

    OPENMS_POSTCONDITION(end_scan != std::numeric_limits<Int>::min(), "RawMSSignalSimulation::samplePeptideModel2D_(): setting RT bounds failed!");

    // new intensity is AREA==SUM of all peaks
    active_feature.setIntensity(intensity_sum);

    // -------------------------
    // --- store convex hull ---
    // -------------------------
    active_feature.getConvexHulls().clear();

    // use isotope model (to determine mass traces)

    DoubleList isotope_intensities;
    for (IsotopeDistribution::iterator iter = iso_dist.begin();
         iter != iso_dist.end(); ++iter)
    {
      const SimTypes::SimCoordinateType mz = mz_mono + double(iter->getMZ() - iso_dist.begin()->getMZ()) / q; // this is only an approximated trace' m/z position (as we do assume 1Da space between them)

      SimTypes::SimCoordinateType rt_min =  std::numeric_limits<SimTypes::SimCoordinateType>::max();
      SimTypes::SimCoordinateType rt_max = -std::numeric_limits<SimTypes::SimCoordinateType>::max();
      bool has_data = false;

      // for each trace, sample the model again and see how far it extends
      SimTypes::SimCoordinateType rt(0);
      for (exp_iter = exp_start; rt < rt_end && exp_iter != experiment.end(); ++exp_iter)
      {
        rt = exp_iter->getRT();
        double distortion = double(exp_iter->getMetaValue("distortion"));
        ProductModel<2>::IntensityType intensity = pm.getIntensity(DPosition<2>(rt, mz)) * distortion;
        if (intensity <= 0.0)
          continue; // intensity cutoff (below that we don't want to see a signal)

        // update min&max
        if (rt_min > rt)
          rt_min = rt;
        if (rt_max < rt)
          rt_max = rt;
        has_data = true;
      }
      if (!has_data)
        continue;

      // add four edge points of mass trace
      ConvexHull2D hull;
      std::vector<DPosition<2> > points;
      points.push_back(DPosition<2>(rt_min, mz - 0.001));
      points.push_back(DPosition<2>(rt_min, mz + 0.001));
      points.push_back(DPosition<2>(rt_max, mz - 0.001));
      points.push_back(DPosition<2>(rt_max, mz + 0.001));
      hull.addPoints(points);
      active_feature.getConvexHulls().push_back(hull);

      isotope_intensities.push_back(iter->getIntensity());
    }

    active_feature.setMetaValue("isotope_intensities", isotope_intensities);

  }

  void RawMSSignalSimulation::chooseElutionProfile_(EGHModel* const elutionmodel, Feature& feature, const double scale, const double rt_sampling_rate, const SimTypes::MSSimExperiment& experiment)
  {
    SimTypes::SimCoordinateType f_rt = feature.getRT();

    Param p;
    // WARNING: step used to be 'rt_sampling_rate / 3.0', but distortion is not part of RT sim, and thus only
    //          modeled 1:1
    p.setValue("interpolation_step", rt_sampling_rate);
    p.setValue("statistics:variance", 1.0);
    p.setValue("statistics:mean", f_rt);

    p.setValue("egh:height", scale);
    p.setValue("egh:retention", f_rt);

    if (feature.metaValueExists("RT_width_gaussian")) // this is for contaminants only (we want the gaussian distribution width A+B at 5% of maximal height)
    {
      p.setValue("egh:alpha", 0.05);
      p.setValue("egh:A", double(feature.getMetaValue("RT_width_gaussian")) / 2.0 * 0.9); // make width a little smaller as this is only the 5% height cutoff
      p.setValue("egh:B", double(feature.getMetaValue("RT_width_gaussian")) / 2.0 * 0.9);
    }
    else if (feature.metaValueExists("RT_egh_variance") && feature.metaValueExists("RT_egh_tau"))
    {
      // for CE we want wider profiles with higher MT
      double width_factor(1); // default for HPLC
      if (feature.metaValueExists("RT_CE_width_factor"))
        width_factor = feature.getMetaValue("RT_CE_width_factor");

      p.setValue("egh:guess_parameter", "false");
      p.setValue("egh:tau", (double) feature.getMetaValue("RT_egh_tau"));
      p.setValue("egh:sigma_square", ((double) feature.getMetaValue("RT_egh_variance")) * width_factor);
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Elution profile shape cannot be created. Wrong meta-values!", "");
    }

    elutionmodel->setParameters(p); // does the calculation
    //----------------------------------------------------------------------

    SimTypes::SimCoordinateType rt_em_start = elutionmodel->getInterpolation().supportMin();
    SimTypes::SimCoordinateType rt_em_end = elutionmodel->getInterpolation().supportMax();

    // find scan in experiment at which our elution starts
    SimTypes::MSSimExperiment::ConstIterator exp_it = experiment.RTBegin(rt_em_start);
    if (exp_it == experiment.end())
      --exp_it; // we need the last valid RT below, so .end() is not useful

    DoubleList elution_intensities;
    DoubleList elution_bounds;
    elution_bounds.resize(4); // store min and max RT (in seconds and index terms)
    elution_bounds[0] = std::distance(experiment.begin(), exp_it);
    elution_bounds[1] = exp_it->getRT();
    elution_bounds[2] = elution_bounds[0];
    elution_bounds[3] = elution_bounds[1];

    for (; (exp_it != experiment.end()) && (exp_it->getRT() <= rt_em_end); ++exp_it) // .. and disturb values by (an already smoothed) distortion diced in RTSimulation
    {
      double intensity = (double) exp_it->getMetaValue("distortion") * elutionmodel->getInterpolation().value(exp_it->getRT());
      // store elution profile in feature MetaValue
      elution_intensities.push_back(intensity);
      elution_bounds[2] = std::distance(experiment.begin(), exp_it);
      elution_bounds[3] = exp_it->getRT();
    }
    // set elution profile details in feature -> used for MS^E precursor selection in tandem MS later
    feature.setMetaValue("elution_profile_intensities", elution_intensities);
    feature.setMetaValue("elution_profile_bounds", elution_bounds);
  }

  void RawMSSignalSimulation::createContaminants_(SimTypes::FeatureMapSim& c_map, SimTypes::MSSimExperiment& exp, SimTypes::MSSimExperiment& exp_ct)
  {
    if (exp.size() == 1)
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); // not implemented for 1D yet
    }

    if (!contaminants_loaded_)
      loadContaminants();

    IONIZATIONMETHOD this_im = (String)param_.getValue("ionization_type") == "ESI" ? IM_ESI : IM_MALDI;
    c_map.clear(true);

    Size out_of_range_RT(0), out_of_range_MZ(0);
    SimTypes::SimCoordinateType minimal_mz_measurement_limit = exp[0].getInstrumentSettings().getScanWindows()[0].begin;
    SimTypes::SimCoordinateType maximal_mz_measurement_limit = exp[0].getInstrumentSettings().getScanWindows()[0].end;

    for (Size i = 0; i < contaminants_.size(); ++i)
    {
      if (contaminants_[i].im != IM_ALL && contaminants_[i].im != this_im)
        continue;

      if (exp.getMinRT() > contaminants_[i].rt_end || contaminants_[i].rt_start > exp.getMaxRT())
      {
        ++out_of_range_RT;
        continue;
      }
      // ... create contaminants...
      SimTypes::FeatureMapSim::FeatureType feature;
      feature.setRT((contaminants_[i].rt_end + contaminants_[i].rt_start) / 2);
      feature.setMZ((contaminants_[i].sf.getMonoWeight() / contaminants_[i].q) + Constants::PROTON_MASS_U); // m/z (incl. protons)
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
        feature.setMetaValue("RT_width_gaussian", contaminants_[i].rt_end - contaminants_[i].rt_start);
      }
      feature.setMetaValue("sum_formula", contaminants_[i].sf.toString()); // formula without adducts or charges
      feature.setCharge(contaminants_[i].q);
      feature.setMetaValue("charge_adducts", "H" + String(contaminants_[i].q)); // adducts separately
      add2DSignal_(feature, exp, exp_ct);
      c_map.push_back(feature);
    }

    c_map.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
    OPENMS_LOG_INFO << "Contaminants out-of-RT-range: " << out_of_range_RT << " / " << contaminants_.size() << std::endl;
    OPENMS_LOG_INFO << "Contaminants out-of-MZ-range: " << out_of_range_MZ << " / " << contaminants_.size() << std::endl;

  }

  void RawMSSignalSimulation::addShotNoise_(SimTypes::MSSimExperiment& experiment, SimTypes::SimCoordinateType minimal_mz_measurement_limit, SimTypes::SimCoordinateType maximal_mz_measurement_limit)
  {
    // we model the amount of (background) noise as Poisson process
    // i.e. the number of noise data points per unit m/z interval follows a Poisson
    // distribution. Noise intensity is assumed to be exponentially-distributed.
    double rate = param_.getValue("noise:shot:rate");
    double intensity_mean = param_.getValue("noise:shot:intensity-mean");

    // avoid sampling 0 values
    if (rate == 0.0 || intensity_mean == 0.0)
      return;

    const SimTypes::SimCoordinateType window_size = 100.0;
    SimTypes::SimCoordinateType mz_lw = minimal_mz_measurement_limit;
    SimTypes::SimCoordinateType mz_up = window_size + minimal_mz_measurement_limit;

    // we distribute the rate in 100 Th windows
    double scaled_rate = rate * window_size;
    SimTypes::SimPointType shot_noise_peak;

    //distributions to sample from
    boost::random::poisson_distribution<UInt, double> pdist(scaled_rate);
    boost::uniform_real<SimTypes::SimCoordinateType> udist(mz_lw, mz_up);
    boost::random::exponential_distribution<SimTypes::SimCoordinateType> edist(intensity_mean);

    OPENMS_LOG_INFO << "Adding shot noise to spectra ..." << std::endl;
    Size num_intervals = std::ceil((maximal_mz_measurement_limit - minimal_mz_measurement_limit) / window_size);

    for (SimTypes::MSSimExperiment::Iterator spectrum_it = experiment.begin(); spectrum_it != experiment.end(); ++spectrum_it)
    {
      for (Size j = 0; j < num_intervals; ++j)
      {
        UInt counts = pdist(rnd_gen_->getTechnicalRng());
        for (UInt c = 0; c < counts; ++c)
        {
          SimTypes::SimCoordinateType mz = udist(rnd_gen_->getTechnicalRng());
          SimTypes::SimCoordinateType intensity = edist(rnd_gen_->getTechnicalRng());

          // we only add points if they have an intensity>0 and are inside of the measurement range
          if (mz < maximal_mz_measurement_limit)
          {
            shot_noise_peak.setIntensity(intensity);
            shot_noise_peak.setMZ(mz);
            spectrum_it->push_back(shot_noise_peak);
          }
        }

        mz_up += window_size;
      }
    } // end of each scan

    experiment.updateRanges();

  }

  void RawMSSignalSimulation::addBaseLine_(SimTypes::MSSimExperiment& experiment, SimTypes::SimCoordinateType minimal_mz_measurement_limit)
  {
    double scale = param_.getValue("baseline:scaling");
    double shape = param_.getValue("baseline:shape");

    if (scale == 0.0)
      return;

    // TODO: switch to iterator
    for (Size i = 0; i < experiment.size(); ++i)
    {
      for (Size j = 0; j < experiment[i].size(); ++j)
      {
        SimTypes::SimCoordinateType x = (experiment[i][j].getMZ() - minimal_mz_measurement_limit);
        //if (x >= 1000.0) continue; // speed-up TODO: revise this ..

        boost::math::exponential_distribution<double> ed(shape);
        double bx = boost::math::pdf(ed, x);
        bx *= scale;
        experiment[i][j].setIntensity(experiment[i][j].getIntensity() + bx);
      }
    }
  }

  void RawMSSignalSimulation::addWhiteNoise_(SimTypes::MSSimExperiment& experiment)
  {
    OPENMS_LOG_INFO << "Adding white noise to spectra ..." << std::endl;

    // get white noise parameters
    double white_noise_mean = param_.getValue("noise:white:mean");
    double white_noise_stddev = param_.getValue("noise:white:stddev");

    if (white_noise_mean == 0.0 && white_noise_stddev == 0.0)
    {
      return;
    }

    boost::normal_distribution<SimTypes::SimIntensityType> ndist(white_noise_mean, white_noise_stddev);

    for (SimTypes::MSSimExperiment::iterator spectrum_it = experiment.begin(); spectrum_it != experiment.end(); ++spectrum_it)
    {
      SimTypes::MSSimExperiment::SpectrumType new_spec = (*spectrum_it);
      new_spec.clear(false);

      for (SimTypes::MSSimExperiment::SpectrumType::iterator peak_it = (*spectrum_it).begin(); peak_it != (*spectrum_it).end(); ++peak_it)
      {
        SimTypes::SimIntensityType intensity = peak_it->getIntensity() + ndist(rnd_gen_->getTechnicalRng());
        if (intensity > 0.0)
        {
          peak_it->setIntensity(intensity);
          new_spec.push_back(*peak_it);
        }
      }

      *spectrum_it = new_spec;
    }
  }

  void RawMSSignalSimulation::addDetectorNoise_(SimTypes::MSSimExperiment& experiment)
  {
    OPENMS_LOG_INFO << "Adding detector noise to spectra ..." << std::endl;

    // get white noise parameters
    double detector_noise_mean = param_.getValue("noise:detector:mean");
    double detector_noise_stddev = param_.getValue("noise:detector:stddev");

    if (detector_noise_mean == 0.0 && detector_noise_stddev == 0.0)
    {
      OPENMS_LOG_INFO << "Detector noise was disabled." << std::endl;
      return;
    }

    boost::normal_distribution<SimTypes::SimIntensityType> ndist(detector_noise_mean, detector_noise_stddev);
    for (SimTypes::MSSimExperiment::iterator spectrum_it = experiment.begin(); spectrum_it != experiment.end(); ++spectrum_it)
    {
      SimTypes::MSSimExperiment::SpectrumType new_spec = (*spectrum_it);
      new_spec.clear(false);

      std::vector<SimTypes::SimCoordinateType>::iterator grid_it = grid_.begin();
      SimTypes::MSSimExperiment::SpectrumType::iterator peak_it = spectrum_it->begin();
      for (; grid_it != grid_.end(); ++grid_it)
      {
        // if peak is in grid
        if (peak_it != spectrum_it->end() && *grid_it == peak_it->getMZ())
        {
          SimTypes::SimIntensityType intensity = peak_it->getIntensity() + ndist(rnd_gen_->getTechnicalRng());
          if (intensity > 0.0)
          {
            peak_it->setIntensity(intensity);
            new_spec.push_back(*peak_it);
          }
          ++peak_it;
        }
        else // we have no point here, generate one if noise is above 0
        {
          SimTypes::SimIntensityType intensity = ndist(rnd_gen_->getTechnicalRng());
          if (intensity > 0.0)
          {
            SimTypes::MSSimExperiment::SpectrumType::PeakType noise_peak;
            noise_peak.setMZ(*grid_it);
            noise_peak.setIntensity(intensity);
            new_spec.push_back(noise_peak);
          }
        }
      }

      *spectrum_it = new_spec;
    }

  }

  void RawMSSignalSimulation::getSamplingGrid_(std::vector<SimTypes::SimCoordinateType>& grid, const SimTypes::SimCoordinateType mz_min, const SimTypes::SimCoordinateType mz_max, const Int step_Da)
  {
    if (fabs(mz_max - mz_min) < step_Da)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sampling grid seems very small. This cannot be computed!");
    }
    grid.clear();
    SimTypes::SimCoordinateType mz = mz_min;
    double sampling_rate(0);
    while (mz <= mz_max)
    {
      SimTypes::SimCoordinateType fwhm = getPeakWidth_(mz, param_.getValue("peak_shape") == "Gaussian");
      sampling_rate = (fwhm / sampling_points_per_FWHM_);
      SimTypes::SimCoordinateType mz_cell_end = std::min(mz + step_Da, mz_max);
      for (; mz <= mz_cell_end; mz += sampling_rate)
      {
        grid.push_back(mz);
      }
    }
    grid.push_back(mz + sampling_rate); // one more point after mz_max, for binary search later
    return;
  }

  // TODO: add instrument specific sampling technique
  void RawMSSignalSimulation::compressSignals_(SimTypes::MSSimExperiment& experiment)
  {
    if (experiment.size() < 1 || experiment[0].getInstrumentSettings().getScanWindows().size() < 1)
    {
      throw Exception::IllegalSelfOperation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    SimTypes::SimCoordinateType min_mz = experiment[0].getInstrumentSettings().getScanWindows()[0].begin;
    SimTypes::SimCoordinateType max_mz = experiment[0].getInstrumentSettings().getScanWindows()[0].end;

    if (min_mz >= max_mz)
    {
      OPENMS_LOG_WARN << "No data to compress." << std::endl;
      return;
    }

    typedef std::vector<SimTypes::SimCoordinateType>::iterator GridTypeIt;
    std::vector<SimTypes::SimCoordinateType> grid;
    getSamplingGrid_(grid, min_mz, max_mz, 5); // every 5 Da we adjust the sampling width by local FWHM

    if (grid.size() < 3)
    {
      OPENMS_LOG_WARN << "Data spacing is weird - either you selected a very small interval or a very low resolution - or both. Not compressing." << std::endl;
      return;
    }

    Size point_count_before(0), point_count_after(0);
    SimTypes::SimPointType p;
    for (Size i = 0; i < experiment.size(); ++i)
    {
      if (experiment[i].size() <= 1)
        continue;

      if (experiment[i].isSorted() == false) // this should be true - however we check
      {
        experiment[i].sortByPosition();
      }

      // copy Spectrum and remove Peaks ..
      SimTypes::MSSimExperiment::SpectrumType cont = experiment[i];
      cont.clear(false);

      GridTypeIt grid_pos = grid.begin();
      GridTypeIt grid_pos_next(grid_pos + 1);

      double int_sum(0);
      bool break_scan(false);
      // match points to closest grid point
      for (Size j = 0; j < experiment[i].size(); ++j)
      {
        Size advance_by_binary_search = 3;
        while (fabs((*grid_pos_next) - experiment[i][j].getMZ()) < fabs((*grid_pos) - experiment[i][j].getMZ()))
        {
          if (int_sum > 0) // we collected some points before --> save them
          {
            p.setIntensity(int_sum);
            p.setMZ(*grid_pos);
            cont.push_back(p);
            int_sum = 0; // reset
          }

          if (--advance_by_binary_search == 0)
          {
            // advance using binary search
            grid_pos_next = std::lower_bound(grid_pos, grid.end(), experiment[i][j].getMZ());
            grid_pos = grid_pos_next - 1; // this should always work, since we ran at least 3 steps forward before
            advance_by_binary_search = 10; // just so we do not run into here again
          }
          else
          {
            // advance to next grid element
            ++grid_pos;
            ++grid_pos_next;
          }

          if (grid_pos_next == grid.end())
          {
            break_scan = true;
            break;
          }
        }
        if (break_scan)
          break; // skip remaining points of the scan (we reached the end of the grid)

        int_sum += experiment[i][j].getIntensity();

      } // end of scan

      if (int_sum > 0) // don't forget the last one
      {
        p.setIntensity(int_sum);
        p.setMZ(*grid_pos);
        cont.push_back(p);
      }

      point_count_before += experiment[i].size(); // stats
      experiment[i] = cont;
      point_count_after += experiment[i].size();
    }

    if (point_count_before != 0)
    {
      OPENMS_LOG_INFO << "Compressed data to grid ... " <<  point_count_before << " --> " << point_count_after << " (" << (point_count_after * 100 / point_count_before) << "%)\n";
    }
    else
    {
      OPENMS_LOG_INFO << "Not enough points in map .. did not compress!\n";
    }

    return;
  }

  SimTypes::SimIntensityType RawMSSignalSimulation::getFeatureScaledIntensity_(const SimTypes::SimIntensityType feature_intensity, const SimTypes::SimIntensityType natural_scaling_factor)
  {
    SimTypes::SimIntensityType intensity = feature_intensity * natural_scaling_factor * intensity_scale_;

    // add some noise
    // TODO: German comment
    // TODO: variables model f??r den intensit??ts-einfluss
    // e.g. sqrt(intensity) || ln(intensity)
    boost::normal_distribution<SimTypes::SimIntensityType> ndist(0, intensity_scale_stddev_ * intensity);
    intensity += ndist(rnd_gen_->getTechnicalRng());

    return intensity;
  }

}
