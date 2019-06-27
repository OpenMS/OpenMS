// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow, Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>
#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>

namespace OpenMS
{
  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation() :
    DefaultParamHandler("RawTandemMSSignalSimulation"),
    rnd_gen_(new SimTypes::SimRandomNumberGenerator())
  {
    initParam_();
  }

  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr rng) :
    DefaultParamHandler("RawTandemMSSignalSimulation"),
    rnd_gen_(rng)
  {
    initParam_();
  }

  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(const RawTandemMSSignalSimulation& source) :
    DefaultParamHandler(source)
  {
    setParameters(source.getParameters());
    rnd_gen_ = source.rnd_gen_;
  }

  RawTandemMSSignalSimulation& RawTandemMSSignalSimulation::operator=(const RawTandemMSSignalSimulation& source)
  {
    DefaultParamHandler::operator=(source);
    setParameters(source.getParameters());
    rnd_gen_ = source.rnd_gen_;
    return *this;
  }

  void RawTandemMSSignalSimulation::initParam_()
  {
    // Tandem MS params
    defaults_.setValue("status", "disabled", "Create Tandem-MS scans?");
    defaults_.setValidStrings("status", ListUtils::create<String>("disabled,precursor,MS^E"));

    subsections_.push_back("Precursor:");
    defaults_.insert("Precursor:", OfflinePrecursorIonSelection().getDefaults());
    defaults_.remove("Precursor:peptides_per_protein");
    defaults_.setValue("Precursor:charge_filter", ListUtils::create<Int>("2,3"), "Charges considered for MS2 fragmentation.");
    defaults_.setMinInt("Precursor:charge_filter", 1);
    defaults_.setMaxInt("Precursor:charge_filter", 5);

    defaults_.setValue("MS_E:add_single_spectra", "false", "If true, the MS2 spectra for each peptide signal are included in the output (might be a lot). They will have a meta value 'MSE_DebugSpectrum' attached, so they can be filtered out. Native MS_E spectra will have 'MSE_Spectrum' instead.");
    defaults_.setValidStrings("MS_E:add_single_spectra", ListUtils::create<String>("true,false"));
    defaults_.setValue("tandem_mode", 0, "Algorithm to generate the tandem-MS spectra. 0 - fixed intensities, 1 - SVC prediction (abundant/missing), 2 - SVR prediction of peak intensity \n");
    defaults_.setMinInt("tandem_mode", 0);
    defaults_.setMaxInt("tandem_mode", 2);
    defaults_.setValue("svm_model_set_file", "examples/simulation/SvmModelSet.model", "File containing the filenames of SVM Models for different charge variants");

    subsections_.push_back("TandemSim:");
    defaults_.insert("TandemSim:Simple:", TheoreticalSpectrumGenerator().getDefaults());
    Param svm_par = SvmTheoreticalSpectrumGenerator().getDefaults();
    svm_par.remove("svm_mode");
    svm_par.remove("model_file_name");
    defaults_.insert("TandemSim:SVM:", svm_par);

    // sync'ed Param (also appears in IonizationSimulation)
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    defaults_.setValidStrings("ionization_type", ListUtils::create<String>("MALDI,ESI"));

    defaultsToParam_();
  }

  RawTandemMSSignalSimulation::~RawTandemMSSignalSimulation()
  {
  }

  void RawTandemMSSignalSimulation::generateMSESpectra_(const SimTypes::FeatureMapSim& features, const SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& ms2)
  {
    //get tandem mode
    Size tandem_mode = param_.getValue("tandem_mode");

    Param simple_gen_params = param_.copy("TandemSim:Simple:", true);

    TheoreticalSpectrumGenerator simple_generator;
    simple_generator.setParameters(simple_gen_params);

    SvmTheoreticalSpectrumGeneratorSet svm_spec_gen_set;
    //this set will hold the precursor charges that have an Svm model
    std::set<Size> svm_model_charges;

    // if SVR or SVC shall be used
    if (tandem_mode)
    {
      String svm_filename = param_.getValue("svm_model_set_file");
      if (!File::readable(svm_filename)) // look in OPENMS_DATA_PATH
      {
        svm_filename = File::find(svm_filename);
      }

      svm_spec_gen_set.load(svm_filename);
      svm_spec_gen_set.getSupportedCharges(svm_model_charges);

      //set the parameters for each model
      Param svm_gen_params = param_.copy("TandemSim:SVM:", true);
      svm_gen_params.setValue("svm_mode", tandem_mode - 1);
      std::set<Size>::iterator it;
      for (it = svm_model_charges.begin(); it != svm_model_charges.end(); ++it)
      {
        svm_spec_gen_set.getSvmModel(*it).setParameters(svm_gen_params);
      }
    }

    Param p;
    p.setValue("block_method:rt_block_size", features.size()); // merge all single spectra
    p.setValue("block_method:ms_levels", ListUtils::create<Int>("2"));
    SpectraMerger sm;
    sm.setParameters(p);

    double sampling_rate = 1;
    //guess sampling rate from two adjacent full scans:
    if (experiment.size() >= 2)
      sampling_rate = experiment[1].getRT() - experiment[0].getRT();

    SimTypes::MSSimExperiment precomputed_MS2;
    precomputed_MS2.resize(features.size());


    // preparation & validation of input
    for (Size i_f = 0; i_f < features.size(); ++i_f)
    {
      // sample MS2 spectra for each feature
      AASequence seq = features[i_f].getPeptideIdentifications()[0].getHits()[0].getSequence();
      //TODO: work around RichPeak1D restriction
      PeakSpectrum tmp_spec;
      Int prec_charge = features[i_f].getCharge();

      if (tandem_mode && svm_model_charges.count(prec_charge))
      {
        svm_spec_gen_set.simulate(tmp_spec, seq, rnd_gen_->getBiologicalRng(), prec_charge);
      }
      else
      {
        simple_generator.getSpectrum(tmp_spec, seq, 1, prec_charge);
      }
      for (Size peak = 0; peak < tmp_spec.size(); ++peak)
      {
        Peak1D p = tmp_spec[peak];
        precomputed_MS2[i_f].push_back(p);
      }

      precomputed_MS2[i_f].setMSLevel(2);
      Precursor prec;
      prec.setMZ(features[i_f].getMZ());
      precomputed_MS2[i_f].setPrecursors(std::vector<Precursor>(1, prec));
      precomputed_MS2[i_f].setMetaValue("FeatureID", static_cast<String>(features[i_f].getUniqueId()));


      // validate features meta values exist and are valid:
      if (!features[i_f].metaValueExists("elution_profile_bounds")
         ||
          !features[i_f].metaValueExists("elution_profile_intensities"))
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MetaValue:elution_profile_***");
      }
      // check if values fit the experiment:

#ifdef OPENMS_ASSERTIONS
      const DoubleList& elution_bounds = features[i_f].getMetaValue("elution_profile_bounds");
      const DoubleList& elution_ints   = features[i_f].getMetaValue("elution_profile_intensities");

      OPENMS_PRECONDITION(elution_bounds[0] < experiment.size(), "Elution profile out of bounds (left)");
      OPENMS_PRECONDITION(elution_bounds[2] < experiment.size(), "Elution profile out of bounds (right)");
      OPENMS_PRECONDITION(experiment[elution_bounds[0]].getRT() == elution_bounds[1], "Elution profile RT shifted (left)");
      OPENMS_PRECONDITION(experiment[elution_bounds[2]].getRT() == elution_bounds[3], "Elution profile RT shifted (right)");
      OPENMS_PRECONDITION(elution_bounds[2] - elution_bounds[0] + 1 == elution_ints.size(), "Elution profile size does not match bounds");
#endif
    }

    // creating the MS^E scan:
    bool add_debug_spectra = static_cast<String>(param_.getValue("MS_E:add_single_spectra")) == "true";

    for (Size i = 0; i < experiment.size(); ++i) // create MS2 for every MS scan
    { // check which features elute in the current MS scan
      std::vector<Size> features_fragmented;
      for (Size i_f = 0; i_f < features.size(); ++i_f)
      {
        const DoubleList& elution_bounds = features[i_f].getMetaValue("elution_profile_bounds");
        if ((elution_bounds[1] <= experiment[i].getRT()) && (experiment[i].getRT() <= elution_bounds[3]))
        {
          features_fragmented.push_back(i_f);
        }
      }

      if (features_fragmented.empty())
        continue;

      // now we have all features that elute in this scan -> create MS2 scans
      SimTypes::MSSimExperiment MS2_spectra;
      MS2_spectra.resize(features_fragmented.size());

      StringList feature_seq;
      IntList feature_intensities;
      for (Size index = 0; index < features_fragmented.size(); ++index)
      {
        Size i_f = features_fragmented[index];
        // create spectrum
        MS2_spectra[index] = precomputed_MS2[i_f];
        MS2_spectra[index].setRT(experiment[i].getRT() + sampling_rate * (double(index + 1) / double(features_fragmented.size() + 2)));
        // adjust intensity of single MS2 spectra by feature intensity
        const DoubleList& elution_bounds = features[i_f].getMetaValue("elution_profile_bounds");
        const DoubleList& elution_ints   = features[i_f].getMetaValue("elution_profile_intensities");
        double factor = elution_ints[i - elution_bounds[0]] * features[i_f].getIntensity();
        for (SimTypes::MSSimExperiment::SpectrumType::iterator it = MS2_spectra[index].begin(); it != MS2_spectra[index].end(); ++it)
        {
          it->setIntensity(it->getIntensity() * factor);
        }
        // store true sequence of spectrum
        feature_seq.push_back(features[i_f].getPeptideIdentifications()[0].getHits()[0].getSequence().toString() + ":" + features[i_f].getCharge());
        MS2_spectra[index].setMetaValue("sequence", feature_seq.back());
        feature_intensities.push_back(factor);
        MS2_spectra[index].setMetaValue("intensity factor", factor);
      }

      // debug: also add single spectra
      if (add_debug_spectra)
      {
        for (Size ii = 0; ii < MS2_spectra.size(); ++ii)
        {
          ms2.addSpectrum(MS2_spectra[ii]); // DEBUG
          ms2.getSpectra().back().setMetaValue("MSE_DebugSpectrum", "true");
        }
      }

      // merge all MS2 spectra
      sm.mergeSpectraBlockWise(MS2_spectra);
      if (MS2_spectra.size() != 1)
        throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, MS2_spectra.size());
      // store merged spectrum
      MS2_spectra[0].setMetaValue("MSE_Spectrum", "true");
      MS2_spectra[0].setMetaValue("MSE_sequences", feature_seq);
      MS2_spectra[0].setMetaValue("MSE_intensities", feature_intensities);
      ms2.addSpectrum(MS2_spectra[0]);

    }

  }

  void RawTandemMSSignalSimulation::generatePrecursorSpectra_(const SimTypes::FeatureMapSim& features, const SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& ms2)
  {
    IntList qs = param_.getValue("Precursor:charge_filter");
    std::set<Int> qs_set(qs.begin(), qs.end());

    //** precursor selection **//
    OfflinePrecursorIonSelection ps;
    Param param = param_.copy("Precursor:", true);
    param.remove("charge_filter");
    ps.setParameters(param);
    // different selection strategies for MALDI and ESI
    bool is_MALDI = static_cast<String>(param_.getValue("ionization_type")) == "MALDI";

    // fill 'ms2' with precursor information, but leave spectra empty
    ps.makePrecursorSelectionForKnownLCMSMap(features, experiment, ms2, qs_set, is_MALDI);


    //** actual MS2 signal **//
    std::cout << "MS2 features selected: " << ms2.size() << "\n";

    TheoreticalSpectrumGenerator simple_generator;
    Param simple_gen_params = param_.copy("TandemSim:Simple:", true);
    simple_generator.setParameters(simple_gen_params);

    //get tandem mode
    Size tandem_mode = param_.getValue("tandem_mode");

    SvmTheoreticalSpectrumGeneratorSet svm_spec_gen_set;
    //this set will hold the precursor charges that have an SVM model
    std::set<Size> svm_model_charges;

    // if SVR or SVC shall be used
    if (tandem_mode)
    {
      String svm_filename = param_.getValue("svm_model_set_file");
      if (!File::readable(svm_filename)) // look in OPENMS_DATA_PATH
      {
        svm_filename = File::find(svm_filename);
      }

      svm_spec_gen_set.load(svm_filename);
      svm_spec_gen_set.getSupportedCharges(svm_model_charges);

      //set the parameters for each model
      Param svm_gen_params = param_.copy("TandemSim:SVM:", true);
      svm_gen_params.setValue("svm_mode", tandem_mode - 1);
      std::set<Size>::iterator it;
      for (it = svm_model_charges.begin(); it != svm_model_charges.end(); ++it)
      {
        svm_spec_gen_set.getSvmModel(*it).setParameters(svm_gen_params);
      }
    }

    for (Size i = 0; i < ms2.size(); ++i)
    {
      IntList ids = ms2[i].getMetaValue("parent_feature_ids");

      // std::cerr << "precursor: " << i << " (#ids: " << ids.size() << ")\n";
      // std::cerr << ms2[i].getRT() << " " << ms2[i].getPrecursors()[0].getMZ() << "\n";

      OPENMS_POSTCONDITION(ids.size() == ms2[i].getPrecursors().size(), "#parent features should be equal to # of precursors")

      SimTypes::MSSimExperiment tmp_spectra;
      tmp_spectra.resize(ids.size());

      for (Size id = 0; id < ids.size(); ++id)
      {
        double prec_intens = ms2[i].getPrecursors()[id].getIntensity();
        AASequence seq = features[ids[id]].getPeptideIdentifications()[0].getHits()[0].getSequence();
        PeakSpectrum tmp_spec;

        Int prec_charge = features[ids[id]].getCharge();

        if (tandem_mode && svm_model_charges.count(prec_charge))
        {
          svm_spec_gen_set.simulate(tmp_spec, seq, rnd_gen_->getBiologicalRng(), prec_charge);
        }
        else
        {
          simple_generator.getSpectrum(tmp_spec, seq, 1, prec_charge);
        }

        // scale intensity and copy 
        for (Size peak = 0; peak < tmp_spec.size(); ++peak)
        {
          Peak1D p = tmp_spec[peak];
          p.setIntensity(p.getIntensity() * prec_intens);
          tmp_spectra[id].push_back(p);
        }
        tmp_spectra[id].sortByPosition();
      }

      if (tmp_spectra.size() > 1)
      {
        Param p;
        p.setValue("block_method:rt_block_size", ids.size()); // merge all single spectra from (multiple) precursors
        p.setValue("block_method:ms_levels", ListUtils::create<Int>("2"));
        SpectraMerger sm;
        sm.setParameters(p);
        sm.mergeSpectraBlockWise(tmp_spectra);
      }
      // preserve precursor information etc and just insert peaks
      ms2[i].insert(ms2[i].begin(), tmp_spectra[0].begin(), tmp_spectra[0].end());
    }
  }

  void RawTandemMSSignalSimulation::generateRawTandemSignals(const SimTypes::FeatureMapSim& features, SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& experiment_ct)
  {
    OPENMS_LOG_INFO << "Tandem MS Simulation ... ";

    // will hold the MS2 scans
    SimTypes::MSSimExperiment ms2;

    if (param_.getValue("status") == "disabled")
    {
      OPENMS_LOG_INFO << "disabled" << std::endl;
      return;
    }
    else if (param_.getValue("status") == "precursor")
    {
      OPENMS_LOG_INFO << "precursor" << std::endl;
      generatePrecursorSpectra_(features, experiment, ms2);
    }
    else // MS^E
    {
      OPENMS_LOG_INFO << "MS^E" << std::endl;
      generateMSESpectra_(features, experiment, ms2);
    }

    // append MS2 to experiment
    experiment.getSpectra().insert(experiment.end(), ms2.begin(), ms2.end());
    // .. and to picked experiment (if we ever simulate MS2 data with raw peaks, this needs to be adapted)
    experiment_ct.getSpectra().insert(experiment_ct.end(), ms2.begin(), ms2.end());

  }

}
