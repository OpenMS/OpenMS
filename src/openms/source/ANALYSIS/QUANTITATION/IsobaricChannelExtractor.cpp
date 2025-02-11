// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>

#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/MATH/StatisticFunctions.h>

// #define ISOBARIC_CHANNEL_EXTRACTOR_DEBUG
// #undef ISOBARIC_CHANNEL_EXTRACTOR_DEBUG

namespace OpenMS
{

  // Maximum allowed search window for TMT-10 reporter ions. The channels are only 0.006 Th apart.
  // Allowing anything larger will result in wrong quantifications for empty channels.
  // Also used for TMT_11PLEX
  double TMT_10AND11PLEX_CHANNEL_TOLERANCE = 0.003;

  /// small quality control class, holding temporary data for reporting
  struct ChannelQC
  {
    // C'tor
    ChannelQC() :
      mz_deltas()
      
    {}

    std::vector<double> mz_deltas; ///< m/z distance between expected and observed reporter ion closest to expected position
    int signal_not_unique{0};  ///< counts if more than one peak was found within the search window of each reporter position
  };


  IsobaricChannelExtractor::PuritySate_::PuritySate_(const PeakMap& targetExp) :
    baseExperiment(targetExp)
  {
    // initialize precursorScan with end(), it will be updated later on
    // from the calling method
    precursorScan = baseExperiment.end();

    // find the first ms1 scan in the experiment
    followUpScan = baseExperiment.begin();
    while (followUpScan != baseExperiment.end() && followUpScan->getMSLevel() != 1)
    {
      ++followUpScan;
    }

    // check if we found one
    hasFollowUpScan = followUpScan != baseExperiment.end();
  }


  void IsobaricChannelExtractor::PuritySate_::advanceFollowUp(const double rt)
  {
    // advance follow up scan until we found a ms1 scan with a bigger RT
    if (followUpScan != baseExperiment.end()) ++followUpScan;
    while (followUpScan != baseExperiment.end())
    {
      if (followUpScan->getMSLevel() == 1 && followUpScan->getRT() > rt)
      {
        break;
      }
      ++followUpScan;
    }

    // check if we found one
    hasFollowUpScan = followUpScan != baseExperiment.end();
  }

  bool IsobaricChannelExtractor::PuritySate_::followUpValid(const double rt) const
  {
    return hasFollowUpScan ? rt < followUpScan->getRT() : true;
  }

  IsobaricChannelExtractor::IsobaricChannelExtractor(const IsobaricQuantitationMethod* const quant_method) :
    DefaultParamHandler("IsobaricChannelExtractor"),
    quant_method_(quant_method),
    selected_activation_("any"),
    reporter_mass_shift_(0.1),
    min_precursor_intensity_(1.0),
    keep_unannotated_precursor_(true),
    min_reporter_intensity_(0.0),
    remove_low_intensity_quantifications_(false),
    min_precursor_purity_(0.0),
    max_precursor_isotope_deviation_(10),
    interpolate_precursor_purity_(false)
  {
    setDefaultParams_();
  }

  IsobaricChannelExtractor::IsobaricChannelExtractor(const IsobaricChannelExtractor& other) = default;

  IsobaricChannelExtractor& IsobaricChannelExtractor::operator=(const IsobaricChannelExtractor& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    quant_method_ = rhs.quant_method_;
    selected_activation_ = rhs.selected_activation_;
    reporter_mass_shift_ = rhs.reporter_mass_shift_;
    min_precursor_intensity_ = rhs.min_precursor_intensity_;
    keep_unannotated_precursor_ = rhs.keep_unannotated_precursor_;
    min_reporter_intensity_ = rhs.min_reporter_intensity_;
    remove_low_intensity_quantifications_ = rhs.remove_low_intensity_quantifications_;
    min_precursor_purity_ = rhs.min_precursor_purity_;
    max_precursor_isotope_deviation_ = rhs.max_precursor_isotope_deviation_;
    interpolate_precursor_purity_ = rhs.interpolate_precursor_purity_;

    return *this;
  }

  void IsobaricChannelExtractor::setDefaultParams_()
  {
    defaults_.setValue("select_activation", "auto", "Operate only on MSn scans where any of its precursors features a certain activation method. Setting to \"auto\" uses HCD and HCID spectra. Set to empty string if you want to disable filtering.");
    std::vector<std::string> activation_list;
    activation_list.emplace_back("auto");
    activation_list.insert(activation_list.end(), Precursor::NamesOfActivationMethod, Precursor::NamesOfActivationMethod + Precursor::SIZE_OF_ACTIVATIONMETHOD - 1);
    activation_list.emplace_back("any"); // allow disabling this

    defaults_.setValidStrings("select_activation", activation_list);

    defaults_.setValue("reporter_mass_shift", 0.002, "Allowed shift (left to right) in Th from the expected position.");
    defaults_.setMinFloat("reporter_mass_shift", 0.0001); // ~0.7ppm -- no need to allow any lower value; this is more than enough for TMT-10plex (0.006 distance between channels, i.e. 60 times wider)
    defaults_.setMaxFloat("reporter_mass_shift", 0.5);

    defaults_.setValue("min_precursor_intensity", 1.0, "Minimum intensity of the precursor to be extracted. MS/MS scans having a precursor with a lower intensity will not be considered for quantitation.");
    defaults_.setMinFloat("min_precursor_intensity", 0.0);

    defaults_.setValue("keep_unannotated_precursor", "true", "Flag if precursor with missing intensity value or missing precursor spectrum should be included or not.");
    defaults_.setValidStrings("keep_unannotated_precursor", {"true","false"});

    defaults_.setValue("min_reporter_intensity", 0.0, "Minimum intensity of the individual reporter ions to be extracted.");
    defaults_.setMinFloat("min_reporter_intensity", 0.0);

    defaults_.setValue("discard_low_intensity_quantifications", "false", "Remove all reporter intensities if a single reporter is below the threshold given in 'min_reporter_intensity'.");
    defaults_.setValidStrings("discard_low_intensity_quantifications", {"true","false"});

    defaults_.setValue("min_precursor_purity", 0.0, "Minimum fraction of the total intensity in the isolation window of the precursor spectrum attributable to the selected precursor.");
    defaults_.setMinFloat("min_precursor_purity", 0.0);
    defaults_.setMaxFloat("min_precursor_purity", 1.0);

    defaults_.setValue("precursor_isotope_deviation", 10.0, "Maximum allowed deviation (in ppm) between theoretical and observed isotopic peaks of the precursor peak in the isolation window to be counted as part of the precursor.");
    defaults_.setMinFloat("precursor_isotope_deviation", 0.0);
    defaults_.addTag("precursor_isotope_deviation", "advanced");

    defaults_.setValue("purity_interpolation", "true", "If set to true the algorithm will try to compute the purity as a time weighted linear combination of the precursor scan and the following scan. If set to false, only the precursor scan will be used.");
    defaults_.setValidStrings("purity_interpolation", {"true","false"});
    defaults_.addTag("purity_interpolation", "advanced");

    defaultsToParam_();
  }

  void IsobaricChannelExtractor::updateMembers_()
  {
    selected_activation_ = getParameters().getValue("select_activation").toString();
    reporter_mass_shift_ = getParameters().getValue("reporter_mass_shift");
    min_precursor_intensity_ = getParameters().getValue("min_precursor_intensity");
    keep_unannotated_precursor_ = getParameters().getValue("keep_unannotated_precursor") == "true";
    min_reporter_intensity_ = getParameters().getValue("min_reporter_intensity");
    remove_low_intensity_quantifications_ = getParameters().getValue("discard_low_intensity_quantifications") == "true";
    min_precursor_purity_ = getParameters().getValue("min_precursor_purity");
    max_precursor_isotope_deviation_ = getParameters().getValue("precursor_isotope_deviation");
    interpolate_precursor_purity_ = getParameters().getValue("purity_interpolation") == "true";
    Size number_of_channels = quant_method_->getNumberOfChannels();

    /* check for sensible parameters */
    if ((( number_of_channels == 10) || (number_of_channels == 11))
        && reporter_mass_shift_ > TMT_10AND11PLEX_CHANNEL_TOLERANCE)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Both TMT-10plex and TMT-11plex require reporter mass shifts <= 0.003 to avoid channel ambiguity!");
    }
  }

  bool IsobaricChannelExtractor::isValidPrecursor_(const Precursor& precursor) const
  {
    return (!(precursor.getIntensity() > 0.0) && keep_unannotated_precursor_) || !(precursor.getIntensity() < min_precursor_intensity_);
  }

  bool IsobaricChannelExtractor::hasLowIntensityReporter_(const ConsensusFeature& cf) const
  {
    for (ConsensusFeature::const_iterator cf_it = cf.begin();
         cf_it != cf.end();
         ++cf_it)
    {
      if (cf_it->getIntensity() == 0.0)
      {
        return true;
      }
    }

    return false;
  }

  double IsobaricChannelExtractor::computeSingleScanPrecursorPurity_(const PeakMap::ConstIterator& ms2_spec, const PeakMap::SpectrumType& precursor_spec) const
  {

    typedef PeakMap::SpectrumType::ConstIterator const_spec_iterator;

    // compute distance between isotopic peaks based on the precursor charge.
    const double charge_dist = Constants::NEUTRON_MASS_U / static_cast<double>(ms2_spec->getPrecursors()[0].getCharge());

    // the actual boundary values
    const double strict_lower_mz = ms2_spec->getPrecursors()[0].getMZ() - ms2_spec->getPrecursors()[0].getIsolationWindowLowerOffset();
    const double strict_upper_mz = ms2_spec->getPrecursors()[0].getMZ() + ms2_spec->getPrecursors()[0].getIsolationWindowUpperOffset();

    const double fuzzy_lower_mz = strict_lower_mz - (strict_lower_mz * max_precursor_isotope_deviation_ / 1000000);
    const double fuzzy_upper_mz = strict_upper_mz + (strict_upper_mz * max_precursor_isotope_deviation_ / 1000000);

    // first find the actual precursor peak
    Size precursor_peak_idx = precursor_spec.findNearest(ms2_spec->getPrecursors()[0].getMZ());
    const Peak1D& precursor_peak = precursor_spec[precursor_peak_idx];

    // now we get ourselves some border iterators
    const_spec_iterator lower_bound = precursor_spec.MZBegin(fuzzy_lower_mz);
    const_spec_iterator upper_bound = precursor_spec.MZEnd(ms2_spec->getPrecursors()[0].getMZ());

    Peak1D::IntensityType precursor_intensity = precursor_peak.getIntensity();
    Peak1D::IntensityType total_intensity = precursor_peak.getIntensity();

    // ------------------------------------------------------------------------------
    // try to find a match for our isotopic peak on the left side

    double expected_next_mz = precursor_peak.getMZ() - charge_dist;

    while (expected_next_mz > fuzzy_lower_mz)
    {
      // find nearest peak in precursor window
      const_spec_iterator np_it = precursor_spec.MZBegin(lower_bound, expected_next_mz, upper_bound);

      // handle border cases

      // check if next peak has smaller dist
      const_spec_iterator np_it2 = np_it;
      ++np_it;

      if (std::fabs(np_it2->getMZ() - expected_next_mz) < std::fabs(np_it->getMZ() - expected_next_mz))
      {
        np_it = np_it2;
      }

      // compute difference between found peak and expected
      double min_diff = std::fabs(np_it->getMZ() - expected_next_mz)  * 1000000 / expected_next_mz;

      // check if we found an isotopic peak
      if (min_diff < max_precursor_isotope_deviation_)
      {
#ifdef ISOBARIC_CHANNEL_EXTRACTOR_DEBUG
        std::cerr << "Mark peak as isotopic peak POS: " << precursor_spec[min_idx] << " (diff: " << min_diff << " vs " << max_precursor_isotope_deviation_ << ")" << std::endl;
#endif
        if (np_it->getMZ() > strict_lower_mz)
        {
          precursor_intensity += np_it->getIntensity();
        }
        else
        {
          // we're in the fuzzy area, so we will take only 50% of the given intensity
          // since we assume that the isolation window borders are not sharp
          precursor_intensity += 0.5 * np_it->getIntensity();
        }

        // update expected_next_mz
        expected_next_mz = np_it->getMZ() - charge_dist;
      }
      else
      {
#ifdef ISOBARIC_CHANNEL_EXTRACTOR_DEBUG
        std::cerr << "No matching isotopic peak for expected pos: " << expected_next_mz << " (min reached diff: " << min_diff << " vs " << max_precursor_isotope_deviation_ << ")" << std::endl;
#endif
        // update expected_next_mz with theoretical position
        expected_next_mz -= charge_dist;
      }
    }

    // ------------------------------------------------------------------------------
    // try to find a match for our isotopic peak on the right

    // redefine bounds
    lower_bound = precursor_spec.MZBegin(ms2_spec->getPrecursors()[0].getMZ());
    upper_bound = precursor_spec.MZEnd(fuzzy_upper_mz);

    expected_next_mz = precursor_peak.getMZ() + charge_dist;

    while (expected_next_mz < fuzzy_upper_mz)
    {
      // find nearest peak in precursor window
      const_spec_iterator np_it = precursor_spec.MZBegin(lower_bound, expected_next_mz, upper_bound);

      // handle border cases

      // check if next peak has smaller dist
      const_spec_iterator np_it2 = np_it;
      ++np_it;

      if (std::fabs(np_it2->getMZ() - expected_next_mz) < std::fabs(np_it->getMZ() - expected_next_mz))
      {
        np_it = np_it2;
      }

      // compute difference between found peak and expected
      double min_diff = std::fabs(np_it->getMZ() - expected_next_mz)  * 1000000 / expected_next_mz;

      // check if we found an isotopic peak
      if (min_diff < max_precursor_isotope_deviation_)
      {
#ifdef ISOBARIC_CHANNEL_EXTRACTOR_DEBUG
        std::cerr << "Mark peak as isotopic peak POS: " << precursor_spec[min_idx] << " (diff: " << min_diff << " vs " << max_precursor_isotope_deviation_ << ")" << std::endl;
#endif
        if (np_it->getMZ() < strict_upper_mz)
        {
          precursor_intensity += np_it->getIntensity();
        }
        else
        {
          // we're in the fuzzy area, so we will take only 50% of the given intensity
          // since we assume that the isolation window borders are not sharp
          precursor_intensity += 0.5 * np_it->getIntensity();
        }

        // update expected_next_mz
        expected_next_mz = np_it->getMZ() + charge_dist;
      }
      else
      {
#ifdef ISOBARIC_CHANNEL_EXTRACTOR_DEBUG
        std::cerr << "No matching isotopic peak for expected pos: " << expected_next_mz << " (min reached diff: " << min_diff << " vs " << max_precursor_isotope_deviation_ << ")" << std::endl;
#endif
        // update expected_next_mz with theoretical position
        expected_next_mz += charge_dist;
      }
    }

    // ------------------------------------------------------------------------------
    // compute total intensity
    int idx = static_cast<int>(precursor_peak_idx) - 1;
    while (idx >= 0 && precursor_spec[idx].getMZ() > fuzzy_lower_mz)
    {
      if (precursor_spec[idx].getMZ() > strict_lower_mz)
      {
        total_intensity += precursor_spec[idx].getIntensity();
      }
      else
      {
        // we're in the fuzzy area, so we will take only 50% of the given intensity
        // since we assume that the isolation window borders are not sharp
        total_intensity += 0.5 * precursor_spec[idx].getIntensity();
      }
      --idx;
    }

    idx = static_cast<int>(precursor_peak_idx) + 1;
    while (idx < static_cast<int>(precursor_spec.size()) && precursor_spec[idx].getMZ() < fuzzy_upper_mz)
    {
      if (precursor_spec[idx].getMZ() < strict_upper_mz)
      {
        total_intensity += precursor_spec[idx].getIntensity();
      }
      else
      {
        // we're in the fuzzy area, so we will take only 50% of the given intensity
        // since we assume that the isolation window borders are not sharp
        total_intensity += 0.5 * precursor_spec[idx].getIntensity();
      }
      ++idx;
    }

    return precursor_intensity / total_intensity;
  }

  double IsobaricChannelExtractor::computePrecursorPurity_(const PeakMap::ConstIterator& ms2_spec, const PuritySate_& pState) const
  {
    // we cannot analyze precursors without a charge
    if (ms2_spec->getPrecursors()[0].getCharge() == 0)
    {
      return 1.0;
    }
    else
    {
#ifdef ISOBARIC_CHANNEL_EXTRACTOR_DEBUG
      std::cerr << "------------------ analyzing " << ms2_spec->getNativeID() << std::endl;
#endif

      // compute purity of preceding ms1 scan
      double early_scan_purity = computeSingleScanPrecursorPurity_(ms2_spec, *(pState.precursorScan));

      if (pState.hasFollowUpScan && interpolate_precursor_purity_)
      {
        double late_scan_purity = computeSingleScanPrecursorPurity_(ms2_spec, *(pState.followUpScan));

        // calculating the extrapolated, S2I value as a time weighted linear combination of the two scans
        // see: Savitski MM, Sweetman G, Askenazi M, Marto JA, Lang M, Zinn N, et al. (2011).
        // Analytical chemistry 83: 8959–67. http://www.ncbi.nlm.nih.gov/pubmed/22017476
        // std::fabs is applied to compensate for potentially negative RTs
        return std::fabs(ms2_spec->getRT() - pState.precursorScan->getRT()) *
               ((late_scan_purity - early_scan_purity) / std::fabs(pState.followUpScan->getRT() - pState.precursorScan->getRT()))
               + early_scan_purity;
      }
      else
      {
        return early_scan_purity;
      }
    }
  }

  void IsobaricChannelExtractor::extractChannels(const PeakMap& ms_exp_data, ConsensusMap& consensus_map)
  {
    if (ms_exp_data.empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.\n";
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Experiment has no scans!");
    }

    // check if RT is sorted (we rely on it)
    if (!ms_exp_data.isSorted(false))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Spectra are not sorted in RT! Please sort them first!");
    }

    // clear the output map
    consensus_map.clear(false);
    consensus_map.setExperimentType("labeled_MS2");

    // create predicate for spectrum checking
    OPENMS_LOG_INFO << "Selecting scans with activation mode: " << selected_activation_ << std::endl;
    
    // Select the two possible HCD activation modes according to PSI-MS ontology: HCID and HCD
    if (selected_activation_ == "auto") 
    {
      selected_activation_ = Precursor::NamesOfActivationMethod[Precursor::HCID] + "," + Precursor::NamesOfActivationMethod[Precursor::HCD];
    }

    HasActivationMethod<PeakMap::SpectrumType> isValidActivation(ListUtils::create<String>(selected_activation_));

    // walk through spectra and count the number of scans with valid activation method per MS-level
    // only the highest level will be used for quantification (e.g. MS3, if present)
    std::map<UInt, UInt> ms_level;
    std::map<String, int> activation_modes;
    for (PeakMap::ConstIterator it = ms_exp_data.begin(); it != ms_exp_data.end(); ++it)
    {
      if (it->getMSLevel() == 1) continue; // never report MS1
      ++activation_modes[getActivationMethod_(*it)]; // count HCD, CID, ...
      if (selected_activation_ == "any" || isValidActivation(*it))
      {
        ++ms_level[it->getMSLevel()];
      }
    }
    if (ms_level.empty())
    {
      OPENMS_LOG_WARN << "Filtering by MS/MS(/MS) and activation mode: no spectra pass activation mode filter!\n"
               << "Activation modes found:\n";
      for (std::map<String, int>::const_iterator it = activation_modes.begin(); it != activation_modes.end(); ++it)
      {
        OPENMS_LOG_WARN << "  mode " << (it->first.empty() ? "<none>" : it->first) << ": " << it->second << " scans\n";
      }
      OPENMS_LOG_WARN << "Result will be empty!" << std::endl;
      return;
    }
    OPENMS_LOG_INFO << "Filtering by MS/MS(/MS) and activation mode:\n";
    for (std::map<UInt, UInt>::const_iterator it = ms_level.begin(); it != ms_level.end(); ++it)
    {
      OPENMS_LOG_INFO << "  level " << it->first << ": " << it->second << " scans\n";
    }
    UInt quant_ms_level = ms_level.rbegin()->first;
    OPENMS_LOG_INFO << "Using MS-level " << quant_ms_level << " for quantification." << std::endl;

    // now we have picked data
    // --> assign peaks to channels
    UInt64 element_index(0);

    // remember the current precursor spectrum
    PuritySate_ pState(ms_exp_data);

    typedef std::map<String, ChannelQC > ChannelQCSet;
    ChannelQCSet channel_mz_delta;
    const double qc_dist_mz = 0.5; // fixed! Do not change!

    Size number_of_channels = quant_method_->getNumberOfChannels();

    PeakMap::ConstIterator it_last_MS2 = ms_exp_data.end(); // remember last MS2 spec, to get precursor in MS1 (also if quant is in MS3)
    bool ms3 = false;
    for (PeakMap::ConstIterator it = ms_exp_data.begin(); it != ms_exp_data.end(); ++it)
    {
      // remember the last MS1 spectra as we assume it to be the precursor spectrum
      if (it->getMSLevel() ==  1)
      {
        // remember potential precursor and continue
        pState.precursorScan = it;
        // reset last MS2 -- we expect to see a new one soon and the old one should not be used for the following MS3 (if any)
        it_last_MS2 = ms_exp_data.end();
        continue;
      }

      if (it->getMSLevel() != quant_ms_level) continue;
      if ((*it).empty()) continue; // skip empty spectra
      if (!(selected_activation_ == "any" || isValidActivation(*it))) continue;

      // find following ms1 scan (needed for purity computation)
      if (!pState.followUpValid(it->getRT()))
      {
        // advance iterator
        pState.advanceFollowUp(it->getRT());
      }

      // check precursor constraints
      if (!isValidPrecursor_(it->getPrecursors()[0]))
      {
        OPENMS_LOG_DEBUG << "Skip spectrum " << it->getNativeID() << ": Precursor doesn't fulfill all constraints." << std::endl;
        continue;
      }

      // check precursor purity if we have a valid precursor ..
      double precursor_purity = -1.0;
      if (pState.precursorScan != ms_exp_data.end())
      {
        precursor_purity = computePrecursorPurity_(it, pState);
        // check if purity is high enough
        if (precursor_purity < min_precursor_purity_)
        {
          OPENMS_LOG_DEBUG << "Skip spectrum " << it->getNativeID() << ": Precursor purity is below the threshold. [purity = " << precursor_purity << "]" << std::endl;
          continue;
        }
      }
      else
      {
        OPENMS_LOG_INFO << "No precursor available for spectrum: " << it->getNativeID() << std::endl;
      }

      if (it->getMSLevel() == 3)
      {
        ms3 = true;
        // we cannot save just the last MS2 but need to compare to the precursor info stored in the (potential MS3 spectrum)
        it_last_MS2 = ms_exp_data.getPrecursorSpectrum(it);

        if (it_last_MS2 == ms_exp_data.end())
        { // this only happens if an MS3 spec does not have a preceding MS2
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("No MS2 precursor information given for MS3 scan native ID ") + it->getNativeID() + " with RT " + String(it->getRT()));
        }
      }
      else
      {
        it_last_MS2 = it;
      }

      // check if MS1 precursor info is available
      if (it_last_MS2->getPrecursors().empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("No precursor information given for scan native ID ") + it->getNativeID() + " with RT " + String(it->getRT()));
      }

      // store RT of MS2 scan and MZ of MS1 precursor ion as centroid of ConsensusFeature
      ConsensusFeature cf;
      cf.setUniqueId();
      cf.setRT(it_last_MS2->getRT());
      cf.setMZ(it_last_MS2->getPrecursors()[0].getMZ());

      Peak2D channel_value;
      channel_value.setRT(it->getRT());
      // for each each channel
      UInt64 map_index = 0;
      Peak2D::IntensityType overall_intensity = 0;

      for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator cl_it = quant_method_->getChannelInformation().begin();
            cl_it != quant_method_->getChannelInformation().end();
            ++cl_it)
      {
        // set mz-position of channel
        channel_value.setMZ(cl_it->center);
        // reset intensity
        channel_value.setIntensity(0);

        // as every evaluation requires time, we cache the MZEnd iterator
        const PeakMap::SpectrumType::ConstIterator mz_end = it->MZEnd(cl_it->center + qc_dist_mz);

        // search for the non-zero signal closest to theoretical position
        // & check for closest signal within reasonable distance (0.5 Da) -- might find neighbouring TMT channel, but that should not confuse anyone
        int peak_count(0); // count peaks in user window -- should be only one, otherwise Window is too large
        PeakMap::SpectrumType::ConstIterator idx_nearest(mz_end);
        for (PeakMap::SpectrumType::ConstIterator mz_it = it->MZBegin(cl_it->center - qc_dist_mz);
              mz_it != mz_end;
              ++mz_it)
        {
          if (mz_it->getIntensity() == 0) continue; // ignore 0-intensity shoulder peaks -- could be detrimental when de-calibrated
          double dist_mz = fabs(mz_it->getMZ() - cl_it->center);
          if (dist_mz < reporter_mass_shift_) ++peak_count;
          if (idx_nearest == mz_end // first peak
              || ((dist_mz < fabs(idx_nearest->getMZ() - cl_it->center)))) // closer to best candidate
          {
            idx_nearest = mz_it;
          }
        }
        if (idx_nearest != mz_end)
        {
          double mz_delta = cl_it->center - idx_nearest->getMZ();
          // stats: we don't care what shift the user specified
          channel_mz_delta[cl_it->name].mz_deltas.push_back(mz_delta);
          if (peak_count > 1) ++channel_mz_delta[cl_it->name].signal_not_unique;
          // pass user threshold
          if (std::fabs(mz_delta) < reporter_mass_shift_)
          {
            channel_value.setIntensity(idx_nearest->getIntensity());
          }
        }

        // discard contribution of this channel as it is below the required intensity threshold
        if (channel_value.getIntensity() < min_reporter_intensity_)
        {
          channel_value.setIntensity(0);
        }

        overall_intensity += channel_value.getIntensity();
        // add channel to ConsensusFeature
        cf.insert(map_index, channel_value, element_index);
        ++map_index;
      } // ! channel_iterator

      // check if we keep this feature or if it contains low-intensity quantifications
      if (remove_low_intensity_quantifications_ && hasLowIntensityReporter_(cf))
      {
        continue;
      }

      // check featureHandles are not empty
      if (overall_intensity <= 0)
      {
        cf.setMetaValue("all_empty", String("true"));
      }
      // add purity information if we could compute it
      if (precursor_purity > 0.0)
      {
        cf.setMetaValue("precursor_purity", precursor_purity);
      }

      // embed the id of the scan from which the quantitative information was extracted
      cf.setMetaValue("scan_id", it->getNativeID());
      // embed the id of the scan from which the ID information should be extracted
      // helpful for mapping later
      if (ms3)
      {
        cf.setMetaValue("id_scan_id", it_last_MS2->getNativeID());
      }
      // ...as well as additional meta information
      cf.setMetaValue("precursor_intensity", it->getPrecursors()[0].getIntensity());

      cf.setCharge(it_last_MS2->getPrecursors()[0].getCharge());
      cf.setIntensity(overall_intensity);
      consensus_map.push_back(cf);

      // the tandem-scan in the order they appear in the experiment
      ++element_index;
    } // ! Experiment iterator

    // print stats about m/z calibration / presence of signal
    OPENMS_LOG_INFO << "Calibration stats: Median distance of observed reporter ions m/z to expected position (up to " << qc_dist_mz << " Th):\n";
    bool impurities_found(false);
    for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator cl_it = quant_method_->getChannelInformation().begin();
      cl_it != quant_method_->getChannelInformation().end();
      ++cl_it)
    {
      OPENMS_LOG_INFO << "  ch " << String(cl_it->name).fillRight(' ', 4) << " (~" << String(cl_it->center).substr(0, 7).fillRight(' ', 7) << "): ";
      if (channel_mz_delta.find(cl_it->name) != channel_mz_delta.end())
      {
        // sort
        double median = Math::median(channel_mz_delta[cl_it->name].mz_deltas.begin(), channel_mz_delta[cl_it->name].mz_deltas.end(), false);
        if (((number_of_channels == 10) || (number_of_channels == 11)) &&
            (fabs(median) > TMT_10AND11PLEX_CHANNEL_TOLERANCE) &&
            (int(cl_it->center) != 126 && int(cl_it->center) != 131)) // these two channels have ~1 Th spacing.. so they do not suffer from the tolerance problem
        { // the channel was most likely empty, and we picked up the neighbouring channel's data (~0.006 Th apart). So reporting median here is misleading.
          OPENMS_LOG_INFO << "<invalid data (>" << TMT_10AND11PLEX_CHANNEL_TOLERANCE << " Th channel tolerance)>\n";
        }
        else
        {
          OPENMS_LOG_INFO << median << " Th";
          if (channel_mz_delta[cl_it->name].signal_not_unique > 0) 
          {
            OPENMS_LOG_INFO << " [MSn impurity (within " << reporter_mass_shift_ << " Th): " << channel_mz_delta[cl_it->name].signal_not_unique << " windows|spectra]";
            impurities_found = true;
          }
          OPENMS_LOG_INFO << "\n";
        }
      }
      else
      {
        OPENMS_LOG_INFO << "<no data>\n";
      }
    }
    if (impurities_found) OPENMS_LOG_INFO << "\nImpurities within the allowed reporter mass shift " << reporter_mass_shift_ << " Th have been found." 
                                   << "They can be ignored if the spectra are m/z calibrated (see above), since only the peak closest to the theoretical position is used for quantification!";
    OPENMS_LOG_INFO << std::endl;


    /// add meta information to the map
    registerChannelsInOutputMap_(consensus_map);
  }

  void IsobaricChannelExtractor::registerChannelsInOutputMap_(ConsensusMap& consensus_map)
  {
    // register the individual channels in the output consensus map
    Int index = 0;
    for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator cl_it = quant_method_->getChannelInformation().begin();
         cl_it != quant_method_->getChannelInformation().end();
         ++cl_it)
    {
      ConsensusMap::ColumnHeader channel_as_map;
      // label is the channel + description provided in the Params
      channel_as_map.label = quant_method_->getMethodName() + "_" + cl_it->name;

      // TODO(aiche): number of features need to be set later
      channel_as_map.size = consensus_map.size();

      // add some more MetaInfo
      channel_as_map.setMetaValue("channel_name", cl_it->name);
      channel_as_map.setMetaValue("channel_id", cl_it->id);
      channel_as_map.setMetaValue("channel_description", cl_it->description);
      channel_as_map.setMetaValue("channel_center", cl_it->center);
      consensus_map.getColumnHeaders()[index] = channel_as_map;
      ++index;
    }
  }

} // namespace
