// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Lars Nilse, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/MultiplexResolverAlgorithm.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>

using namespace std;

namespace OpenMS
{
  namespace
  {
    // names of the integer data arrays with the owners of blacklisted peaks (upper and lower 32 bits of the unique id)
    const std::string OWNER_HIGH = "owner_unique_id_high";
    const std::string OWNER_LOW = "owner_unique_id_low";

    // Was the blacklist spectrum acquired at the FAIMS CV of the multiplet? (Peaks at another CV cannot occlude it.)
    bool sameCompensationVoltage(const MSSpectrum& spectrum, const ConsensusFeature& consensus)
    {
      const bool feature_has_cv = consensus.metaValueExists(Constants::UserParam::FAIMS_CV);
      const bool spectrum_has_cv
        = spectrum.metaValueExists(Constants::UserParam::FAIMS_CV) || spectrum.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;
      if (feature_has_cv != spectrum_has_cv) { return false; }
      if (!feature_has_cv) { return true; }
      const double spectrum_cv = spectrum.metaValueExists(Constants::UserParam::FAIMS_CV)
                                   ? static_cast<double>(spectrum.getMetaValue(Constants::UserParam::FAIMS_CV))
                                   : spectrum.getDriftTime();
      return static_cast<double>(consensus.getMetaValue(Constants::UserParam::FAIMS_CV)) == spectrum_cv;
    }

    // a channel added for an identification (MS1LabeledWorkflow): not detected, intensity 0
    bool isAddedChannel(const FeatureHandle& handle)
    {
      return handle.getIntensity() == 0.0 && !handle.hasValidUniqueId();
    }

    // the label modifications of the samples, as they appear in a sequence string, e.g. "(Label:13C(6)15N(2))"
    std::vector<std::string> labelModifications(MultiplexDeltaMassesGenerator& generator)
    {
      std::vector<std::string> label_mods;
      for (const std::vector<std::string>& sample : generator.getSamplesLabelsList())
      {
        for (const std::string& label : sample)
        {
          if (label == "no_label") { continue; }
          const std::string name = generator.getLabelLong(label);
          if (name.empty()) { continue; }
          const std::string mod = "(" + name + ")";
          if (std::find(label_mods.begin(), label_mods.end(), mod) == label_mods.end()) { label_mods.push_back(mod); }
        }
      }
      return label_mods;
    }

    // the peptide of a hit at a charge: its sequence without the label modifications, and the charge
    std::string peptideOf(const std::string& sequence, int charge, const std::vector<std::string>& label_mods)
    {
      std::string peptide = sequence;
      for (const std::string& mod : label_mods) { StringUtils::substitute(peptide, mod, ""); }
      return peptide + "/" + std::to_string(charge);
    }
  }

  MultiplexResolverAlgorithm::MultiplexResolverAlgorithm() :
    DefaultParamHandler("MultiplexResolverAlgorithm")
  {
    defaults_.setValue("algorithm:labels", "[][Lys8,Arg10]", "Labels used for labelling the samples. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL");
    defaults_.setValue("algorithm:max_nr_labelled_aas", 0, "Maximum number of labelled amino acids per peptide, minus one. The algorithm searches for peptides with up to (this value + 1) labelled amino acids. For SILAC with trypsin digestion, this parameter corresponds to the maximum number of missed cleavages.");
    defaults_.setMinInt("algorithm:max_nr_labelled_aas", 0);
    defaults_.setValue("algorithm:mass_tolerance", 0.1, "Mass tolerance in Da for matching the mass shifts in the detected peptide multiplet to the theoretical mass shift pattern.", {"advanced"});
    // Continuous quantities, and read into double members: registering them as integers rejected
    // a fractional value on the command line and, worse, let an INI file carrying one through as 0.
    defaults_.setValue("algorithm:mz_tolerance", 10.0, "m/z tolerance in ppm for checking if dummy feature vicinity was blacklisted.", {"advanced"});
    defaults_.setMinFloat("algorithm:mz_tolerance", 0.0);
    defaults_.setValue("algorithm:rt_tolerance", 5.0, "Retention time tolerance in seconds for checking if dummy feature vicinity was blacklisted.", {"advanced"});
    defaults_.setMinFloat("algorithm:rt_tolerance", 0.0);
    defaults_.setValue("algorithm:use_blacklist_owners", "false",
                       "For a blacklist that records which feature reserved each peak (FeatureFinderMultiplex -out_blacklist).\n"
                       "true: a missing channel is not quantifiable (NaN) only if a peak reserved by another multiplet lies at its expected position, "
                       "and channels added for an identification (intensity 0, no feature id) are checked the same way. Peaks of the multiplet's own "
                       "channels (their isotopes) or of no reported feature (a pattern the detection dropped; FeatureFinderMultiplex records only the "
                       "first and the last pattern that reserved a peak) leave it at 0, which then means 'not occupied', not 'absent': use this only "
                       "if such channels are measured afterwards (as MS1LabeledWorkflow re-quantifies them).\n"
                       "false: the owners are ignored, every blacklisted peak counts.", {"advanced"});
    defaults_.setValidStrings("algorithm:use_blacklist_owners", {"true", "false"});
    defaults_.setSectionDescription("algorithm", "Parameters for the algorithm.");

    // one advanced parameter per known label: its mass shift
    MultiplexDeltaMassesGenerator generator;
    Param p = generator.getParameters();
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      const std::string label_name = "labels:" + it->name;
      defaults_.setValue(label_name, it->value, it->description, {"advanced"});
      defaults_.setMinFloat(label_name, 0.0);
    }
    defaults_.setSectionDescription("labels", "Isotopic labels that can be specified in section 'algorithm:labels'.");

    defaultsToParam_();
  }

  void MultiplexResolverAlgorithm::updateMembers_()
  {
    mass_tolerance_ = param_.getValue("algorithm:mass_tolerance");
    mz_tolerance_ = param_.getValue("algorithm:mz_tolerance");
    rt_tolerance_ = param_.getValue("algorithm:rt_tolerance");
    use_blacklist_owners_ = param_.getValue("algorithm:use_blacklist_owners").toBool();
  }

  double MultiplexResolverAlgorithm::deltaMassFromMapIndex_(const ConsensusFeature::HandleSetType& feature_handles, unsigned idx)
  {
    const double first_mass = feature_handles.begin()->getMZ() * feature_handles.begin()->getCharge();

    for (const FeatureHandle& handle : feature_handles)
    {
      if (handle.getMapIndex() == idx)
      {
        return handle.getMZ() * handle.getCharge() - first_mass;
      }
    }

    // no handle with this map index
    return numeric_limits<double>::quiet_NaN();
  }

  double MultiplexResolverAlgorithm::matchLabelSet_(const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                                    const MultiplexDeltaMasses::LabelSet& label_set, int& index_label_set)
  {
    for (auto it_mass_shift = pattern.begin(); it_mass_shift != pattern.end(); ++it_mass_shift)
    {
      if (it_mass_shift->label_set == label_set)
      {
        index_label_set = static_cast<int>(it_mass_shift - pattern.begin());
        return it_mass_shift->delta_mass;
      }
    }

    // no delta mass carries this label set
    return numeric_limits<double>::quiet_NaN();
  }

  bool MultiplexResolverAlgorithm::matchDeltaMasses_(const ConsensusFeature& consensus,
                                                     double detected_delta_mass_at_label_set,
                                                     const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                                     double theoretical_delta_mass_at_label_set,
                                                     std::vector<bool>& delta_mass_matched) const
  {
    const double first_mass = consensus.getFeatures().begin()->getMZ() * consensus.getFeatures().begin()->getCharge();

    // every detected feature needs a theoretical counterpart, relative to the identified feature
    for (const FeatureHandle& handle : consensus.getFeatures())
    {
      const double mass_shift_detected = (handle.getMZ() * handle.getCharge() - first_mass) - detected_delta_mass_at_label_set;
      bool matched = false;

      for (auto it_mass_shift = pattern.begin(); it_mass_shift != pattern.end(); ++it_mass_shift)
      {
        const double mass_shift_theoretical = it_mass_shift->delta_mass - theoretical_delta_mass_at_label_set;

        if (std::abs(mass_shift_detected - mass_shift_theoretical) < mass_tolerance_)
        {
          delta_mass_matched[it_mass_shift - pattern.begin()] = true;
          matched = true;
          break;
        }
      }

      if (!matched)
      {
        return false;
      }
    }

    return true;
  }

  bool MultiplexResolverAlgorithm::matchesResidues_(const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern, const AASequence& sequence,
                                                    const std::vector<std::vector<std::string>>& samples)
  {
    // the residue a metabolic label sits on ('\0' for labels that are not counted)
    auto residue_of = [](const std::string& label) -> char
    {
      if (label.rfind("Lys", 0) == 0) { return 'K'; }
      if (label.rfind("Arg", 0) == 0) { return 'R'; }
      if (label.rfind("Leu", 0) == 0) { return 'L'; }
      return '\0';
    };

    std::map<char, Size> residues;
    for (Size i = 0; i < sequence.size(); ++i)
    {
      ++residues[sequence[i].getOneLetterCode()[0]];
    }

    for (const MultiplexDeltaMasses::DeltaMass& delta_mass : pattern)
    {
      // labels per residue kind in this channel
      std::map<char, Size> labels;
      for (const std::string& label : delta_mass.label_set)
      {
        if (const char residue = residue_of(label); residue != '\0') { ++labels[residue]; }
      }
      // A kind the channel's sample (the one whose labels include the channel's) labels but the channel does not carry
      // must be absent from the peptide: the Lys8 channel of sample [Lys8,Arg10] fits a peptide with one K and no R.
      const auto sample = std::find_if(samples.begin(), samples.end(), [&delta_mass](const std::vector<std::string>& sample_labels)
      {
        return std::all_of(delta_mass.label_set.begin(), delta_mass.label_set.end(), [&sample_labels](const std::string& label)
        {
          return std::find(sample_labels.begin(), sample_labels.end(), label) != sample_labels.end();
        });
      });
      if (sample != samples.end())
      {
        for (const std::string& label : *sample)
        {
          if (const char residue = residue_of(label); residue != '\0') { labels.try_emplace(residue, 0); }
        }
      }
      for (const auto& [residue, count] : labels)
      {
        if (residues[residue] != count) { return false; }
      }
    }
    return true;
  }

  int MultiplexResolverAlgorithm::findMatchingPattern_(const ConsensusFeature& consensus,
                                                       double detected_delta_mass_at_label_set,
                                                       const AASequence& sequence,
                                                       const MultiplexDeltaMasses::LabelSet& label_set,
                                                       const std::vector<MultiplexDeltaMasses>& theoretical_patterns,
                                                       const std::vector<std::vector<std::string>>& samples,
                                                       std::vector<bool>& delta_mass_matched, int& index_label_set) const
  {
    for (auto it_pattern = theoretical_patterns.begin(); it_pattern != theoretical_patterns.end(); ++it_pattern)
    {
      const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern = it_pattern->getDeltaMasses();

      const double shift = matchLabelSet_(pattern, label_set, index_label_set);
      if (!std::isnan(shift) && matchesResidues_(pattern, sequence, samples))
      {
        // the label set occurs in this pattern; now all detected mass shifts have to fit as well
        delta_mass_matched.assign(delta_mass_matched.size(), false);

        if (matchDeltaMasses_(consensus, detected_delta_mass_at_label_set, pattern, shift, delta_mass_matched))
        {
          return static_cast<int>(it_pattern - theoretical_patterns.begin());
        }
      }
    }

    return -1;
  }

  double MultiplexResolverAlgorithm::findNewMZ_(double mz, int charge, const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                                const std::vector<bool>& delta_mass_matched)
  {
    // the first matched delta mass tells which pattern position the detected lightest feature has
    auto it_mass_shift = pattern.begin();
    auto it_delta_mass_matched = delta_mass_matched.begin();
    for (; it_mass_shift != pattern.end() && it_delta_mass_matched != delta_mass_matched.end();
         ++it_mass_shift, ++it_delta_mass_matched)
    {
      if (*it_delta_mass_matched)
      {
        return (mz * charge - it_mass_shift->delta_mass) / charge;
      }
    }

    // should never happen: at least the identified feature was matched
    return mz;
  }

  bool MultiplexResolverAlgorithm::isBlacklisted_(const MSExperiment& blacklist,
                                                  double rt,
                                                  double mz,
                                                  size_t charge,
                                                  const ConsensusFeature& consensus) const
  {
    const double mz_tolerance = mz_tolerance_ * mz / 1000000; // m/z tolerance in Da

    const MSExperiment::ConstIterator it_rt_begin = blacklist.RTBegin(rt - rt_tolerance_);
    const MSExperiment::ConstIterator it_rt_end = blacklist.RTEnd(rt + rt_tolerance_);

    // loop over range of relevant spectra
    for (MSExperiment::ConstIterator it_rt = it_rt_begin; it_rt < it_rt_end; ++it_rt)
    {
      if (!sameCompensationVoltage(*it_rt, consensus)) { continue; }
      // check the first three isotopes of the dummy feature
      for (size_t isotope = 0; isotope < 3; ++isotope)
      {
        const double mz_isotope = mz + isotope * Constants::C13C12_MASSDIFF_U / charge;

        const MSSpectrum::ConstIterator it_mz = it_rt->MZBegin(mz_isotope);
        if (it_mz == it_rt->end())
        {
          continue;
        }

        if ((std::abs(it_mz->getMZ() - mz_isotope)) < mz_tolerance)
        {
          // a blacklisted peak close-by
          return true;
        }
      }
    }

    // none of the first three isotopes has a blacklisted peak near-by
    return false;
  }

  bool MultiplexResolverAlgorithm::isReservedByOther_(const MSExperiment& blacklist, double rt, double mz, size_t charge,
                                                      const ConsensusFeature& consensus, const OwnerPeptides& peptides,
                                                      OwnerCheckCounts_& counts) const
  {
    const double mz_tolerance = mz_tolerance_ * mz / 1000000; // m/z tolerance in Da
    bool other_peptide = false, unidentified = false, own = false, same_peptide = false, orphan = false;
    const MSExperiment::ConstIterator it_rt_end = blacklist.RTEnd(rt + rt_tolerance_);
    for (auto it_rt = blacklist.RTBegin(rt - rt_tolerance_); it_rt < it_rt_end; ++it_rt)
    {
      if (!sameCompensationVoltage(*it_rt, consensus)) { continue; }
      for (size_t isotope = 0; isotope < 3; ++isotope)
      {
        const double mz_isotope = mz + isotope * Constants::C13C12_MASSDIFF_U / charge;
        for (auto it_mz = it_rt->MZBegin(mz_isotope - mz_tolerance); it_mz != it_rt->end() && it_mz->getMZ() < mz_isotope + mz_tolerance; ++it_mz)
        {
          if (std::abs(it_mz->getMZ() - mz_isotope) >= mz_tolerance) { continue; }
          switch (classifyBlacklistOwner(blacklistOwner(*it_rt, it_mz - it_rt->begin()), consensus, peptides))
          {
            case OwnerRelation::OTHER_PEPTIDE: other_peptide = true; break;
            case OwnerRelation::SAME_PEPTIDE: same_peptide = true; break;
            case OwnerRelation::UNIDENTIFIED: unidentified = true; break;
            case OwnerRelation::OWN: own = true; break;
            case OwnerRelation::ORPHAN: orphan = true; break;
          }
        }
      }
    }
    // the same rule as isForeign()
    if (other_peptide) { ++counts.other_peptide; return true; }
    if (same_peptide) { ++counts.same_peptide; return true; }
    if (unidentified) { ++counts.unidentified; return true; }
    if (own || orphan)
    {
      ++counts.usable;
      counts.own += own;
      counts.orphan += orphan;
    }
    return false;
  }

  bool MultiplexResolverAlgorithm::hasBlacklistOwners(const MSExperiment& blacklist)
  {
    bool peaks = false;
    for (const MSSpectrum& spectrum : blacklist)
    {
      if (spectrum.empty()) { continue; }
      const MSSpectrum::IntegerDataArrays& arrays = spectrum.getIntegerDataArrays();
      const auto high = getDataArrayByName(arrays, OWNER_HIGH);
      const auto low = getDataArrayByName(arrays, OWNER_LOW);
      if (high == arrays.end() || low == arrays.end() || high->size() != spectrum.size() || low->size() != spectrum.size()) { return false; }
      peaks = true;
    }
    return peaks;
  }

  UInt64 MultiplexResolverAlgorithm::blacklistOwner(const MSSpectrum& spectrum, Size peak_index)
  {
    const MSSpectrum::IntegerDataArrays& arrays = spectrum.getIntegerDataArrays();
    const auto high = getDataArrayByName(arrays, OWNER_HIGH);
    const auto low = getDataArrayByName(arrays, OWNER_LOW);
    if (high == arrays.end() || low == arrays.end() || peak_index >= high->size() || peak_index >= low->size())
    {
      return UniqueIdInterface::INVALID;
    }
    return (static_cast<UInt64>(static_cast<UInt32>((*high)[peak_index])) << 32) | static_cast<UInt32>((*low)[peak_index]);
  }

  void MultiplexResolverAlgorithm::setBlacklistOwners(MSSpectrum& spectrum, const std::vector<UInt64>& owners)
  {
    MSSpectrum::IntegerDataArrays& arrays = spectrum.getIntegerDataArrays();
    for (const std::string* name : {&OWNER_HIGH, &OWNER_LOW})
    {
      if (getDataArrayByName(arrays, *name) == arrays.end())
      {
        arrays.emplace_back();
        arrays.back().setName(*name);
      }
    }
    auto high = getDataArrayByName(arrays, OWNER_HIGH);
    auto low = getDataArrayByName(arrays, OWNER_LOW);
    high->resize(owners.size());
    low->resize(owners.size());
    for (Size i = 0; i < owners.size(); ++i)
    {
      (*high)[i] = static_cast<Int>(static_cast<UInt32>(owners[i] >> 32));
      (*low)[i] = static_cast<Int>(static_cast<UInt32>(owners[i]));
    }
  }

  MultiplexResolverAlgorithm::OwnerPeptides MultiplexResolverAlgorithm::ownerPeptides(const ConsensusMap& multiplets) const
  {
    MultiplexDeltaMassesGenerator generator = createGenerator_();
    const std::vector<std::string> label_mods = labelModifications(generator);

    OwnerPeptides peptides;
    for (const ConsensusFeature& cf : multiplets)
    {
      const auto linked = cf.getLinkedIdentifications(multiplets.getIdentificationData());
      if (linked.empty() || linked[0].matches.empty()) { continue; }
      const std::string peptide = peptideOf(linked[0].matches[0]->representation, cf.getCharge(), label_mods);
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        if (handle.hasValidUniqueId()) { peptides.emplace(handle.getUniqueId(), peptide); }
      }
    }
    return peptides;
  }

  std::string MultiplexResolverAlgorithm::peptideKey(const std::string& sequence, int charge) const
  {
    MultiplexDeltaMassesGenerator generator = createGenerator_();
    return peptideOf(sequence, charge, labelModifications(generator));
  }

  MultiplexResolverAlgorithm::OwnerRelation MultiplexResolverAlgorithm::classifyBlacklistOwner(UInt64 owner, const std::string& peptide,
                                                                                            const OwnerPeptides& peptides)
  {
    if (owner == UniqueIdInterface::INVALID) { return OwnerRelation::ORPHAN; }
    const auto it_owner = peptides.find(owner);
    if (it_owner == peptides.end()) { return OwnerRelation::UNIDENTIFIED; }
    return it_owner->second == peptide ? OwnerRelation::SAME_PEPTIDE : OwnerRelation::OTHER_PEPTIDE;
  }

  MultiplexResolverAlgorithm::OwnerRelation MultiplexResolverAlgorithm::classifyBlacklistOwner(UInt64 owner, const ConsensusFeature& consensus,
                                                                                            const OwnerPeptides& peptides)
  {
    if (owner == UniqueIdInterface::INVALID) { return OwnerRelation::ORPHAN; }
    for (const FeatureHandle& handle : consensus.getFeatures())
    {
      if (handle.hasValidUniqueId() && handle.getUniqueId() == owner) { return OwnerRelation::OWN; }
    }
    const auto it_owner = peptides.find(owner);
    if (it_owner == peptides.end()) { return OwnerRelation::UNIDENTIFIED; }
    // the peptide of the multiplet is that of its detected channels
    for (const FeatureHandle& handle : consensus.getFeatures())
    {
      if (!handle.hasValidUniqueId()) { continue; }
      const auto it_own = peptides.find(handle.getUniqueId());
      if (it_own != peptides.end())
      {
        return it_own->second == it_owner->second ? OwnerRelation::SAME_PEPTIDE : OwnerRelation::OTHER_PEPTIDE;
      }
    }
    return OwnerRelation::OTHER_PEPTIDE;
  }

  bool MultiplexResolverAlgorithm::isForeign(OwnerRelation relation)
  {
    return relation == OwnerRelation::UNIDENTIFIED || relation == OwnerRelation::SAME_PEPTIDE || relation == OwnerRelation::OTHER_PEPTIDE;
  }

  ConsensusFeature MultiplexResolverAlgorithm::completeConsensus_(const ConsensusFeature& consensus,
                                                                  const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                                                  const std::vector<bool>& delta_mass_matched, int index_label_set,
                                                                  const MSExperiment& blacklist, const OwnerPeptides* peptides,
                                                                  OwnerCheckCounts_* counts) const
  {
    OwnerCheckCounts_ no_counts;
    if (counts == nullptr) { counts = &no_counts; }
    // with owners, a channel added for an identification is not quantifiable if another multiplet occupies its position
    auto checked = [&](FeatureHandle handle)
    {
      if (peptides != nullptr && isAddedChannel(handle))
      {
        ++counts->added;
        if (isReservedByOther_(blacklist, consensus.getRT(), handle.getMZ(), consensus.getCharge(), consensus, *peptides, *counts))
        {
          ++counts->added_vetoed;
          handle.setIntensity(std::numeric_limits<double>::quiet_NaN());
        }
      }
      return handle;
    };

    // nothing to complete: the detected consensus is already complete
    if (consensus.size() == pattern.size())
    {
      ConsensusFeature complete(consensus);
      if (peptides != nullptr)
      {
        ConsensusFeature::HandleSetType handles;
        for (const FeatureHandle& handle : consensus.getFeatures()) { handles.insert(checked(handle)); }
        complete.setFeatures(std::move(handles));
      }
      return complete;
    }

    if (pattern.size() != delta_mass_matched.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, delta_mass_matched.size(), "pattern size does not match delta_mass_matched size");
    }

    ConsensusFeature consensus_complete;
    // Completion creates a new feature identity, but its acquisition CV remains the same.
    if (consensus.metaValueExists(Constants::UserParam::FAIMS_CV))
    {
      consensus_complete.setMetaValue(Constants::UserParam::FAIMS_CV, consensus.getMetaValue(Constants::UserParam::FAIMS_CV));
    }

    const int charge = consensus.getCharge();
    const double rt = consensus.getRT();
    const double mz = consensus.getMZ();

    // m/z of the lightest peptide of the complete multiplet
    const double mz_complete = findNewMZ_(mz, charge, pattern, delta_mass_matched);

    consensus_complete.setMZ(mz_complete);
    consensus_complete.setRT(consensus.getRT());
    consensus_complete.setCharge(consensus.getCharge());
    consensus_complete.setIntensity(consensus.getIntensity()); // alternatively, reduce intensity due to new zero-intensity dummy features
    consensus_complete.setQuality(consensus.getQuality());
    // the identifications of the multiplet (resolve() records the new map index in the first match)
    consensus_complete.getIDQueries() = consensus.getIDQueries();
    consensus_complete.getIDMatches() = consensus.getIDMatches();

    // walk the theoretical pattern: copy detected features, construct dummy features for the missing ones
    auto it_mass_shift = pattern.begin();
    auto it_delta_mass_matched = delta_mass_matched.begin();
    auto it_feature = consensus.getFeatures().begin();
    for (; it_mass_shift != pattern.end() && it_delta_mass_matched != delta_mass_matched.end();
         ++it_mass_shift, ++it_delta_mass_matched)
    {
      if (*it_delta_mass_matched)
      {
        // copy feature from incomplete consensus
        FeatureHandle feature_handle = checked(*it_feature);
        feature_handle.setMapIndex(it_mass_shift - pattern.begin());
        consensus_complete.insert(feature_handle);

        if (it_feature != consensus.getFeatures().end())
        {
          ++it_feature;
        }
      }
      else
      {
        // construct dummy feature
        FeatureHandle feature_handle;
        const double mz_dummy = mz_complete + it_mass_shift->delta_mass / charge;
        feature_handle.setMZ(mz_dummy);
        feature_handle.setRT(rt);
        bool blacklisted = false;
        if (peptides != nullptr)
        {
          ++counts->dummies;
          counts->anonymous += isBlacklisted_(blacklist, rt, mz_dummy, charge, consensus);
          blacklisted = isReservedByOther_(blacklist, rt, mz_dummy, charge, consensus, *peptides, *counts);
          counts->dummies_vetoed += blacklisted;
        }
        else
        {
          blacklisted = isBlacklisted_(blacklist, rt, mz_dummy, charge, consensus);
        }
        if (blacklisted)
        {
          // Some peaks close-by were blacklisted during feature detection (by another multiplet, if owners are used), i.e. another peptide feature overlaps with the dummy feature.
          // Consequently, we better report NaN i.e. not quantifiable.
          feature_handle.setIntensity(std::numeric_limits<double>::quiet_NaN());
        }
        else
        {
          // There is no blacklisted peak near-by, i.e. there is no peptide feature in the vicinity: the peptide is absent.
          // (With owners: no other multiplet's, so the caller measures the channel; see algorithm:use_blacklist_owners.)
          feature_handle.setIntensity(0.0);
        }
        feature_handle.setCharge(charge);
        feature_handle.setMapIndex(it_mass_shift - pattern.begin());
        consensus_complete.insert(feature_handle);
      }
    }

    return consensus_complete;
  }

  MultiplexDeltaMassesGenerator MultiplexResolverAlgorithm::createGenerator_() const
  {
    // mass shift of every label, as configured in section 'labels'
    std::map<std::string, double> label_mass_shift;
    Param label_params = param_.copy("labels:", true);
    for (Param::ParamIterator it = label_params.begin(); it != label_params.end(); ++it)
    {
      label_mass_shift.insert(make_pair(it->name, static_cast<double>(it->value)));
    }

    return MultiplexDeltaMassesGenerator(param_.getValue("algorithm:labels").toString(),
                                         param_.getValue("algorithm:max_nr_labelled_aas"), label_mass_shift);
  }

  std::vector<bool> MultiplexResolverAlgorithm::checkIdentifications(const ConsensusFeature& consensus, const IdentificationData& data) const
  {
    const std::vector<IdentificationData::QueryMatches> identifications = consensus.getLinkedIdentifications(data);
    std::vector<bool> consistent(identifications.size(), false);
    if (consensus.getFeatures().empty()) { return consistent; }

    MultiplexDeltaMassesGenerator generator = createGenerator_();
    const std::vector<MultiplexDeltaMasses> theoretical_masses = generator.getDeltaMassesList();
    const std::vector<std::vector<std::string>> samples = generator.getSamplesLabelsList();
    const size_t multiplicity = theoretical_masses[0].getDeltaMasses().size();

    for (Size i = 0; i < identifications.size(); ++i)
    {
      const IdentificationData::QueryMatches& identification = identifications[i];
      if (identification.matches.empty() || !identification.query->metaValueExists("map_index")) { continue; }
      const double detected_delta_mass = deltaMassFromMapIndex_(consensus.getFeatures(), identification.query->getMetaValue("map_index"));
      if (std::isnan(detected_delta_mass)) { continue; }

      const AASequence sequence = AASequence::fromString(identification.matches[0]->representation);
      const MultiplexDeltaMasses::LabelSet label_set = generator.extractLabelSet(sequence);
      std::vector<bool> delta_mass_matched(multiplicity, false);
      int index_label_set = -1;
      consistent[i] = findMatchingPattern_(consensus, detected_delta_mass, sequence, label_set, theoretical_masses, samples,
                                           delta_mass_matched, index_label_set) >= 0;
    }
    return consistent;
  }

  std::vector<double> MultiplexResolverAlgorithm::expectedChannelMZs(const AASequence& sequence, int charge, Size& position) const
  {
    std::vector<MultiplexDeltaMasses::LabelSet> label_sets;
    return expectedChannelMZs(sequence, charge, position, label_sets);
  }

  std::vector<double> MultiplexResolverAlgorithm::expectedChannelMZs(const AASequence& sequence, int charge, Size& position,
                                                                     std::vector<MultiplexDeltaMasses::LabelSet>& label_sets) const
  {
    label_sets.clear();
    if (charge == 0) { return {}; }
    MultiplexDeltaMassesGenerator generator = createGenerator_();
    const MultiplexDeltaMasses::LabelSet label_set = generator.extractLabelSet(sequence);
    const std::vector<std::vector<std::string>> samples = generator.getSamplesLabelsList();
    for (const MultiplexDeltaMasses& theoretical : generator.getDeltaMassesList())
    {
      const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern = theoretical.getDeltaMasses();
      int index_label_set = -1;
      const double shift = matchLabelSet_(pattern, label_set, index_label_set);
      if (std::isnan(shift) || !matchesResidues_(pattern, sequence, samples)) { continue; }

      const double mz = sequence.getMZ(charge);
      std::vector<double> mzs;
      for (const MultiplexDeltaMasses::DeltaMass& delta_mass : pattern)
      {
        mzs.push_back(mz + (delta_mass.delta_mass - shift) / std::abs(charge));
        label_sets.push_back(delta_mass.label_set);
      }
      position = static_cast<Size>(index_label_set);
      return mzs;
    }
    return {};
  }

  namespace
  {
    using ID = IdentificationData;

    /**
      @brief Remove from the identification data of @p map what only consensus features of @p original that @p map
      does not have link

      Those identifications were assigned in @p original, so they are no unassigned ones of @p map. Matches that
      were unassigned in @p original stay.
    */
    void removeOtherFeatureIdentifications(ConsensusMap& map, const ConsensusMap& original)
    {
      ID& data = map.getIdentificationData();
      std::set<ID::MatchReference> kept_matches;
      std::set<ID::QueryReference> kept_queries;
      for (const ConsensusFeature& feature : map)
      {
        kept_matches.insert(feature.getIDMatches().begin(), feature.getIDMatches().end());
        const auto queries = feature.getLinkedIDQueries(data);
        kept_queries.insert(queries.begin(), queries.end());
      }
      std::set<ID::MatchReference> removed_matches;
      std::set<ID::QueryReference> removed_queries;
      for (const ConsensusFeature& feature : original)
      {
        for (const auto& match : feature.getIDMatches())
        {
          if (!kept_matches.contains(match)) removed_matches.insert(match);
        }
        for (const auto& query : feature.getLinkedIDQueries(original.getIdentificationData()))
        {
          if (!kept_queries.contains(query)) removed_queries.insert(query);
        }
      }
      if (removed_matches.empty() && removed_queries.empty()) return;
      data.eraseMatches([&](const ID::Run& run, const ID::Identification&, const ID::Match& match) {
        return removed_matches.contains({run.getUuid(), match.getId()});
      });
      // identifications that have no unassigned matches left
      data.eraseIdentifications([&](const ID::Run& run, const ID::Identification& query) {
        return query.getMatches().empty() && removed_queries.contains({run.getUuid(), query.getId()});
      });
    }
  } // namespace

  void MultiplexResolverAlgorithm::resolve(const ConsensusMap& map_in, ConsensusMap& map_out, ConsensusMap& map_conflicts,
                                           const MSExperiment& blacklist, bool keep_conflicting_identifications) const
  {
    MultiplexDeltaMassesGenerator generator = createGenerator_();

    const ConsensusMap& in = map_in;
    const ID& data = in.getIdentificationData();

    // both outputs start as empty copies of the input, i.e. with its meta data
    map_out = in;
    map_conflicts = in;
    map_out.resize(0);
    map_conflicts.resize(0);

    const std::vector<MultiplexDeltaMasses> theoretical_masses = generator.getDeltaMassesList();
    const std::vector<std::vector<std::string>> samples = generator.getSamplesLabelsList();
    const size_t multiplicity = theoretical_masses[0].getDeltaMasses().size();

    // with the owners of the blacklisted peaks, only peaks of other multiplets make a channel not quantifiable
    const bool owners = use_blacklist_owners_ && hasBlacklistOwners(blacklist);
    OwnerPeptides peptides;
    if (owners) { peptides = ownerPeptides(in); }
    OwnerCheckCounts_ counts;

    for (const ConsensusFeature& cf : in)
    {
      // consensus features without sequence annotations are written unchanged to the conflict output
      const auto linked = cf.getLinkedIdentifications(data);
      if (linked.empty() || linked[0].matches.empty())
      {
        map_conflicts.push_back(cf);
        continue;
      }

      // extract the label set from the attached peptide sequence (only the first one is considered)
      const ID::QueryMatches& identification = linked[0];
      const ID::Match& match = *identification.matches[0];
      const AASequence sequence = AASequence::fromString(match.representation);
      const MultiplexDeltaMasses::LabelSet label_set = generator.extractLabelSet(sequence);
      if (!identification.query->metaValueExists("map_index"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The meta value 'map_index' is missing in the input data. In the IDMapper tool, please set the advanced parameter consensus:annotate_ids_with_subelements = true.");
      }
      const double detected_delta_mass = deltaMassFromMapIndex_(cf.getFeatures(), identification.query->getMetaValue("map_index"));
      if (std::isnan(detected_delta_mass))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No delta mass with this map_index could be found.", "");
      }
      std::vector<bool> delta_mass_matched(multiplicity, false);
      int index_label_set = -1;

      const int index = findMatchingPattern_(cf, detected_delta_mass, sequence, label_set, theoretical_masses, samples,
                                             delta_mass_matched, index_label_set);

      if (index >= 0)
      {
        const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern = theoretical_masses[index].getDeltaMasses();
        map_out.push_back(completeConsensus_(cf, pattern, delta_mass_matched, index_label_set, blacklist, owners ? &peptides : nullptr, &counts));
        if (cf.size() != pattern.size())
        {
          // the new map index of the identified feature
          ID::Run& run = map_out.getIdentificationData().getRun(identification.run->getIdentifier());
          ID::MatchData annotated = match.getData();
          annotated.setMetaValue("map_index", index_label_set);
          run.replaceMatch(match.getId(), annotated);
        }
      }
      else
      {
        map_conflicts.push_back(cf);
      }
    }

    // the identifications of the consensus features of the other output are not unassigned ones (unless those of the
    // conflicting ones are to be kept)
    if (!keep_conflicting_identifications) removeOtherFeatureIdentifications(map_out, in);
    removeOtherFeatureIdentifications(map_conflicts, in);

    if (owners && counts.dummies + counts.added > 0)
    {
      OPENMS_LOG_INFO << "Blacklist with owners: " << counts.dummies_vetoed << " of " << counts.dummies << " missing channel(s) and "
                      << counts.added_vetoed << " of " << counts.added << " channel(s) added for identifications are not quantifiable ("
                      << counts.other_peptide << " occupied by another identified peptide, " << counts.same_peptide
                      << " by another multiplet of the same peptide, " << counts.unidentified << " only by unidentified multiplets); "
                      << counts.usable << " lie on peaks reserved only by their own multiplet's channels (" << counts.own
                      << ") or by no reported feature (" << counts.orphan << "). Without owners, " << counts.anonymous
                      << " of the missing channels and none of the added ones would be not quantifiable." << endl;
    }

    // update map sizes
    for (unsigned map_index = 0; map_index < multiplicity; ++map_index)
    {
      map_out.getColumnHeaders()[map_index].size = map_out.size();
    }

    map_out.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    map_conflicts.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  }
} // namespace OpenMS
