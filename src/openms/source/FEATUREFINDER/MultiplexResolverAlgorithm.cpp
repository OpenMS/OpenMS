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
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <cmath>
#include <limits>

using namespace std;

namespace OpenMS
{
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
                                                     const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                                     double theoretical_delta_mass_at_label_set,
                                                     std::vector<bool>& delta_mass_matched) const
  {
    const double first_mass = consensus.getFeatures().begin()->getMZ() * consensus.getFeatures().begin()->getCharge();
    const PeptideIdentification& id = consensus.getPeptideIdentifications()[0];
    if (!id.metaValueExists("map_index"))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The meta value 'map_index' is missing in the input data. In the IDMapper tool, please set the advanced parameter consensus:annotate_ids_with_subelements = true.");
    }
    const double detected_delta_mass_at_label_set = deltaMassFromMapIndex_(consensus.getFeatures(), id.getMetaValue("map_index"));
    if (std::isnan(detected_delta_mass_at_label_set))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No delta mass with this map_index could be found.", "");
    }

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

  int MultiplexResolverAlgorithm::findMatchingPattern_(const ConsensusFeature& consensus,
                                                       const MultiplexDeltaMasses::LabelSet& label_set,
                                                       const std::vector<MultiplexDeltaMasses>& theoretical_patterns,
                                                       std::vector<bool>& delta_mass_matched, int& index_label_set) const
  {
    for (auto it_pattern = theoretical_patterns.begin(); it_pattern != theoretical_patterns.end(); ++it_pattern)
    {
      const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern = it_pattern->getDeltaMasses();

      const double shift = matchLabelSet_(pattern, label_set, index_label_set);
      if (!std::isnan(shift))
      {
        // the label set occurs in this pattern; now all detected mass shifts have to fit as well
        delta_mass_matched.assign(delta_mass_matched.size(), false);

        if (matchDeltaMasses_(consensus, pattern, shift, delta_mass_matched))
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

  bool MultiplexResolverAlgorithm::isBlacklisted_(const MSExperiment& blacklist, double rt, double mz, size_t charge) const
  {
    const double mz_tolerance = mz_tolerance_ * mz / 1000000; // m/z tolerance in Da

    const MSExperiment::ConstIterator it_rt_begin = blacklist.RTBegin(rt - rt_tolerance_);
    const MSExperiment::ConstIterator it_rt_end = blacklist.RTEnd(rt + rt_tolerance_);

    // loop over range of relevant spectra
    for (MSExperiment::ConstIterator it_rt = it_rt_begin; it_rt < it_rt_end; ++it_rt)
    {
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

  ConsensusFeature MultiplexResolverAlgorithm::completeConsensus_(const ConsensusFeature& consensus,
                                                                  const std::vector<MultiplexDeltaMasses::DeltaMass>& pattern,
                                                                  const std::vector<bool>& delta_mass_matched, int index_label_set,
                                                                  const MSExperiment& blacklist) const
  {
    // nothing to do: the detected consensus is already complete
    if (consensus.size() == pattern.size())
    {
      return ConsensusFeature(consensus);
    }

    if (pattern.size() != delta_mass_matched.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, delta_mass_matched.size(), "pattern size does not match delta_mass_matched size");
    }

    ConsensusFeature consensus_complete;

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
    consensus_complete.setPeptideIdentifications(consensus.getPeptideIdentifications());
    consensus_complete.getPeptideIdentifications()[0].getHits()[0].setMetaValue("map_index", index_label_set);

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
        FeatureHandle feature_handle(*it_feature);
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
        if (isBlacklisted_(blacklist, rt, mz_dummy, charge))
        {
          // Some peaks close-by were blacklisted during feature detection, i.e. another peptide feature overlaps with the dummy feature.
          // Consequently, we better report NaN i.e. not quantifiable.
          feature_handle.setIntensity(std::numeric_limits<double>::quiet_NaN());
        }
        else
        {
          // There is no blacklisted peak near-by, i.e. there is no peptide feature in the vicinity.
          // Consequently, we can confidently report zero, i.e. the peptide is absent.
          feature_handle.setIntensity(0.0);
        }
        feature_handle.setCharge(charge);
        feature_handle.setMapIndex(it_mass_shift - pattern.begin());
        consensus_complete.insert(feature_handle);
      }
    }

    return consensus_complete;
  }

  void MultiplexResolverAlgorithm::resolve(const ConsensusMap& map_in, ConsensusMap& map_out, ConsensusMap& map_conflicts,
                                           const MSExperiment& blacklist) const
  {
    // mass shift of every label, as configured in section 'labels'
    std::map<std::string, double> label_mass_shift;
    Param label_params = param_.copy("labels:", true);
    for (Param::ParamIterator it = label_params.begin(); it != label_params.end(); ++it)
    {
      label_mass_shift.insert(make_pair(it->name, static_cast<double>(it->value)));
    }

    MultiplexDeltaMassesGenerator generator(param_.getValue("algorithm:labels").toString(),
                                            param_.getValue("algorithm:max_nr_labelled_aas"), label_mass_shift);

    // both outputs start as empty copies of the input, i.e. with its meta data
    map_out = map_in;
    map_conflicts = map_in;
    map_out.resize(0);
    map_conflicts.resize(0);

    const std::vector<MultiplexDeltaMasses> theoretical_masses = generator.getDeltaMassesList();
    const size_t multiplicity = theoretical_masses[0].getDeltaMasses().size();

    for (const ConsensusFeature& cf : map_in)
    {
      // consensus features without sequence annotations are written unchanged to the conflict output
      if (cf.getPeptideIdentifications().empty() || cf.getPeptideIdentifications()[0].getHits().empty())
      {
        map_conflicts.push_back(cf);
        continue;
      }

      // extract the label set from the attached peptide sequence (only the first one is considered)
      const AASequence& sequence = cf.getPeptideIdentifications()[0].getHits()[0].getSequence();
      const MultiplexDeltaMasses::LabelSet label_set = generator.extractLabelSet(sequence);
      std::vector<bool> delta_mass_matched(multiplicity, false);
      int index_label_set = -1;

      const int index = findMatchingPattern_(cf, label_set, theoretical_masses, delta_mass_matched, index_label_set);

      if (index >= 0)
      {
        map_out.push_back(completeConsensus_(cf, theoretical_masses[index].getDeltaMasses(), delta_mass_matched, index_label_set, blacklist));
      }
      else
      {
        map_conflicts.push_back(cf);
      }
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
