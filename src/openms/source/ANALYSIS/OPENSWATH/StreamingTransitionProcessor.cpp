// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/StreamingTransitionProcessor.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <chrono>

namespace OpenMS
{

  StreamingTransitionProcessor::StreamingTransitionProcessor() :
    DefaultParamHandler("StreamingTransitionProcessor"),
    ProgressLogger(),
    rng_(std::chrono::system_clock::now().time_since_epoch().count())
  {
    defaults_.setValue("product_mz_threshold", 0.05, "MZ threshold for product ion comparison");
    defaults_.setValue("precursor_mz_threshold", 0.05, "MZ threshold for precursor ion comparison");
    defaults_.setValue("min_product_mz", 350.0, "Minimum product m/z");
    defaults_.setValue("max_product_mz", 2000.0, "Maximum product m/z");
    defaults_.setValue("max_transitions", 6, "Maximum number of transitions per peptide");
    defaults_.setValue("generate_decoys", "true", "Generate decoy transitions");
    defaults_.setValidStrings("generate_decoys", {"true", "false"});
    defaults_.setValue("decoy_tag", "DECOY_", "Prefix tag for decoy identifiers");
    defaults_.setValue("precursor_mz_shift", 0.0, "m/z shift for decoy precursor");
    defaults_.setValue("product_mz_shift", 20.0, "m/z shift for decoy product ions");

    defaultsToParam_();
  }

  StreamingTransitionProcessor::~StreamingTransitionProcessor() = default;

  void StreamingTransitionProcessor::updateMembers_()
  {
    product_mz_threshold_ = param_.getValue("product_mz_threshold");
    precursor_mz_threshold_ = param_.getValue("precursor_mz_threshold");
    min_product_mz_ = param_.getValue("min_product_mz");
    max_product_mz_ = param_.getValue("max_product_mz");
    max_transitions_ = param_.getValue("max_transitions");
    generate_decoys_ = param_.getValue("generate_decoys").toString() == "true";
    decoy_tag_ = param_.getValue("decoy_tag").toString();
    precursor_mz_shift_ = param_.getValue("precursor_mz_shift");
    product_mz_shift_ = param_.getValue("product_mz_shift");
  }

  void StreamingTransitionProcessor::process(const String& input_file, const String& output_file)
  {
    updateMembers_();

    stats_ = Statistics();
    decoy_collision_index_.clear();

    // Pass 1: Build index
    StreamingTSVReader reader;
    OPENMS_LOG_INFO << "=== Pass 1: Building index ===" << std::endl;
    auto index = reader.buildIndex(input_file);

    stats_.input_transitions = index.total_transitions;
    stats_.input_groups = index.total_groups;

    OPENMS_LOG_INFO << "Index built: " << index.collision_index.size()
                    << " unique peptides for collision detection" << std::endl;

    // Pass 2: Process groups
    StreamingTSVWriter writer;
    writer.open(output_file);

    OPENMS_LOG_INFO << "=== Pass 2: Processing groups ===" << std::endl;

    reader.processGroups(input_file, index,
      [&](const std::string& group_id,
          std::vector<OpenSwath::LightTransition>& transitions,
          const OpenSwath::LightCompound& compound,
          const std::vector<std::string>& protein_ids)
      {
        std::vector<OpenSwath::LightTransition> output_transitions;
        std::vector<OpenSwath::LightTransition> output_decoys;

        processGroup_(group_id, transitions, compound, protein_ids,
                      index.collision_index, output_transitions, output_decoys);

        if (!output_transitions.empty())
        {
          writer.writeGroup(group_id, output_transitions, output_decoys,
                            compound, protein_ids);
          stats_.output_groups++;
        }
      });

    writer.close();

    stats_.output_transitions = writer.getTransitionsWritten();

    OPENMS_LOG_INFO << "=== Processing complete ===" << std::endl;
    OPENMS_LOG_INFO << "Input: " << stats_.input_transitions << " transitions, "
                    << stats_.input_groups << " groups" << std::endl;
    OPENMS_LOG_INFO << "Output: " << stats_.output_transitions << " transitions, "
                    << stats_.output_groups << " groups" << std::endl;
    OPENMS_LOG_INFO << "Decoys: " << stats_.decoys_generated << " generated, "
                    << stats_.decoys_collided << " collisions avoided" << std::endl;
  }

  void StreamingTransitionProcessor::processGroup_(
    const std::string& /* group_id */,
    std::vector<OpenSwath::LightTransition>& transitions,
    const OpenSwath::LightCompound& compound,
    const std::vector<std::string>& /* protein_ids */,
    const std::unordered_map<std::string, std::string>& collision_index,
    std::vector<OpenSwath::LightTransition>& output_transitions,
    std::vector<OpenSwath::LightTransition>& output_decoys)
  {
    // Skip decoy transitions in input
    transitions.erase(
      std::remove_if(transitions.begin(), transitions.end(),
                     [](const OpenSwath::LightTransition& tr) { return tr.decoy; }),
      transitions.end());

    if (transitions.empty())
    {
      return;
    }

    // Filter by m/z bounds
    filterTransitions_(transitions);

    if (transitions.empty())
    {
      return;
    }

    // Select top N transitions
    selectTopTransitions_(transitions, max_transitions_);

    // Copy to output
    output_transitions = transitions;

    // Generate decoys
    if (generate_decoys_)
    {
      size_t collisions = generateDecoys_(compound, output_transitions,
                                          collision_index, output_decoys);
      stats_.decoys_collided += collisions;
    }
  }

  void StreamingTransitionProcessor::filterTransitions_(
    std::vector<OpenSwath::LightTransition>& transitions)
  {
    transitions.erase(
      std::remove_if(transitions.begin(), transitions.end(),
                     [this](const OpenSwath::LightTransition& tr) {
                       return tr.product_mz < min_product_mz_ ||
                              tr.product_mz > max_product_mz_;
                     }),
      transitions.end());
  }

  void StreamingTransitionProcessor::selectTopTransitions_(
    std::vector<OpenSwath::LightTransition>& transitions,
    size_t max_transitions)
  {
    if (transitions.size() <= max_transitions)
    {
      return;
    }

    // Sort by intensity (descending)
    std::sort(transitions.begin(), transitions.end(),
              [](const OpenSwath::LightTransition& a, const OpenSwath::LightTransition& b) {
                return a.library_intensity > b.library_intensity;
              });

    // Keep only top N
    transitions.resize(max_transitions);
  }

  std::string StreamingTransitionProcessor::pseudoReverseSequence_(const std::string& sequence)
  {
    if (sequence.size() <= 1)
    {
      return sequence;
    }

    // Keep C-terminal amino acid fixed (for tryptic peptides)
    std::string reversed = sequence.substr(0, sequence.size() - 1);
    std::reverse(reversed.begin(), reversed.end());
    return reversed + sequence.back();
  }

  bool StreamingTransitionProcessor::checkCollision_(
    const std::string& collision_key,
    const std::unordered_map<std::string, std::string>& target_index)
  {
    // Check against target index
    if (target_index.count(collision_key) > 0)
    {
      return true;
    }

    // Check against previously generated decoys
    if (decoy_collision_index_.count(collision_key) > 0)
    {
      return true;
    }

    return false;
  }

  size_t StreamingTransitionProcessor::generateDecoys_(
    const OpenSwath::LightCompound& compound,
    const std::vector<OpenSwath::LightTransition>& target_transitions,
    const std::unordered_map<std::string, std::string>& collision_index,
    std::vector<OpenSwath::LightTransition>& decoy_transitions)
  {
    size_t collisions = 0;

    // Create decoy compound
    OpenSwath::LightCompound decoy_compound = compound;
    decoy_compound.id = decoy_tag_ + compound.id;
    decoy_compound.sequence = pseudoReverseSequence_(compound.sequence);

    // Build collision key for decoy
    std::string collision_key = StreamingTSVReader::buildCollisionKey(decoy_compound);

    // Check for collision
    if (checkCollision_(collision_key, collision_index))
    {
      ++collisions;
      // Try shuffled sequence as fallback
      std::string shuffled = compound.sequence;
      if (shuffled.size() > 2)
      {
        std::shuffle(shuffled.begin(), shuffled.end() - 1, rng_);
        decoy_compound.sequence = shuffled;
        collision_key = StreamingTSVReader::buildCollisionKey(decoy_compound);

        if (checkCollision_(collision_key, collision_index))
        {
          // Still collides - skip this decoy
          return collisions;
        }
      }
      else
      {
        return collisions;
      }
    }

    // Add to decoy collision index
    decoy_collision_index_[collision_key] = decoy_compound.id;

    // Generate decoy transitions
    for (const auto& target_tr : target_transitions)
    {
      OpenSwath::LightTransition decoy_tr = target_tr;
      decoy_tr.transition_name = decoy_tag_ + target_tr.transition_name;
      decoy_tr.peptide_ref = decoy_compound.id;
      decoy_tr.decoy = true;

      // Shift product m/z
      decoy_tr.product_mz += product_mz_shift_;

      // Shift precursor m/z if configured
      if (precursor_mz_shift_ != 0.0)
      {
        decoy_tr.precursor_mz += precursor_mz_shift_;
      }

      decoy_transitions.push_back(decoy_tr);
      ++stats_.decoys_generated;
    }

    return collisions;
  }

} // namespace OpenMS
