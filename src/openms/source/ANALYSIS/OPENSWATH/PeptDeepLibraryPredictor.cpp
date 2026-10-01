// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeptDeepLibraryPredictor.h>

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h>

#include <cstdint>
#include <unordered_set>
#include <utility>

namespace OpenMS
{
  namespace
  {
    constexpr Size PEPTDEEP_MS2_CHANNELS = 8;
    constexpr const char* DEFAULT_RT_MODEL = "peptdeep_rt_dynamic.onnx";
    constexpr const char* DEFAULT_CCS_MODEL = "peptdeep_ccs_dynamic.onnx";
    constexpr const char* DEFAULT_MS2_MODEL = "peptdeep_ms2_dynamic.onnx";

    std::string resolveModelPath_(const std::string& configured_path, const char* model_name)
    {
      if (!configured_path.empty())
      {
        return configured_path;
      }

      // PeptDeep models are downloaded into the build tree's share/OpenMS/models
      // directory and installed into the normal OpenMS data directory. Search the
      // build-tree data directory first, then let File::find() fall back to the
      // resolved installed/source OpenMS data path.
      return File::find(
        std::string("models/") + model_name,
        {std::string(OPENMS_BINARY_PATH) + "/share/OpenMS"});
    }

    void addModifications_(const AASequence& peptide, OpenSwath::LightCompound& compound)
    {
      if (peptide.hasNTerminalModification())
      {
        compound.modifications.push_back({-1, peptide.getNTerminalModification()->getUniModRecordId()});
      }

      for (Size i = 0; i < peptide.size(); ++i)
      {
        if (peptide[i].isModified())
        {
          compound.modifications.push_back(
            {static_cast<int>(i), peptide.getResidue(i).getModification()->getUniModRecordId()});
        }
      }

      if (peptide.hasCTerminalModification())
      {
        compound.modifications.push_back(
          {static_cast<int>(peptide.size()), peptide.getCTerminalModification()->getUniModRecordId()});
      }
    }
  } // namespace

  PeptDeepLibraryPredictor::PeptDeepLibraryPredictor() : PeptDeepLibraryPredictor(Config{})
  {
  }

  PeptDeepLibraryPredictor::PeptDeepLibraryPredictor(const Config& config) :
    rt_predictor_(resolveModelPath_(config.rt_model_path, DEFAULT_RT_MODEL), config.intra_op_threads, config.batch_size),
    ms2_predictor_(resolveModelPath_(config.ms2_model_path, DEFAULT_MS2_MODEL), config.intra_op_threads, config.batch_size)
  {
    if (config.predict_ccs)
    {
      ccs_predictor_ = std::make_unique<PeptDeepCCSInference>(
        resolveModelPath_(config.ccs_model_path, DEFAULT_CCS_MODEL), config.intra_op_threads, config.batch_size);
    }
  }

  OpenSwath::LightTargetedExperiment PeptDeepLibraryPredictor::predict(
    const std::vector<PeptDeepLibraryPrecursor>& precursors)
  {
    if (precursors.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor batch cannot be empty.");
    }

    std::vector<std::string> peptide_strings;
    std::vector<float> charges;
    std::vector<float> nces;
    std::vector<int64_t> instrument_indices;
    peptide_strings.reserve(precursors.size());
    charges.reserve(precursors.size());
    nces.reserve(precursors.size());
    instrument_indices.reserve(precursors.size());

    Size transition_count = 0;
    std::unordered_set<std::string> precursor_ids;
    precursor_ids.reserve(precursors.size());

    for (const auto& precursor : precursors)
    {
      if (precursor.id.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor ID cannot be empty.");
      }
      if (!precursor_ids.insert(precursor.id).second)
      {
        throw Exception::IllegalArgument(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor IDs must be unique within a prediction batch: " + precursor.id);
      }
      if (precursor.charge <= 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor charge must be positive.");
      }
      if (precursor.peptide.size() < 2)
      {
        throw Exception::IllegalArgument(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide sequences must contain at least two residues.");
      }
      for (const auto& protein_ref : precursor.protein_refs)
      {
        if (protein_ref.empty())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Protein IDs cannot be empty.");
        }
      }

      peptide_strings.push_back(precursor.peptide.toString());
      charges.push_back(static_cast<float>(precursor.charge));
      nces.push_back(precursor.nce);
      instrument_indices.push_back(static_cast<int64_t>(precursor.instrument_index));
      transition_count += (precursor.peptide.size() - 1) * 4;
    }

    const std::vector<float> predicted_rt = rt_predictor_.predictRT(peptide_strings);
    const std::vector<std::vector<float>> predicted_ms2 =
      ms2_predictor_.predictMS2(peptide_strings, charges, nces, instrument_indices);

    std::vector<float> predicted_ccs;
    if (ccs_predictor_)
    {
      predicted_ccs = ccs_predictor_->predictCCS(peptide_strings, charges);
    }

    if (predicted_rt.size() != precursors.size() || predicted_ms2.size() != precursors.size() ||
        (ccs_predictor_ && predicted_ccs.size() != precursors.size()))
    {
      throw Exception::InvalidValue(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "prediction count mismatch", "PeptDeep returned an unexpected number of predictions.");
    }

    OpenSwath::LightTargetedExperiment experiment;
    experiment.compounds.reserve(precursors.size());
    experiment.transitions.reserve(transition_count);

    std::unordered_set<std::string> protein_ids;
    for (Size precursor_index = 0; precursor_index < precursors.size(); ++precursor_index)
    {
      const auto& precursor = precursors[precursor_index];
      const Size cleavage_count = precursor.peptide.size() - 1;
      const Size expected_ms2_values = cleavage_count * PEPTDEEP_MS2_CHANNELS;
      if (predicted_ms2[precursor_index].size() != expected_ms2_values)
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          std::to_string(predicted_ms2[precursor_index].size()),
          "Expected PeptDeep MS2 output size " + std::to_string(expected_ms2_values) + " for precursor " + precursor.id + ".");
      }

      const double precursor_mz = precursor.peptide.getMZ(precursor.charge);
      double precursor_im = -1.0;
      if (ccs_predictor_)
      {
        precursor_im = IMTypes::ccsToOneOverK0(predicted_ccs[precursor_index], precursor_mz, precursor.charge);
      }

      OpenSwath::LightCompound compound;
      compound.id = precursor.id;
      compound.sequence = precursor.peptide.toUniModString();
      compound.charge = precursor.charge;
      compound.rt = predicted_rt[precursor_index];
      compound.drift_time = precursor_im;
      compound.protein_refs = precursor.protein_refs;
      addModifications_(precursor.peptide, compound);
      experiment.compounds.push_back(std::move(compound));

      for (const auto& protein_ref : precursor.protein_refs)
      {
        if (protein_ids.insert(protein_ref).second)
        {
          OpenSwath::LightProtein protein;
          protein.id = protein_ref;
          experiment.proteins.push_back(std::move(protein));
        }
      }

      const auto& ms2 = predicted_ms2[precursor_index];
      for (Size row = 0; row < cleavage_count; ++row)
      {
        const Size cleavage = row + 1;
        const Size y_ordinal = precursor.peptide.size() - cleavage;
        const AASequence b_fragment = precursor.peptide.getPrefix(cleavage);
        const AASequence y_fragment = precursor.peptide.getSuffix(y_ordinal);
        const Size offset = row * PEPTDEEP_MS2_CHANNELS;

        auto add_transition = [&](OpenSwath::FragmentIonType fragment_type,
                                  Size ordinal,
                                  Int fragment_charge,
                                  double product_mz,
                                  float intensity)
        {
          OpenSwath::LightTransition transition;
          transition.transition_name = precursor.id + "_" + OpenSwath::fragmentIonTypeToString(fragment_type) +
                                       std::to_string(ordinal) + "^" + std::to_string(fragment_charge);
          transition.peptide_ref = precursor.id;
          transition.library_intensity = intensity;
          transition.product_mz = product_mz;
          transition.precursor_mz = precursor_mz;
          transition.precursor_im = precursor_im;
          transition.fragment_charge = static_cast<int8_t>(fragment_charge);
          transition.fragment_nr = static_cast<int16_t>(ordinal);
          transition.fragment_type = fragment_type;
          experiment.transitions.push_back(std::move(transition));
        };

        add_transition(OpenSwath::FragmentIonType::BIon, cleavage, 1, b_fragment.getMZ(1, Residue::BIon), ms2[offset]);
        add_transition(OpenSwath::FragmentIonType::BIon, cleavage, 2, b_fragment.getMZ(2, Residue::BIon), ms2[offset + 1]);
        add_transition(OpenSwath::FragmentIonType::YIon, y_ordinal, 1, y_fragment.getMZ(1, Residue::YIon), ms2[offset + 2]);
        add_transition(OpenSwath::FragmentIonType::YIon, y_ordinal, 2, y_fragment.getMZ(2, Residue::YIon), ms2[offset + 3]);
      }
    }

    return experiment;
  }
} // namespace OpenMS
