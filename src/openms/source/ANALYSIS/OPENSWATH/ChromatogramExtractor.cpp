// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <unordered_map>

namespace OpenMS
{
  // Populate the precursor (MS1) extraction coordinate for one compound.
  // Returns false if the compound has no associated transitions, in which case
  // the caller should record an empty coordinate and continue.
  bool populateMS1Transition(const std::unordered_map<String, std::vector<const OpenSwath::LightTransition*>>& pep2tr,
                             const OpenSwath::LightCompound & pep,
                             ChromatogramExtractor::ExtractionCoordinates & coord)
  {
    // default values for RT window (negative range)
    coord.rt_end = -1;
    coord.rt_start = 0;

    // Catch cases where a compound has no transitions
    auto it = pep2tr.find(pep.id);
    if (it == pep2tr.end())
    {
      OPENMS_LOG_INFO << "Warning: no transitions found for compound " << pep.id << std::endl;
      coord.id = OpenSwathHelper::computePrecursorId(pep.id, 0);
      return false;
    }

    // The m/z of the precursor is *not* stored in the precursor object but only
    // in the transition object itself. We get the first transition to look it up.
    coord.mz = it->second[0]->getPrecursorMZ();

    // Set chromatogram reference id: even though we use the peptide id here,
    // it is possible that these ids overlap with the transition ids, leading
    // to bad downstream consequences (e.g. ambiguity which chromatograms are
    // precursor and which ones are fragment chromatograms). This is especially
    // problematic with pqp files where peptide precursors and transitions are
    // simply numbered and are guaranteed to overlap.
    coord.id = OpenSwathHelper::computePrecursorId(pep.id, 0);
    return true;
  }

  // Populate the fragment (MS2) extraction coordinate from a single transition.
  void populateMS2Transition(const OpenSwath::LightTransition & transition,
                             ChromatogramExtractor::ExtractionCoordinates & coord)
  {
    // default values for RT window (negative range)
    coord.rt_end = -1;
    coord.rt_start = 0;

    coord.mz = transition.getProductMZ();
    coord.mz_precursor = transition.getPrecursorMZ();
    coord.id = transition.getNativeID();
  }

  void ChromatogramExtractor::prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
                                                  std::vector< ExtractionCoordinates > & coordinates,
                                                  const OpenSwath::LightTargetedExperiment & transition_exp_used,
                                                  const double rt_extraction_window,
                                                  const bool ms1,
                                                  const int ms1_isotopes)
  {
    std::unordered_map<String, std::vector<const OpenSwath::LightTransition*> > pep2tr;
    for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
    {
      String ref = transition_exp_used.getTransitions()[i].getPeptideRef();
      pep2tr[ref].push_back(&transition_exp_used.getTransitions()[i]);
    }
    std::unordered_map<String, const OpenSwath::LightCompound*> tr2pep;
    tr2pep.reserve(transition_exp_used.getCompounds().size());
    for (const auto & p : transition_exp_used.getCompounds()) {tr2pep[p.id] = &p;}

    // Determine iteration size:
    // When extracting MS1/precursor transitions, we iterate over compounds.
    // Otherwise (for SWATH/fragment ions), we iterate over the transitions.
    Size itersize;
    if (ms1)
    {
      itersize = transition_exp_used.getCompounds().size();
    }
    else
    {
      itersize = transition_exp_used.getTransitions().size();
    }

    // Pre-allocate vectors for better performance
    Size expected_size = itersize * (1 + (ms1 && ms1_isotopes > 0 ? ms1_isotopes : 0));
    output_chromatograms.reserve(output_chromatograms.size() + expected_size);
    coordinates.reserve(coordinates.size() + expected_size);

    for (Size i = 0; i < itersize; i++)
    {
      OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
      output_chromatograms.push_back(s);

      ChromatogramExtractor::ExtractionCoordinates coord;
      OpenSwath::LightCompound pep;
      OpenSwath::LightTransition transition;

      if (ms1)
      {
        pep = transition_exp_used.getCompounds()[i];
        if (!populateMS1Transition(pep2tr, pep, coord))
        {
          // Catch cases where a compound has no transitions
          coordinates.push_back(coord);
          continue;
        }
      }
      else
      {
        transition = transition_exp_used.getTransitions()[i];
        pep = (*tr2pep[transition.getPeptideRef()]);
        populateMS2Transition(transition, coord);
      }

      if (rt_extraction_window >= 0)
      {
        // Non-negative rt_extraction_window: use pep.rt as window center.
        // pep.rt is NaN when the source library had no retention-time info;
        // throw to match the pre-unification heavyweight overload's behavior
        // (silently using NaN here would produce an empty extraction window).
        if (std::isnan(pep.rt))
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error: Peptide " + pep.id +
              " does not have retention time information which is necessary to perform an RT-limited extraction");
        }
        double rt = pep.rt;
        coord.rt_start = rt - rt_extraction_window / 2.0;
        coord.rt_end = rt + rt_extraction_window / 2.0;
      }
      else if (std::isnan(rt_extraction_window))
      {
        // NaN means: use the per-compound RT range encoded on the LightCompound
        // (populated from a heavyweight TargetedExperiment::Peptide/Compound whose
        // `rts` vector had exactly two entries — see
        // OpenSwathDataAccessHelper::convertTargetedExp).
        if (std::isnan(pep.rt_start) || std::isnan(pep.rt_end))
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error: rt_extraction_window=NaN requires the compound '" + pep.id +
              "' to carry a retention-time range (rt_start, rt_end), but none is set.");
        }
        coord.rt_start = pep.rt_start;
        coord.rt_end   = pep.rt_end;
      }
      coord.ion_mobility = pep.getDriftTime();
      coordinates.push_back(coord);

      if (ms1 && ms1_isotopes > 0)
      {
        for (int k = 1; k <= ms1_isotopes; k++)
        {
          OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
          output_chromatograms.push_back(s);
          ChromatogramExtractor::ExtractionCoordinates coord_new = coord;
          coord_new.id = OpenSwathHelper::computePrecursorId(pep.id, k);
          coord_new.mz = coord.mz + k * Constants::C13C12_MASSDIFF_U;
          coordinates.push_back(coord_new);
        }
      }
    }

    // sort result, use stable_sort to ensure that ordering is preserved
    std::stable_sort(coordinates.begin(), coordinates.end(), ChromatogramExtractor::ExtractionCoordinates::SortExtractionCoordinatesByMZ);
  }

  bool ChromatogramExtractor::outsideExtractionWindow_(const ReactionMonitoringTransition& transition, double current_rt,
                                 const TransformationDescription& trafo, double rt_extraction_window)
  {
    if (rt_extraction_window < 0)
    {
      return false;
    }

    // Get the expected retention time, apply the RT-transformation
    // (which describes the normalization) and then take the difference.
    // Note that we inverted the transformation in the beginning because
    // we want to transform from normalized to real RTs here and not the
    // other way round.
    double expected_rt = PeptideRTMap_[transition.getPeptideRef()];
    double de_normalized_experimental_rt = trafo.apply(expected_rt);
    if (current_rt < de_normalized_experimental_rt - rt_extraction_window / 2.0 ||
        current_rt > de_normalized_experimental_rt + rt_extraction_window / 2.0 )
    {
      return true;
    }
    return false;
  }

  int ChromatogramExtractor::getFilterNr_(const String& filter)
  {
    if (filter == "tophat")
    {
      return 1;
    }
    else if (filter == "bartlett")
    {
      return 2;
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Filter either needs to be tophat or bartlett");
    }
  }

  void ChromatogramExtractor::populatePeptideRTMap_(OpenMS::TargetedExperiment& transition_exp, double rt_extraction_window)
  {
      // Store the peptide retention times in an intermediate map
      PeptideRTMap_.clear();
      for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
      {
        const TargetedExperiment::Peptide& pep = transition_exp.getPeptides()[i];
        if (!pep.hasRetentionTime())
        {
          // we don't have retention times -> this is only a problem if we actually
          // wanted to use the RT limit feature.
          if (rt_extraction_window >= 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                             "Error: Peptide " + pep.id + " does not have retention time information which is necessary to perform an RT-limited extraction");
          }
          continue;
        }
        PeptideRTMap_[pep.id] = pep.getRetentionTime();
      }
  }

}
