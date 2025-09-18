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

#define IMPLIES(a, b) !(a) || (b)

namespace OpenMS
{
  template <typename MapT, typename PepT>
  bool populateMS1Transition(MapT & pep2tr,
                             const PepT & pep,
                             ChromatogramExtractor::ExtractionCoordinates & coord)
  {
    // default values for RT window (negative range)
    coord.rt_end = -1;
    coord.rt_start = 0;

    // Catch cases where a compound has no transitions
    if (pep2tr.count(pep.id) == 0)
    {
      OPENMS_LOG_INFO << "Warning: no transitions found for compound " << pep.id << std::endl;
      coord.id = OpenSwathHelper::computePrecursorId(pep.id, 0);
      return false;
    }

    // This is slightly awkward but the m/z of the precursor is *not*
    // stored in the precursor object but only in the transition object
    // itself. So we have to get the first transition to look it up.
    auto transition = (*pep2tr[pep.id][0]);
    coord.mz = transition.getPrecursorMZ();

    // Set chromatogram reference id: even though we use the peptide id
    // here, it is possible that these ids overlap with the transition
    // ids, leading to bad downstream consequences (e.g. ambiguity which
    // chromatograms are precursor and which ones are fragment
    // chromatograms). This is especially problematic with pqp files
    // where peptide precursors and transitions are simply numbered and
    // are guaranteed to overlap.
    coord.id = OpenSwathHelper::computePrecursorId(pep.id, 0);
    return true;
  }

  template <typename TransitionT>
  void populateMS2Transition(const TransitionT & transition,
                             ChromatogramExtractor::ExtractionCoordinates & coord)
  {
    // default values for RT window (negative range)
    coord.rt_end = -1;
    coord.rt_start = 0;

    coord.mz = transition.getProductMZ();
    coord.mz_precursor = transition.getPrecursorMZ();
    coord.id = transition.getNativeID();
  }


  const TargetedExperimentHelper::PeptideCompound* getPeptideHelperMS2_(const OpenMS::TargetedExperiment& transition_exp_used,
                                                                        const OpenMS::ReactionMonitoringTransition& transition,
                                                                        bool do_peptides)
  {
    OPENMS_PRECONDITION(IMPLIES(do_peptides, !transition.getPeptideRef().empty()), "PeptideRef cannot be empty for peptides")
    OPENMS_PRECONDITION(IMPLIES(!do_peptides, !transition.getCompoundRef().empty()), "CompoundRef cannot be empty for compounds")

    if (do_peptides)
    {
      return &transition_exp_used.getPeptideByRef(transition.getPeptideRef()); 
    }
    else
    {
      return &transition_exp_used.getCompoundByRef(transition.getCompoundRef()); 
    }
  }

  const TargetedExperimentHelper::PeptideCompound* getPeptideHelperMS1_(const OpenMS::TargetedExperiment & transition_exp_used,
                                                                        Size i,
                                                                        bool do_peptides)
  {
    OPENMS_PRECONDITION(IMPLIES(do_peptides, i < transition_exp_used.getPeptides().size()), "Index i must be smaller than the number of peptides")
    OPENMS_PRECONDITION(IMPLIES(!do_peptides, i < transition_exp_used.getCompounds().size()), "Index i must be smaller than the number of compounds")

    if (do_peptides)
    {
      return &transition_exp_used.getPeptides()[i];
    }
    else
    {
      return &transition_exp_used.getCompounds()[i];
    }
  }

  void ChromatogramExtractor::prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
                                                  std::vector< ExtractionCoordinates > & coordinates,
                                                  const OpenSwath::LightTargetedExperiment & transition_exp_used,
                                                  const double rt_extraction_window,
                                                  const bool ms1,
                                                  const int ms1_isotopes)
  {
    // hash of the peptide reference containing all transitions
    std::map<String, std::vector<const OpenSwath::LightTransition*> > pep2tr;
    for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
    {
      String ref = transition_exp_used.getTransitions()[i].getPeptideRef();
      pep2tr[ref].push_back(&transition_exp_used.getTransitions()[i]);
    }
    std::map<String, const OpenSwath::LightCompound*> tr2pep;
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
        // if 'rt_extraction_window' is non-zero, just use the (first) RT value
        double rt = pep.rt;
        coord.rt_start = rt - rt_extraction_window / 2.0;
        coord.rt_end = rt + rt_extraction_window / 2.0;
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

  void ChromatogramExtractor::prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
                                                  std::vector< ExtractionCoordinates > & coordinates,
                                                  const OpenMS::TargetedExperiment & transition_exp_used,
                                                  const double rt_extraction_window,
                                                  const bool ms1,
                                                  const int ms1_isotopes)
  {
    // hash of the peptide reference containing all transitions
    typedef std::map<String, std::vector<const ReactionMonitoringTransition*> > PeptideTransitionMapType;
    PeptideTransitionMapType pep2tr;
    for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
    {
      String ref = transition_exp_used.getTransitions()[i].getPeptideRef();
      if (ref.empty()) ref = transition_exp_used.getTransitions()[i].getCompoundRef();
      pep2tr[ref].push_back(&transition_exp_used.getTransitions()[i]);
    }

    // std::map<String, const TargetedExperimentHelper::PeptideCompound* > tr2pep;
    // for (const auto & p : transition_exp_used.getPeptides()) {tr2pep[p.id] = &p;}
    // for (const auto & c : transition_exp_used.getCompounds()) {tr2pep[c.id] = &c;}

    bool have_peptides = (!transition_exp_used.getPeptides().empty());

    // Determine iteration size (nr peptides or nr transitions)
    Size itersize;
    if (ms1)
    {
      if (have_peptides)
      {
        itersize = transition_exp_used.getPeptides().size();
      }
      else
      {
        itersize = transition_exp_used.getCompounds().size();
      }
    }
    else
    {
      itersize = transition_exp_used.getTransitions().size();
    }

    for (Size i = 0; i < itersize; i++)
    {
      OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
      output_chromatograms.push_back(s);

      ChromatogramExtractor::ExtractionCoordinates coord;
      const TargetedExperimentHelper::PeptideCompound* pep;
      OpenMS::ReactionMonitoringTransition transition;

      if (ms1) 
      {
        pep = getPeptideHelperMS1_(transition_exp_used, i, have_peptides);
        if (!populateMS1Transition(pep2tr, *pep, coord))
        {
          // Catch cases where a compound has no transitions
          coordinates.push_back(coord);
          continue;
        }
      }
      else
      {
        transition = transition_exp_used.getTransitions()[i];
        pep = getPeptideHelperMS2_(transition_exp_used, transition, have_peptides);
        populateMS2Transition(transition, coord);
      }

      if (rt_extraction_window < 0) {} // construct for NAN (see below)
      else
      {
        if (!pep->hasRetentionTime())
        {
          // we don't have retention times -> this is only a problem if we actually
          // wanted to use the RT limit feature.
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "Error: Peptide " + pep->id + " does not have retention time information which is necessary to perform an RT-limited extraction");
        }
        else if (std::isnan(rt_extraction_window)) // if 'rt_extraction_window' is NAN, we assume that RT start/end is encoded in the data
        {
          // TODO: better use a single RT entry with start/end
          if (pep->rts.size() != 2)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "Error: Expected exactly two retention time entries for peptide '" + pep->id + "', found " + String(pep->rts.size()));
          }
          coord.rt_start = pep->rts[0].getRT();
          coord.rt_end = pep->rts[1].getRT();
        }
        else // if 'rt_extraction_window' is non-zero, just use the (first) RT value
        {
          double rt = pep->getRetentionTime();
          coord.rt_start = rt - rt_extraction_window / 2.0;
          coord.rt_end = rt + rt_extraction_window / 2.0;
        }
      }
      coord.ion_mobility = pep->getDriftTime();
      coordinates.push_back(coord);

      if (ms1 && ms1_isotopes > 0 && false)
      {
        for (int k = 1; k <= ms1_isotopes; k++)
        {
          OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
          output_chromatograms.push_back(s);
          ChromatogramExtractor::ExtractionCoordinates coord_new = coord;
          coord_new.id = OpenSwathHelper::computePrecursorId(pep->id, k);
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
