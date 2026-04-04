// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// 
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/DIAChromHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <numeric>
#include <set>

namespace OpenMS
{

std::vector<MSChromatogram> DIAChromHandler::collectIrtChromatogramsForIrt(
  const std::vector< OpenSwath::SwathMap > & swath_maps,
  const OpenSwath::LightTargetedExperiment & irt_transitions,
  const Param & mrm_mapping_param,
  const ChromExtractParams & cp,
  const TransformationDescription& trafo,
  bool pasef,
  bool load_into_memory)
{
  // mrm_mapping_param is not used in DIA iRT collection as extraction is done directly
  (void)mrm_mapping_param;
  std::vector<MSChromatogram> chromatograms;
  TransformationDescription trafo_inverse = trafo;
  trafo_inverse.invert();

  // If this is pasef data, do chromatogram extraction beforehand in unparallel workflow
  std::vector<int> tr_win_map; // maps transition k to dia map i from which it should be extracted, only used if pasef flag is on
  if (pasef)
  {
    // Before calling this function, check to ensure that precursors actually have IM data
    for (Size k = 0; k < irt_transitions.transitions.size(); k++)
    {
      const OpenSwath::LightTransition& tr = irt_transitions.transitions[k];
      if (tr.getPrecursorIM() == -1)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Transition " + tr.getNativeID() +  " does not have a valid IM value, this must be set to use the -pasef flag");
      }
    }
    OpenSwathHelper::selectSwathTransitionsPasef(irt_transitions, tr_win_map, cp.min_upper_edge_dist, swath_maps);
  }

  // Progress reporting removed - not essential for functionality
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
  for (SignedSize map_idx = 0; map_idx < boost::numeric_cast<SignedSize>(swath_maps.size()); ++map_idx)
  {
    std::vector< OpenMS::MSChromatogram > tmp_chromatograms;
    if (!swath_maps[map_idx].ms1) // skip MS1
    {

      OpenSwath::LightTargetedExperiment transition_exp_used;

      if (pasef)
      {
        // Step 1.2: select transitions based on matching PRM/PASEF window (best window)
        std::set<std::string> matching_compounds;
        for (Size k = 0; k < tr_win_map.size(); k++)
        {
          if (tr_win_map[k] == map_idx)
          {
             const OpenSwath::LightTransition& tr = irt_transitions.transitions[k];
             transition_exp_used.transitions.push_back(tr);
             matching_compounds.insert(tr.getPeptideRef());
             OPENMS_LOG_DEBUG << "Adding Precursor with m/z " << tr.getPrecursorMZ() << " and IM of " << tr.getPrecursorIM() <<  " to swath with mz lower of " << swath_maps[map_idx].lower << " m/z upper of " << swath_maps[map_idx].upper << " im lower of " << swath_maps[map_idx].imLower << " and im upper of " << swath_maps[map_idx].imUpper << '\n';
          }
        }

        std::set<std::string> matching_proteins;
        for (Size i = 0; i < irt_transitions.compounds.size(); i++)
        {
          if (matching_compounds.find(irt_transitions.compounds[i].id) != matching_compounds.end())
          {
            transition_exp_used.compounds.push_back( irt_transitions.compounds[i] );
            for (Size j = 0; j < irt_transitions.compounds[i].protein_refs.size(); j++)
            {
              matching_proteins.insert(irt_transitions.compounds[i].protein_refs[j]);
            }
          }
        }
        for (Size i = 0; i < irt_transitions.proteins.size(); i++)
        {
          if (matching_proteins.find(irt_transitions.proteins[i].id) != matching_proteins.end())
          {
            transition_exp_used.proteins.push_back( irt_transitions.proteins[i] );
          }
        }
      }
      else
      {
        OpenSwathHelper::selectSwathTransitions(irt_transitions, transition_exp_used,
        cp.min_upper_edge_dist, swath_maps[map_idx].lower, swath_maps[map_idx].upper);
      }

      // If no transitions were selected for this swath map but the map contains
      // chromatograms only (no spectra), fall back to using all transitions
      // so that we can attempt nativeID / m/z matching against the chromatogram list.
      if (transition_exp_used.getTransitions().empty()
          && swath_maps[map_idx].sptr && swath_maps[map_idx].sptr->getNrSpectra() == 0
          && swath_maps[map_idx].sptr->getNrChromatograms() > 0)
      {
        // Fall back to all transitions for chrom-only matching when no transitions selected
        transition_exp_used = irt_transitions;
      }

      if (!transition_exp_used.getTransitions().empty()) // skip if no transitions found
      {


        std::vector< OpenSwath::ChromatogramPtr > tmp_out;
        std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
        ChromatogramExtractor extractor;

        OpenSwath::SpectrumAccessPtr current_swath_map = swath_maps[map_idx].sptr;
        if (load_into_memory)
        {
          // This creates an InMemory object that keeps all data in memory
          current_swath_map = std::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*current_swath_map) );
        }

        // Handle RT transformation like prepareExtractionCoordinates_
        if (cp.rt_extraction_window < 0)
        {
          ChromatogramExtractor::prepare_coordinates(tmp_out, coordinates, transition_exp_used, cp.rt_extraction_window, /*ms1*/ false, /*ms1_isotopes*/ 0);
        }
        else
        {
          // Use an rt extraction window of 0.0 which will just write the retention time in start / end positions
          // Then correct the start/end positions and add the extra_rt_extract parameter
          ChromatogramExtractor::prepare_coordinates(tmp_out, coordinates, transition_exp_used, 0.0, /*ms1*/ false, /*ms1_isotopes*/ 0);
          for (std::vector< ChromatogramExtractor::ExtractionCoordinates >::iterator it = coordinates.begin(); it != coordinates.end(); ++it)
          {
            it->rt_start = trafo_inverse.apply(it->rt_start) - (cp.rt_extraction_window + cp.extra_rt_extract)/ 2.0;
            it->rt_end = trafo_inverse.apply(it->rt_end) + (cp.rt_extraction_window + cp.extra_rt_extract)/ 2.0;
          }
        }

        extractor.extractChromatograms(current_swath_map, tmp_out, coordinates, cp.mz_extraction_window,
              cp.ppm, cp.im_extraction_window, cp.extraction_function);
        extractor.return_chromatogram(tmp_out, coordinates,
            transition_exp_used, SpectrumSettings(), tmp_chromatograms, false, cp.im_extraction_window);

#ifdef _OPENMP
#pragma omp critical (osw_write_chroms)
#endif
        {
          int nr_empty_chromatograms = 0;
          OPENMS_LOG_DEBUG << "[simple] Extracted "  << tmp_chromatograms.size() << " chromatograms from SWATH map " <<
            map_idx << " with m/z " << swath_maps[map_idx].lower << " to " << swath_maps[map_idx].upper << ":\n";
          for (Size chrom_idx = 0; chrom_idx < tmp_chromatograms.size(); chrom_idx++)
          {
            // Check TIC and remove empty chromatograms (can happen if the
            // extraction window is outside the mass spectrometric acquisition
            // window).
            double tic = std::accumulate(tmp_out[chrom_idx]->getIntensityArray()->data.begin(),
                                         tmp_out[chrom_idx]->getIntensityArray()->data.end(),0.0);
            OPENMS_LOG_DEBUG << "Chromatogram "  << coordinates[chrom_idx].id << " with size "
              << tmp_out[chrom_idx]->getIntensityArray()->data.size() << " and TIC " << tic  << '\n';
            if (tic > 0.0)
            {
              // add the chromatogram to the output
              chromatograms.push_back(tmp_chromatograms[chrom_idx]);
            }
            else
            {
              OPENMS_LOG_DEBUG << " - Warning: Empty chromatogram " << coordinates[chrom_idx].id <<
                " detected. Will skip it!\n";
              nr_empty_chromatograms++;
            }
          }

          if (nr_empty_chromatograms > 0)
          {
            OPENMS_LOG_WARN << " - Warning: Detected " << nr_empty_chromatograms << " empty chromatograms. Will skip them!" << std::endl;
          }
        }
      }
      else
      {
        OPENMS_LOG_DEBUG << "Extracted no transitions from SWATH map " << map_idx << " with m/z " <<
            swath_maps[map_idx].lower << " to " << swath_maps[map_idx].upper << '\n';
      }
    }
  }


  return chromatograms;
}


std::vector<MSChromatogram> DIAChromHandler::extractAndMapChromatogramsForTransitions(
  const std::vector< OpenSwath::SwathMap > & swath_maps,
  const OpenSwath::LightTargetedExperiment & transition_exp,
  const ChromExtractParams & cp,
  const Param & mrm_mapping_param)
{
  // mrm_mapping_param is not used in DIA handler as ChromatogramExtractor assigns native IDs directly
  (void)mrm_mapping_param;
  std::vector<MSChromatogram> collected;
  ChromatogramExtractor extractor;

  for (Size i = 0; i < swath_maps.size(); ++i)
  {
    if (swath_maps[i].ms1) continue;
    OpenSwath::SpectrumAccessPtr current = swath_maps[i].sptr;
    if (!current) continue;

  std::vector< OpenSwath::ChromatogramPtr > chrom_list;
  std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
  const OpenSwath::LightTargetedExperiment& transition_exp_copy = transition_exp;
  ChromatogramExtractor::prepare_coordinates(chrom_list, coordinates, transition_exp_copy, cp.rt_extraction_window, /*ms1*/ false, /*ms1_isotopes*/ 0);

    try
    {
      extractor.extractChromatograms(current, chrom_list, coordinates, cp.mz_extraction_window, cp.ppm, cp.im_extraction_window, cp.extraction_function);
      PeakMap pm;
      extractor.return_chromatogram(chrom_list, coordinates, transition_exp, SpectrumSettings(), pm.getChromatograms(), false, cp.im_extraction_window);

  for (auto & c : pm.getChromatograms()) collected.push_back(c);
  OPENMS_LOG_DEBUG << "DIAChromHandler: extracted " << pm.getChromatograms().size() << " chromatograms from swath " << i << " for transitions" << std::endl;
    }
    catch (const std::exception & e)
    {
      OPENMS_LOG_DEBUG << "DIAChromHandler: extraction failed for swath " << i << " : " << e.what() << std::endl;
      continue;
    }
  }

  // ChromatogramExtractor::return_chromatogram already populated the
  // extracted chromatograms with native IDs corresponding to the transition
  // native IDs (coord.id). Therefore, we can directly return the collected
  // chromatograms without additional mapping.
  OPENMS_LOG_DEBUG << "DIAChromHandler: returning " << collected.size() << " extracted chromatograms for transitions" << std::endl;
  return collected;
}

} // namespace OpenMS
