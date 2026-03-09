// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>

#include <OpenMS/CONCEPT/LogStream.h>
// LightTargetedExperiment is defined in openswathalgo
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

using namespace std;

namespace OpenMS
{

  MRMMapping::MRMMapping() :
      DefaultParamHandler("MRMMapping")
    {
      defaults_.setValue("precursor_tolerance", 0.1, "Precursor tolerance when mapping (in Th)");
      defaults_.setValue("product_tolerance", 0.1, "Product tolerance when mapping (in Th)");

      defaults_.setValue("map_multiple_assays", "false", "Allow to map multiple assays to chromatograms and duplicate these chromatograms in the output.");
      defaults_.setValidStrings("map_multiple_assays", {"true","false"});

      defaults_.setValue("error_on_unmapped", "false", "Treat remaining, unmapped chromatograms as an error");
      defaults_.setValidStrings("error_on_unmapped", {"true","false"});

      // write defaults into Param object param_
      defaultsToParam_();
      updateMembers_();
    }

  void MRMMapping::updateMembers_()
  {
    precursor_tol_ = (double)param_.getValue("precursor_tolerance");
    product_tol_ = (double)param_.getValue("product_tolerance");
    map_multiple_assays_ = (bool)param_.getValue("map_multiple_assays").toBool();
    error_on_unmapped_ = (bool)param_.getValue("error_on_unmapped").toBool();
  }

  void MRMMapping::mapExperiment(const OpenMS::PeakMap& chromatogram_map, 
      const OpenMS::TargetedExperiment& targeted_exp,
      OpenMS::PeakMap& output) const
  {
    // copy all meta data from old MSExperiment
    output = (ExperimentalSettings)chromatogram_map;
    output.clear(false);
    std::vector<MSChromatogram > empty_chromats;
    output.setChromatograms(empty_chromats);

    int notmapped = 0;
    int multi_mapped = 0;
    for (Size i = 0; i < chromatogram_map.getChromatograms().size(); i++)
    {
      // try to find the best matching transition for this chromatogram
      const MSChromatogram& chromatogram = chromatogram_map.getChromatograms()[i];

      bool prec_product_set = !( std::fabs(chromatogram.getPrecursor().getMZ()) < 1e-5 && 
                                 std::fabs(chromatogram.getProduct().getMZ()) < 1e-5);
      if (!prec_product_set)
      {
        if (map_multiple_assays_)
        {
          OPENMS_LOG_DEBUG << "Chromatogram " + String(chromatogram.getNativeID()) +
            " has no precursor or product m/z recorded, mapping may not work." << std::endl;
        }
        else
        {
          OPENMS_LOG_DEBUG << "Skip mapping for chromatogram " +
            String(chromatogram.getNativeID()) + " since no precursor or product m/z was recorded." << std::endl;
          continue;
        }
      }

      std::vector<MSChromatogram > mapped_chroms;
      for (Size j = 0; j < targeted_exp.getTransitions().size(); j++)
      {
        if (fabs(chromatogram.getPrecursor().getMZ() - targeted_exp.getTransitions()[j].getPrecursorMZ()) < precursor_tol_ &&
            fabs(chromatogram.getProduct().getMZ()   - targeted_exp.getTransitions()[j].getProductMZ())   < product_tol_)
        {
          OPENMS_LOG_DEBUG << "Mapping chromatogram " << i << " to transition " << j << " (" << targeted_exp.getTransitions()[j].getNativeID() << ")"
             " with precursor mz " << chromatogram.getPrecursor().getMZ() << " / " <<  targeted_exp.getTransitions()[j].getPrecursorMZ() <<
             " and product mz " << chromatogram.getProduct().getMZ() << " / " <<  targeted_exp.getTransitions()[j].getProductMZ() << std::endl;

          // Create precursor and set the peptide sequence
          MSChromatogram c = chromatogram_map.getChromatograms()[i];
          Precursor precursor = c.getPrecursor();
          String pepref = targeted_exp.getTransitions()[j].getPeptideRef();
          precursor.setMetaValue("peptide_sequence", pepref);
          precursor.setMetaValue("description", targeted_exp.getTransitions()[j].getNativeID());
          for (Size pep_idx = 0; pep_idx < targeted_exp.getPeptides().size(); pep_idx++)
          {
            const OpenMS::TargetedExperiment::Peptide * pep = &targeted_exp.getPeptides()[pep_idx];
            if (pep->id == pepref)
            {
              if (!pep->sequence.empty())
              {
                precursor.setMetaValue("peptide_sequence", pep->sequence);
              }
              break;
            }
          }
          // add precursor to chromatogram
          c.setPrecursor(precursor);

          // Set the id of the chromatogram, using the id of the transition (this gives directly the mapping of the two
          c.setNativeID(targeted_exp.getTransitions()[j].getNativeID());

          mapped_chroms.push_back(c);
        }
      }

      // Check whether we have mapped this chromatogram to at least one transition:
      //  - warn if no mapping occurred
      //  - else append all mapped chromatograms (if we allow multiple mappings)
      //  - else append the first mapped chromatograms (if we don't allow multiple mappings)
      if (mapped_chroms.empty())
      {
        OPENMS_LOG_DEBUG << "Did not find a mapping for chromatogram " + String(i) + " with transition " + String(chromatogram.getPrecursor().getMZ()) + \
          " -> " + String(chromatogram.getProduct().getMZ()) +  "! Maybe try to increase your mapping tolerance." << std::endl;
        notmapped++;
        if (error_on_unmapped_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Did not find a mapping for chromatogram " + String(i) + "! Maybe try to increase your mapping tolerance.");
        }
      }
      else if (map_multiple_assays_)
      {
        for (auto & c : mapped_chroms) output.addChromatogram(c);
        if (mapped_chroms.size() > 1)
        {
          ++multi_mapped;
          OPENMS_LOG_DEBUG << "Chromatogram " + String(chromatogram.getNativeID()) <<
            " with " + String(chromatogram.getPrecursor().getMZ()) <<
            " -> " + String(chromatogram.getProduct().getMZ()) <<
            " maps to multiple assays!" << std::endl;
        }
      }
      else
      {
        if (mapped_chroms.size() == 1) output.addChromatogram(mapped_chroms[0]);
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Chromatogram " + String(chromatogram.getNativeID()) + \
           " with " + String(chromatogram.getPrecursor().getMZ()) + \
            " -> " + String(chromatogram.getProduct().getMZ()) + \
              " maps to multiple assays! Either decrease your mapping tolerance or set map_multiple_assays to true.");
        }
      }
    }

    int total = static_cast<int>(chromatogram_map.getChromatograms().size());
    int mapped = total - notmapped;
    OPENMS_LOG_DEBUG << "MRMMapping: mapped=" << mapped << ", unmapped=" << notmapped << ", multi_mapped=" << multi_mapped << " chromatogram(s)." << std::endl;
    if (notmapped > 0 && error_on_unmapped_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Found " + String(notmapped) + \
          " unmapped chromatograms, disable error_on_unmapped to continue.");
    }


  }


void MRMMapping::mapExperiment(const OpenMS::PeakMap& chromatogram_map,
    const OpenSwath::LightTargetedExperiment& targeted_exp,
    OpenMS::PeakMap& output) const
{
  // Convert LightTargetedExperiment -> TargetedExperiment (minimal conversion)
  OpenMS::TargetedExperiment te;

  // Convert compounds/peptides
  for (const auto & lc : targeted_exp.getCompounds())
  {
    OpenMS::TargetedExperiment::Peptide p;
    p.id = lc.id;
    p.sequence = lc.sequence;
    te.addPeptide(p);
  }

  // Convert proteins
  for (const auto & lp : targeted_exp.getProteins())
  {
    OpenMS::TargetedExperiment::Protein prot;
    prot.id = lp.id;
    prot.sequence = lp.sequence;
    te.addProtein(prot);
  }

  // Convert transitions
  for (const auto & lt : targeted_exp.getTransitions())
  {
    OpenMS::ReactionMonitoringTransition tr;
    tr.setNativeID(lt.getNativeID());
    tr.setPeptideRef(lt.getPeptideRef());
    tr.setPrecursorMZ(lt.precursor_mz);
    tr.setProductMZ(lt.product_mz);
    tr.setLibraryIntensity(lt.library_intensity);
    // Flags: detecting / quantifying / identifying -> map to transition flags
    tr.setDetectingTransition(lt.flags.detecting != 0);
    tr.setQuantifyingTransition(lt.flags.quantifying != 0);
    tr.setIdentifyingTransition(lt.flags.identifying != 0);
    // Decoy flag
    if (lt.getDecoy()) tr.setDecoyTransitionType(OpenMS::ReactionMonitoringTransition::DECOY);
    else tr.setDecoyTransitionType(OpenMS::ReactionMonitoringTransition::TARGET);

    te.addTransition(tr);
  }

  // Now call the existing implementation
  mapExperiment(chromatogram_map, te, output);

}

} //namespace

