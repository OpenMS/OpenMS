// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>

namespace OpenMS
{

  void OpenSwathDataAccessHelper::convertToOpenMSSpectrum(const OpenSwath::SpectrumPtr& sptr, OpenMS::MSSpectrum & spectrum)
  {
    std::vector<double>::const_iterator mz_it = sptr->getMZArray()->data.begin();
    std::vector<double>::const_iterator int_it = sptr->getIntensityArray()->data.begin();

    if (!spectrum.empty()) spectrum.clear(false);

    Peak1D p;
    spectrum.reserve(sptr->getMZArray()->data.size());
    for (; mz_it != sptr->getMZArray()->data.end(); ++mz_it, ++int_it)
    {
      p.setMZ(*mz_it);
      p.setIntensity(*int_it);
      spectrum.push_back(p);
    }
  }

  OpenSwath::SpectrumPtr OpenSwathDataAccessHelper::convertToSpectrumPtr(const OpenMS::MSSpectrum & spectrum)
  {
    OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
    OpenSwath::BinaryDataArrayPtr intensity_array = sptr->getIntensityArray();
    OpenSwath::BinaryDataArrayPtr mz_array = sptr->getMZArray();
    mz_array->data.reserve(spectrum.size());
    intensity_array->data.reserve(spectrum.size());
    for (MSSpectrum::const_iterator it = spectrum.begin(); it != spectrum.end(); ++it)
    {
      mz_array->data.push_back(it->getMZ());
      intensity_array->data.push_back(it->getIntensity());
    }
    return sptr;
  }

  OpenSwath::ChromatogramPtr OpenSwathDataAccessHelper::convertToChromatogramPtr(const OpenMS::MSChromatogram & chromatogram)
  {
    OpenSwath::ChromatogramPtr cptr(new OpenSwath::Chromatogram);
    OpenSwath::BinaryDataArrayPtr intensity_array = cptr->getIntensityArray();
    OpenSwath::BinaryDataArrayPtr rt_array = cptr->getTimeArray();
    rt_array->data.reserve(chromatogram.size());
    intensity_array->data.reserve(chromatogram.size());
    for (MSChromatogram::const_iterator it = chromatogram.begin(); it != chromatogram.end(); ++it)
    {
      rt_array->data.push_back(it->getRT());
      intensity_array->data.push_back(it->getIntensity());
    }
    return cptr;
  }

  void OpenSwathDataAccessHelper::convertToOpenMSChromatogram(const OpenSwath::ChromatogramPtr& cptr, OpenMS::MSChromatogram & chromatogram)
  {
    std::vector<double>::const_iterator rt_it = cptr->getTimeArray()->data.begin();
    std::vector<double>::const_iterator int_it = cptr->getIntensityArray()->data.begin();

    if (!chromatogram.empty()) chromatogram.clear(false);

    ChromatogramPeak peak;
    chromatogram.reserve(cptr->getTimeArray()->data.size());
    for (; rt_it != cptr->getTimeArray()->data.end(); ++rt_it, ++int_it)
    {
      peak.setRT(*rt_it);
      peak.setIntensity(*int_it);
      chromatogram.push_back(peak);
    }
  }

  void OpenSwathDataAccessHelper::convertToOpenMSChromatogramFilter(OpenMS::MSChromatogram & chromatogram,
                                                                    const OpenSwath::ChromatogramPtr& cptr,
                                                                    double rt_min,
                                                                    double rt_max)
  {
    std::vector<double>::const_iterator rt_it = cptr->getTimeArray()->data.begin();
    std::vector<double>::const_iterator int_it = cptr->getIntensityArray()->data.begin();

    ChromatogramPeak peak;
    chromatogram.clear(false);
    chromatogram.reserve(cptr->getTimeArray()->data.size());
    for (; rt_it != cptr->getTimeArray()->data.end(); ++rt_it, ++int_it)
    {
      if (*rt_it < rt_min || *rt_it > rt_max)
      {
        continue;
      }
      peak.setRT(*rt_it);
      peak.setIntensity(*int_it);
      chromatogram.push_back(peak);
    }
  }

  void OpenSwathDataAccessHelper::convertTargetedExp(const OpenMS::TargetedExperiment & transition_exp_, OpenSwath::LightTargetedExperiment & transition_exp)
  {
    //copy proteins
    for (Size i = 0; i < transition_exp_.getProteins().size(); i++)
    {
      OpenSwath::LightProtein p;
      p.id = transition_exp_.getProteins()[i].id;
      transition_exp.proteins.push_back(p);
    }

    //copy peptides and store as compounds
    for (Size i = 0; i < transition_exp_.getPeptides().size(); i++)
    {
      OpenSwath::LightCompound p;
      OpenSwathDataAccessHelper::convertTargetedCompound(transition_exp_.getPeptides()[i], p);
      transition_exp.compounds.push_back(p);
    }

    //copy compounds and store as compounds
    for (Size i = 0; i < transition_exp_.getCompounds().size(); i++)
    {
      OpenSwath::LightCompound c;
      OpenSwathDataAccessHelper::convertTargetedCompound(transition_exp_.getCompounds()[i], c);
      transition_exp.compounds.push_back(c);
    }

    //mapping of transitions
    for (const auto& transition : transition_exp_.getTransitions())
    {
      OpenSwath::LightTransition t;
      t.transition_name = transition.getNativeID();
      t.product_mz = transition.getProductMZ();
      t.precursor_mz = transition.getPrecursorMZ();
      t.library_intensity = transition.getLibraryIntensity();
      t.peptide_ref = transition.getPeptideRef();

      // If compound is a peptide, get the ion mobility information from the compound
      if (!t.peptide_ref.empty())
      {
        OpenSwath::LightCompound p = transition_exp.getPeptideByRef(t.peptide_ref);
        t.precursor_im = p.getDriftTime();
      }
      // try compound ref
      else  // (t.peptide_ref.empty())
      {
        t.peptide_ref = transition.getCompoundRef();
      }

      if (transition.isProductChargeStateSet())
      {
        t.fragment_charge = transition.getProductChargeState();
      }
      t.decoy = false;

      // legacy
#if 1
      const auto& cv_terms = transition.getCVTerms();
      if (cv_terms.find("decoy") != cv_terms.end() && cv_terms.at("decoy")[0].getValue().toString() == "1" )
      {
        t.decoy = true;
      }
      else if (cv_terms.find("MS:1002007") != cv_terms.end())    // target SRM transition
      {
        t.decoy = false;
      }
      else if (cv_terms.find("MS:1002008") != cv_terms.end())    // decoy SRM transition
      {
        t.decoy = true;
      }
      else if (cv_terms.find("MS:1002007") != cv_terms.end() && cv_terms.find("MS:1002008") != cv_terms.end())    // both == illegal
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Transition " + t.transition_name + " cannot be target and decoy at the same time.");
      }
      else
#endif
      if (transition.getDecoyTransitionType() == ReactionMonitoringTransition::UNKNOWN ||
          transition.getDecoyTransitionType() == ReactionMonitoringTransition::TARGET)
      {
        // assume its target
        t.decoy = false;
      }
      else if (transition.getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
      {
        t.decoy = true;
      }

      t.detecting_transition = transition.isDetectingTransition();
      t.identifying_transition = transition.isIdentifyingTransition();
      t.quantifying_transition = transition.isQuantifyingTransition();

      transition_exp.transitions.push_back(t);
    }
  }

  void OpenSwathDataAccessHelper::convertTargetedCompound(const TargetedExperiment::Peptide& pep, OpenSwath::LightCompound & p)
  {
    OpenSwath::LightModification light_mod;

    p.id = pep.id;
    if (pep.hasRetentionTime())
    {
      p.rt = pep.getRetentionTime();
      if (pep.getRetentionTimeUnit() == TargetedExperimentHelper::RetentionTime::RTUnit::MINUTE)
      {
        p.rt = 60 * pep.getRetentionTime();
      }
    }
    p.setDriftTime(pep.getDriftTime());

    if (pep.hasCharge())
    {
      p.charge = pep.getChargeState();
    }

    p.sequence = pep.sequence;
    p.peptide_group_label = pep.getPeptideGroupLabel();

    if (pep.metaValueExists("GeneName"))
    {
      p.gene_name = (std::string)pep.getMetaValue("GeneName");
    }

    // Is it potentially a metabolomics compound
    if (pep.metaValueExists("SumFormula"))
    {
      p.sum_formula = (std::string)pep.getMetaValue("SumFormula");
    }
    if (pep.metaValueExists("CompoundName"))
    {
      p.compound_name = (std::string)pep.getMetaValue("CompoundName");
    }

    p.protein_refs.clear();
    if (!pep.protein_refs.empty())
    {
      p.protein_refs.insert( p.protein_refs.begin(), pep.protein_refs.begin(), pep.protein_refs.end() );
    }

    // Mapping of peptide modifications (don't do this for metabolites...)
    if (p.isPeptide())
    {
      OpenMS::AASequence aa_sequence = TargetedExperimentHelper::getAASequence(pep);

      if (aa_sequence.hasNTerminalModification())
      {
        const ResidueModification& rmod = *(aa_sequence.getNTerminalModification());
        light_mod.location = -1;
        light_mod.unimod_id = rmod.getUniModRecordId();
        p.modifications.push_back(light_mod);
      }
      if (aa_sequence.hasCTerminalModification())
      {
        const ResidueModification& rmod = *(aa_sequence.getCTerminalModification());
        light_mod.location = boost::numeric_cast<int>(aa_sequence.size());
        light_mod.unimod_id = rmod.getUniModRecordId();
        p.modifications.push_back(light_mod);
      }
      for (Size i = 0; i != aa_sequence.size(); i++)
      {
        if (aa_sequence[i].isModified())
        {
          // search the residue in the modification database (if the sequence is valid, we should find it)
          const ResidueModification& rmod = *(aa_sequence.getResidue(i).getModification());
          light_mod.location = boost::numeric_cast<int>(i);
          light_mod.unimod_id = rmod.getUniModRecordId();
          p.modifications.push_back(light_mod);
        }
      }

    }
  }

  void OpenSwathDataAccessHelper::convertTargetedCompound(const TargetedExperiment::Compound& compound, OpenSwath::LightCompound & comp)
  {
    comp.id = compound.id;
    if (compound.hasRetentionTime())
    {
      comp.rt = compound.getRetentionTime();
      if (compound.getRetentionTimeUnit() == TargetedExperimentHelper::RetentionTime::RTUnit::MINUTE)
      {
        comp.rt = 60 * compound.getRetentionTime();
      }
    }
    comp.setDriftTime(compound.getDriftTime());

    if (compound.hasCharge())
    {
      comp.charge = compound.getChargeState();
    }

    comp.sum_formula = (std::string)compound.molecular_formula;
    if (compound.metaValueExists("CompoundName"))
    {
      comp.compound_name = (std::string)compound.getMetaValue("CompoundName");
    }
  }

  void OpenSwathDataAccessHelper::convertPeptideToAASequence(const OpenSwath::LightCompound & peptide, AASequence & aa_sequence)
  {
    OPENMS_PRECONDITION(peptide.isPeptide(), "Function needs peptide, not metabolite")

    aa_sequence = AASequence::fromString(peptide.sequence);
    for (const auto & it : peptide.modifications)
    {
      if (it.unimod_id != -1)
      {
        TargetedExperimentHelper::setModification(it.location,
                                                  int(peptide.sequence.size()),
                                                  "UniMod:" + String(it.unimod_id), aa_sequence);
      }
    }
  }

}
