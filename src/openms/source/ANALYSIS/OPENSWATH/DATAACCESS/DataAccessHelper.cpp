// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

namespace OpenMS
{
  void OpenSwathDataAccessHelper::convertToOpenMSSpectrum(const OpenSwath::SpectrumPtr sptr, OpenMS::MSSpectrum & spectrum)
  {
    spectrum.reserve(sptr->getMZArray()->data.size());

    std::vector<double>::const_iterator mz_it = sptr->getMZArray()->data.begin();
    std::vector<double>::const_iterator int_it = sptr->getIntensityArray()->data.begin();
    for (; mz_it != sptr->getMZArray()->data.end(); ++mz_it, ++int_it)
    {
      Peak1D p;
      p.setMZ(*mz_it);
      p.setIntensity(*int_it);
      spectrum.push_back(p);
    }
  }

  OpenSwath::SpectrumPtr OpenSwathDataAccessHelper::convertToSpectrumPtr(const OpenMS::MSSpectrum & spectrum)
  {
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    for (MSSpectrum::const_iterator it = spectrum.begin(); it != spectrum.end(); ++it)
    {
      mz_array->data.push_back(it->getMZ());
      intensity_array->data.push_back(it->getIntensity());
    }

    OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
    sptr->setMZArray(mz_array);
    sptr->setIntensityArray(intensity_array);
    return sptr;
  }

  void OpenSwathDataAccessHelper::convertToOpenMSChromatogram(const OpenSwath::ChromatogramPtr cptr, OpenMS::MSChromatogram & chromatogram)
  {
    chromatogram.reserve(cptr->getTimeArray()->data.size());

    std::vector<double>::const_iterator rt_it = cptr->getTimeArray()->data.begin();
    std::vector<double>::const_iterator int_it = cptr->getIntensityArray()->data.begin();
    for (; rt_it != cptr->getTimeArray()->data.end(); ++rt_it, ++int_it)
    {
      ChromatogramPeak peak;
      peak.setRT(*rt_it);
      peak.setIntensity(*int_it);
      chromatogram.push_back(peak);
    }
  }

  void OpenSwathDataAccessHelper::convertToOpenMSChromatogramFilter(OpenMS::MSChromatogram & chromatogram, const OpenSwath::ChromatogramPtr cptr,
                                                                    double rt_min, double rt_max)
  {
    chromatogram.reserve(cptr->getTimeArray()->data.size());

    std::vector<double>::const_iterator rt_it = cptr->getTimeArray()->data.begin();
    std::vector<double>::const_iterator int_it = cptr->getIntensityArray()->data.begin();
    for (; rt_it != cptr->getTimeArray()->data.end(); ++rt_it, ++int_it)
    {
      if (*rt_it < rt_min || *rt_it > rt_max)
      {
        continue;
      }
      ChromatogramPeak peak;
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
    for (Size i = 0; i < transition_exp_.getTransitions().size(); i++)
    {
      OpenSwath::LightTransition t;
      t.transition_name = transition_exp_.getTransitions()[i].getNativeID();
      t.product_mz = transition_exp_.getTransitions()[i].getProductMZ();
      t.precursor_mz = transition_exp_.getTransitions()[i].getPrecursorMZ();
      t.library_intensity = transition_exp_.getTransitions()[i].getLibraryIntensity();
      t.peptide_ref = transition_exp_.getTransitions()[i].getPeptideRef();
      // try compound ref
      if (t.peptide_ref.empty())
      {
        t.peptide_ref = transition_exp_.getTransitions()[i].getCompoundRef();
      }
      if (transition_exp_.getTransitions()[i].isProductChargeStateSet())
      {
        t.fragment_charge = transition_exp_.getTransitions()[i].getProductChargeState();
      }
      t.decoy = false;

      // legacy
#if 1
      if (transition_exp_.getTransitions()[i].getCVTerms().has("decoy") &&
          transition_exp_.getTransitions()[i].getCVTerms()["decoy"][0].getValue().toString() == "1" )
      {
        t.decoy = true;
      }
      else if (transition_exp_.getTransitions()[i].getCVTerms().has("MS:1002007"))    // target SRM transition
      {
        t.decoy = false;
      }
      else if (transition_exp_.getTransitions()[i].getCVTerms().has("MS:1002008"))    // decoy SRM transition
      {
        t.decoy = true;
      }
      else if (transition_exp_.getTransitions()[i].getCVTerms().has("MS:1002007") &&
          transition_exp_.getTransitions()[i].getCVTerms().has("MS:1002008"))    // both == illegal
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Transition " + t.transition_name + " cannot be target and decoy at the same time.");
      }
      else
#endif
      if (transition_exp_.getTransitions()[i].getDecoyTransitionType() == ReactionMonitoringTransition::UNKNOWN ||
          transition_exp_.getTransitions()[i].getDecoyTransitionType() == ReactionMonitoringTransition::TARGET)
      {
        // assume its target
        t.decoy = false;
      }
      else if (transition_exp_.getTransitions()[i].getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
      {
        t.decoy = true;
      }

      t.detecting_transition = transition_exp_.getTransitions()[i].isDetectingTransition();
      t.identifying_transition = transition_exp_.getTransitions()[i].isIdentifyingTransition();
      t.quantifying_transition = transition_exp_.getTransitions()[i].isQuantifyingTransition();

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
    }

    if (pep.hasCharge())
    {
      p.charge = pep.getChargeState();
    }

    p.sequence = pep.sequence;
    p.peptide_group_label = pep.getPeptideGroupLabel();

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
    }

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
    for (std::vector<OpenSwath::LightModification>::const_iterator it = peptide.modifications.begin();
        it != peptide.modifications.end(); ++it)
    {
      TargetedExperimentHelper::setModification(it->location, 
                                                boost::numeric_cast<int>(peptide.sequence.size()), 
                                                "UniMod:" + String(it->unimod_id), aa_sequence);
    }
  }


}
