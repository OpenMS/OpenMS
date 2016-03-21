// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

namespace OpenMS
{
  void OpenSwathDataAccessHelper::convertToOpenMSSpectrum(const OpenSwath::SpectrumPtr sptr, OpenMS::MSSpectrum<> & spectrum)
  {
    // recreate a spectrum from the data arrays!
    OpenSwath::BinaryDataArrayPtr mz_arr = sptr->getMZArray();
    OpenSwath::BinaryDataArrayPtr int_arr = sptr->getIntensityArray();
    spectrum.reserve(mz_arr->data.size());
    for (Size i = 0; i < mz_arr->data.size(); i++)
    {
      Peak1D p;
      p.setMZ(mz_arr->data[i]);
      p.setIntensity(int_arr->data[i]);
      spectrum.push_back(p);
    }
  }

  OpenSwath::SpectrumPtr OpenSwathDataAccessHelper::convertToSpectrumPtr(const OpenMS::MSSpectrum<> & spectrum)
  {
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    for (MSSpectrum<>::const_iterator it = spectrum.begin(); it != spectrum.end(); ++it)
    {
      mz_array->data.push_back(it->getMZ());
      intensity_array->data.push_back(it->getIntensity());
    }

    OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
    sptr->setMZArray(mz_array);
    sptr->setIntensityArray(intensity_array);
    return sptr;
  }

  void OpenSwathDataAccessHelper::convertToOpenMSChromatogram(OpenMS::MSChromatogram<> & chromatogram, const OpenSwath::ChromatogramPtr cptr)
  {
    OpenSwath::BinaryDataArrayPtr rt_arr = cptr->getTimeArray();
    OpenSwath::BinaryDataArrayPtr int_arr = cptr->getIntensityArray();
    chromatogram.reserve(rt_arr->data.size());
    for (Size i = 0; i < rt_arr->data.size(); i++)
    {
      ChromatogramPeak p;
      p.setRT(rt_arr->data[i]);
      p.setIntensity(int_arr->data[i]);
      chromatogram.push_back(p);
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

    //copy peptides
    for (Size i = 0; i < transition_exp_.getPeptides().size(); i++)
    {
      OpenSwath::LightPeptide p;
      OpenSwathDataAccessHelper::convertTargetedPeptide(transition_exp_.getPeptides()[i], p);
      transition_exp.peptides.push_back(p);
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
      t.product_charge = transition_exp_.getTransitions()[i].getProduct().getChargeState();
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
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
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

      if (transition_exp_.getTransitions()[i].metaValueExists("detecting_transition"))
      {
        if (!transition_exp_.getTransitions()[i].getMetaValue("detecting_transition").toBool())
        {
          t.detecting_transition = false;
        }
        else if (transition_exp_.getTransitions()[i].getMetaValue("detecting_transition").toBool())
        {
          t.detecting_transition = true;
        }
      }
      else
      {
        t.detecting_transition = true;
      }
      if (transition_exp_.getTransitions()[i].metaValueExists("identifying_transition"))
      {
        if (!transition_exp_.getTransitions()[i].getMetaValue("identifying_transition").toBool())
        {
          t.identifying_transition = false;
        }
        else if (transition_exp_.getTransitions()[i].getMetaValue("identifying_transition").toBool())
        {
          t.identifying_transition = true;
        }
      }
      else
      {
        t.identifying_transition = false;
      }
      if (transition_exp_.getTransitions()[i].metaValueExists("quantifying_transition"))
      {
        if (!transition_exp_.getTransitions()[i].getMetaValue("quantifying_transition").toBool())
        {
          t.quantifying_transition = false;
        }
        else if (transition_exp_.getTransitions()[i].getMetaValue("quantifying_transition").toBool())
        {
          t.quantifying_transition = true;
        }
      }
      else
      {
        t.quantifying_transition = true;
      }

      transition_exp.transitions.push_back(t);
    }
  }

  void OpenSwathDataAccessHelper::convertTargetedPeptide(const TargetedExperiment::Peptide& pep, OpenSwath::LightPeptide & p)
  {
    OpenSwath::LightModification m;
    OpenMS::ModificationsDB* mod_db = OpenMS::ModificationsDB::getInstance();

    p.id = pep.id;
    if (!pep.rts.empty())
    {
      p.rt = pep.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
    }
    p.charge = pep.getChargeState();
    p.sequence = pep.sequence;
    p.peptide_group_label = pep.getPeptideGroupLabel();

    p.protein_refs.clear();
    if (!pep.protein_refs.empty())
    {
      p.protein_refs.insert( p.protein_refs.begin(), pep.protein_refs.begin(), pep.protein_refs.end() ); 
    }

    // Mapping of peptide modifications
    {
      OpenMS::AASequence aa_sequence = TargetedExperimentHelper::getAASequence(pep);
      if ( !aa_sequence.getNTerminalModification().empty())
      {
          ResidueModification rmod = mod_db->getTerminalModification(aa_sequence.getNTerminalModification(), ResidueModification::N_TERM);
          m.location = -1;
          m.unimod_id = rmod.getUniModAccession();
          p.modifications.push_back(m);
      }
      if ( !aa_sequence.getCTerminalModification().empty())
      {
          ResidueModification rmod = mod_db->getTerminalModification(aa_sequence.getCTerminalModification(), ResidueModification::C_TERM);
          m.location = boost::numeric_cast<int>(aa_sequence.size());
          m.unimod_id = rmod.getUniModAccession();
          p.modifications.push_back(m);
      }
      for (Size i = 0; i != aa_sequence.size(); i++)
      {
        if (aa_sequence[i].isModified())
        {
          // search the residue in the modification database (if the sequence is valid, we should find it)
          ResidueModification rmod = mod_db->getModification(aa_sequence.getResidue(i).getOneLetterCode(),
                                                             aa_sequence.getResidue(i).getModification(), ResidueModification::ANYWHERE);
          m.location = boost::numeric_cast<int>(i);
          m.unimod_id = rmod.getUniModAccession();
          p.modifications.push_back(m);

        }
      }

    }
    // transition_exp.peptides.push_back(p);
  }

  void OpenSwathDataAccessHelper::convertPeptideToAASequence(const OpenSwath::LightPeptide & peptide, AASequence & aa_sequence)
  {
      aa_sequence = AASequence::fromString(peptide.sequence);
      for (std::vector<OpenSwath::LightModification>::const_iterator it = peptide.modifications.begin(); it != peptide.modifications.end(); ++it)
      {
        TargetedExperimentHelper::setModification(it->location, boost::numeric_cast<int>(peptide.sequence.size()), it->unimod_id, aa_sequence);
      }
  }


}
