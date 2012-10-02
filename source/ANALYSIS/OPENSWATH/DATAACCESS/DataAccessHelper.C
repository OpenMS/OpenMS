// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

namespace OpenMS
{
  void OpenSwathDataAccessHelper::convertToOpenMSSpectrum(OpenMS::MSSpectrum<> & spectrum, const OpenSwath::SpectrumPtr sptr)
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

  void OpenSwathDataAccessHelper::convertToOpenMSChromatogram(OpenMS::MSChromatogram<> & chromatogram,
                                                              const OpenSwath::ChromatogramPtr cptr)
  {
    // recreate a spectrum from the data arrays!
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

  void OpenSwathDataAccessHelper::convertTargetedExp(OpenMS::TargetedExperiment & transition_exp_, OpenSwath::LightTargetedExperiment & transition_exp)
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
      p.id = transition_exp_.getPeptides()[i].id;
      p.rt = transition_exp_.getPeptides()[i].rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
      p.charge = transition_exp_.getPeptides()[i].getChargeState();
      p.sequence = transition_exp_.getPeptides()[i].sequence;
      p.protein_ref = transition_exp_.getPeptides()[i].protein_refs[0];

      // TODO (hroest/georger) merge with code in MRMDecoyGenerator
      // The workaround expects a TraML with UniMod CVTerms for each peptide. The problem in ModificationsDB has been described in TRAC #458.
      // Mapping of peptide modifications
      const OpenMS::TargetedExperiment::Peptide * pep = &transition_exp_.getPeptides()[i];
      for (std::vector<OpenMS::TargetedExperiment::Peptide::Modification>::const_iterator it =
             pep->mods.begin(); it != pep->mods.end(); ++it)
      {
        const Map<OpenMS::String, std::vector<CVTerm> > cv_terms = it->getCVTerms();
        for (Map<OpenMS::String, std::vector<CVTerm> >::const_iterator li = cv_terms.begin();
             li != cv_terms.end(); ++li)
        {
          std::vector<CVTerm> mods = (*li).second;
          for (std::vector<CVTerm>::iterator mo = mods.begin(); mo != mods.end(); ++mo)
          {
            OpenSwath::LightModification m;
            m.location = it->location;
            m.unimod_id = mo->getAccession().substr(7);
            p.modifications.push_back(m);
          }
        }
      }
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
      t.charge = transition_exp_.getTransitions()[i].getProduct().getChargeState();
      transition_exp.transitions.push_back(t);
    }
  }


}
