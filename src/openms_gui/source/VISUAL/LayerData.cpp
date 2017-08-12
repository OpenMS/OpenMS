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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>

using namespace std;

namespace OpenMS
{
  const std::string LayerData::NamesOfLabelType[] = {"None", "Index", "Label meta data", "Peptide identification", "All peptide identifications"};

  const LayerData::ExperimentType::SpectrumType & LayerData::getCurrentSpectrum() const
  {
    return (*peaks)[current_spectrum_];
  }

  std::ostream & operator<<(std::ostream & os, const LayerData & rhs)
  {
    os << "--LayerData BEGIN--" << std::endl;
    os << "name: " << rhs.name << std::endl;
    os << "visible: " << rhs.visible << std::endl;
    os << "number of peaks: " << rhs.getPeakData()->getSize() << std::endl;
    os << "--LayerData END--" << std::endl;
    return os;
  }

  void LayerData::synchronizePeakAnnotations()
  {
    // LayerData & current_layer = widget_1D->canvas()->getCurrentLayer();
    int spectrum_index = getCurrentSpectrumIndex();

    // Return if no valid peak layer attached
    if (getPeakData()->size() == 0 || type != LayerData::DT_PEAK) { return; }

    MSSpectrum<> & spectrum = (*getPeakData())[spectrum_index];
    int ms_level = spectrum.getMSLevel();

    if (ms_level == 2)
    {
      // store user fragment annotations
      vector<PeptideIdentification>& pep_ids = spectrum.getPeptideIdentifications();

      // no ID selected
      if (peptide_id_index == -1 || peptide_hit_index == -1)
      {
        return;
      }

      if (!pep_ids.empty())
      {
        vector<PeptideHit>& hits = pep_ids[peptide_id_index].getHits();

        if (!hits.empty())
        {
          PeptideHit& hit = hits[peptide_hit_index];
          updatePeptideHitAnnotations_(hit);
        }
        else
        {
          PeptideHit hit;
          updatePeptideHitAnnotations_(hit);
          hits.push_back(hit);
        }
      }
      else // PeptideIdentifications are empty, create new PepIDs and PeptideHits to store the PeakAnnotations
      {
        // copy user annotations to fragment annotation vector
        Annotations1DContainer & las = getAnnotations(current_spectrum_);

        // no annoations so we don't need to synchronize
        bool has_peak_annotation(false);
        for (auto& a : las)
        {
          // only store peak annotations
          Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(a);
          if (pa != nullptr) { has_peak_annotation = true; break; }
        }
        if (has_peak_annotation == false) return;

        PeptideIdentification pep_id;
        pep_id.setIdentifier("Unknown");

        // create a dummy ProteinIdentification for all ID-less PeakAnnotations
        vector<ProteinIdentification>& prot_ids = getPeakData()->getProteinIdentifications();
        if (prot_ids.back().getIdentifier() != String("Unknown"))
        {
          ProteinIdentification prot_id;
          prot_id.setIdentifier("Unknown");
          prot_ids.push_back(prot_id);
        }

        PeptideHit hit;
        if (spectrum.getPrecursors().empty() == false)
        {
          pep_id.setMZ(spectrum.getPrecursors()[0].getMZ());
          hit.setCharge(spectrum.getPrecursors()[0].getCharge());
        }
        pep_id.setRT(spectrum.getRT());

        updatePeptideHitAnnotations_(hit);
        std::vector<PeptideHit> hits;
        hits.push_back(hit);
        pep_id.setHits(hits);
        pep_ids.push_back(pep_id);
      }
    }
  }

  void LayerData::updatePeptideHitAnnotations_(PeptideHit& hit)
  {
    // copy user annotations to fragment annotation vector
    Annotations1DContainer & las = getAnnotations(current_spectrum_);

    vector<PeptideHit::PeakAnnotation> fas = hit.getPeakAnnotations();

    bool annotations_changed(false);

    // regular expression for a charge at the end of the annotation
    QRegExp reg_exp("([\\+|\\-]\\d+)$");

    // for each annotation item on the canvas
    for (auto& a : las)
    {
      // only store peak annotations
      Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(a);
      if (pa == nullptr) { continue; }

      int tmp_charge(0);

      // if already annotated we want to keep mz, intensity, and charge information
      bool already_annotated(false);
      for (auto& tmp_a : fas)
      {
        if (fabs(tmp_a.mz - pa->getPeakPosition()[0]) < 1e-6)
        {
          QString peak_anno = pa->getText();
          int match_pos = reg_exp.indexIn(peak_anno);
          // if a charge was found at the annotation, remove it from the annotation string an fill the charge attribute
          if (match_pos >= 0)
          {
            tmp_charge = reg_exp.cap(1).toInt();
            peak_anno = peak_anno.left(match_pos);
          }

          if (tmp_a.annotation == String(peak_anno))
          {
            already_annotated = true;
            break;
          }
          else // peak annotated but different text (e.g., changed by user)
          {
            tmp_a.annotation = peak_anno;
            // if the new annotation has a charge, use the new one
            if (tmp_charge != 0)
            {
              tmp_a.charge = tmp_charge;
            }
            annotations_changed = true;
            already_annotated = true;
            break;
          }
        }
      }

      // add new fragment annotation if peak not yet annotated
      if (!already_annotated)
      {
        QString peak_anno = pa->getText();
        int match_pos = reg_exp.indexIn(peak_anno);
        if (match_pos >= 0)
        {
          tmp_charge = reg_exp.cap(1).toInt();
          peak_anno = peak_anno.left(match_pos);
        }

        PeptideHit::PeakAnnotation fa;
        fa.charge = tmp_charge;
        fa.mz = pa->getPeakPosition()[0];
        fa.intensity = pa->getPeakPosition()[1];
        fa.annotation = peak_anno;
        fas.push_back(fa);
        annotations_changed = true;
      }
    }
    if (annotations_changed) { hit.setPeakAnnotations(fas); }
  }

} //Namespace
