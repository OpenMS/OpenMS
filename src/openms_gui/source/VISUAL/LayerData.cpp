// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <iostream>

using namespace std;

namespace OpenMS
{
  const std::string LayerData::NamesOfLabelType[] = {"None", "Index", "Label meta data", "Peptide identification", "All peptide identifications"};

  std::ostream & operator<<(std::ostream & os, const LayerData & rhs)
  {
    os << "--LayerData BEGIN--" << std::endl;
    os << "name: " << rhs.name << std::endl;
    os << "visible: " << rhs.visible << std::endl;
    os << "number of peaks: " << rhs.getPeakData()->getSize() << std::endl;
    os << "--LayerData END--" << std::endl;
    return os;
  }

  const LayerData::ConstExperimentSharedPtrType LayerData::getPeakData() const
  {
    return boost::static_pointer_cast<const ExperimentType>(peaks);
  }

  void LayerData::updateRanges()
  {
    peaks->updateRanges();
    features->updateRanges();
    consensus->updateRanges();
    // on_disc_peaks->updateRanges(); // note: this is not going to work since its on disk! We currently don't have a good way to access these ranges
    chromatograms->updateRanges();
    cached_spectrum_.updateRanges();
  }

  void LayerData::updateCache_()
  {
    if (peaks->getNrSpectra() > current_spectrum_ && (*peaks)[current_spectrum_].size() > 0)
    {
      cached_spectrum_ = (*peaks)[current_spectrum_];
    }
    else if (on_disc_peaks->getNrSpectra() > current_spectrum_)
    {
      cached_spectrum_ = on_disc_peaks->getSpectrum(current_spectrum_);
    }
  }

  const LayerData::ExperimentType::SpectrumType & LayerData::getCurrentSpectrum() const
  {
    return cached_spectrum_;
  }

  void LayerData::synchronizePeakAnnotations()
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->empty() || type != LayerData::DT_PEAK) { return; }

    // get mutable access to the spectrum
    MSSpectrum & spectrum = getPeakDataMuteable()->getSpectrum(current_spectrum_);

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
        { // no hits? add empty hit
          PeptideHit hit;
          updatePeptideHitAnnotations_(hit);
          hits.push_back(hit);
        }
      }
      else // PeptideIdentifications are empty, create new PepIDs and PeptideHits to store the PeakAnnotations
      {
        // copy user annotations to fragment annotation vector
        Annotations1DContainer & las = getAnnotations(current_spectrum_);

        // no annotations so we don't need to synchronize
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
        vector<ProteinIdentification>& prot_ids = getPeakDataMuteable()->getProteinIdentifications();
        if (prot_ids.empty() || prot_ids.back().getIdentifier() != String("Unknown"))
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
    Annotations1DContainer & las = getCurrentAnnotations();

    // initialize with an empty vector
    vector<PeptideHit::PeakAnnotation> fas;

    // do not change PeptideHit annotations, if there are no annotations on the spectrum
    bool annotations_changed(false);

    // regular expression for a charge at the end of the annotation
    QRegExp reg_exp("([\\+|\\-]\\d+)$");

    // for each annotation item on the canvas
    for (auto& a : las)
    {
      // only store peak annotations (skip general lables and distance annotations)
      Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(a);
      if (pa == nullptr) { continue; }

      // add new fragment annotation
      QString peak_anno = pa->getText().trimmed();

      // check for newlines in the label and only continue with the first line for charge determination
      QStringList lines = peak_anno.split(QRegExp("[\r\n]"), QString::SkipEmptyParts);
      if (lines.size() > 1)
      {
        peak_anno = lines[0];
      }

      // read charge and text from annotation item string
      // we support two notations for the charge suffix: '2+' or '++'
      // cut and convert the trailing + or - to a proper charge
      int match_pos = reg_exp.indexIn(peak_anno);
      int tmp_charge(0);
      if (match_pos >= 0)
      {
        tmp_charge = reg_exp.cap(1).toInt();
        peak_anno = peak_anno.left(match_pos);
      }
      else
      {
        // count number of + and - in suffix (e.g., to support "++" as charge 2 anotation)
        int plus(0), minus(0);

        for (int p = (int)peak_anno.size() - 1; p >= 0; --p)
        {
          if (peak_anno[p] == '+')
          {
            ++plus;
            continue;
          }
          else if (peak_anno[p] == '-')
          {
            ++minus;
            continue;
          }
          else // not '+' or '-'?
          {
            if (plus > 0 && minus == 0) // found pluses?
            {
              tmp_charge = plus;
              peak_anno = peak_anno.left(peak_anno.size() - plus);
              break;
            }
            else if (minus > 0 && plus == 0)  // found minuses?
            {
              tmp_charge = -minus;
              peak_anno = peak_anno.left(peak_anno.size() - minus);
              break;
            }
            break;
          }
        }
      }

      PeptideHit::PeakAnnotation fa;
      fa.charge = tmp_charge;
      fa.mz = pa->getPeakPosition()[0];
      fa.intensity = pa->getPeakPosition()[1];
      if (lines.size() > 1)
      {
        peak_anno.append("\n").append(lines[1]);
      }
      fa.annotation = peak_anno;

      fas.push_back(fa);
      annotations_changed = true;
    }

    if (annotations_changed)
    {
      hit.setPeakAnnotations(fas);
    }
  }

  void LayerData::removePeakAnnotationsFromPeptideHit(const std::vector<Annotation1DItem*>& selected_annotations)
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->empty() || type != LayerData::DT_PEAK) { return; }

    // no ID selected
    if (peptide_id_index == -1 || peptide_hit_index == -1) { return; }

    // get mutable access to the spectrum
    MSSpectrum & spectrum = getPeakDataMuteable()->getSpectrum(current_spectrum_);
    int ms_level = spectrum.getMSLevel();

    // wrong MS level
    if (ms_level < 2) { return; }

    // extract PeptideIdentification and PeptideHit if possible.
    // that this function returns prematurely is unlikely,
    // since we are deleting existing annotations,
    // that have to be somewhere, but better make sure
    vector<PeptideIdentification>& pep_ids = spectrum.getPeptideIdentifications();
    if (pep_ids.empty()) { return; }
    vector<PeptideHit>& hits = pep_ids[peptide_id_index].getHits();
    if (hits.empty()) { return; }
    PeptideHit& hit = hits[peptide_hit_index];
    vector<PeptideHit::PeakAnnotation> fas = hit.getPeakAnnotations();
    if (fas.empty()) { return; }

    // all requirements fulfilled, PH in hit and annotations in selected_annotations
    vector<PeptideHit::PeakAnnotation> to_remove;
    bool annotations_changed(false);

    // collect annotations, that have to be removed
    for (auto const& tmp_a : fas)
    {
      for (auto const& it : selected_annotations)
      {
        Annotation1DPeakItem* pa = dynamic_cast<Annotation1DPeakItem*>(it);
        // only search for peak annotations
        if (pa == nullptr) { continue; }
        if (fabs(tmp_a.mz - pa->getPeakPosition()[0]) < 1e-6)
        {
          if (String(pa->getText()).hasPrefix(tmp_a.annotation))
          {
            to_remove.push_back(tmp_a);
            annotations_changed = true;
          }
        }
      }
    }
    // remove the collected annotations from the PeptideHit annotations
    for (auto const& tmp_a : to_remove)
    {
      fas.erase(std::remove(fas.begin(), fas.end(), tmp_a), fas.end());
    }
    if (annotations_changed) { hit.setPeakAnnotations(fas); }
  }

} //Namespace
