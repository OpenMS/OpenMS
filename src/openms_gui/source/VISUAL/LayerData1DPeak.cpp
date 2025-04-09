// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerData1DPeak.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/Painter1DBase.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

#include <QMenu>
using namespace std;

namespace OpenMS
{

  std::unique_ptr<LayerStoreData> LayerData1DPeak::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = std::make_unique<LayerStoreDataPeakMapVisible>();
    ret->storeVisibleSpectrum(getCurrentSpectrum(), visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerStoreData> LayerData1DPeak::storeFullData() const
  {
    return LayerDataPeak::storeFullData(); // just forward
  }

  QMenu* LayerData1DPeak::getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint)
  {
    auto* context_menu = new QMenu("Peak1D", nullptr);
    context_menu->addAction("Edit", [annot_item, &need_repaint, this]() { // this capture is tricky! Copy 'annot_item' since its a local variable and will be out of scope when the menu is evaluated!
      annot_item->editText();
      synchronizePeakAnnotations();
      need_repaint = true;
    });
    context_menu->addAction("Delete", [annot_item, &need_repaint, this]() { // this capture is tricky! Copy 'annot_item' since its a local variable and will be out of scope when the menu is evaluated!
      vector<Annotation1DItem*> as;
      as.push_back(annot_item);
      removePeakAnnotationsFromPeptideHit(as);
      getCurrentAnnotations().removeSelectedItems();
      need_repaint = true;
    });

    return context_menu;
  }

  PeakIndex LayerData1DPeak::findClosestDataPoint(const RangeAllType& area) const
  {
    Peak1D peak_lt(area.getMinMZ(), area.getMinIntensity()), peak_rb(area.getMaxMZ(), area.getMaxIntensity());
    // reference to the current data
    const auto& spectrum = getCurrentSpectrum();
    const Size spectrum_index = getCurrentIndex();

    // get iterator on first peak with lower position than interval_start
    auto left_it = lower_bound(spectrum.begin(), spectrum.end(), peak_lt, PeakType::PositionLess());

    // get iterator on first peak with higher position than interval_end
    auto right_it = lower_bound(left_it, spectrum.end(), peak_rb, PeakType::PositionLess());

    if (left_it == right_it) // both are equal => no peak falls into this interval
    {
      return PeakIndex();
    }

    if (left_it == right_it - 1)
    {
      return PeakIndex(spectrum_index, left_it - spectrum.begin());
    }

    auto nearest_it = left_it;
    const auto center_intensity = (peak_lt.getIntensity() + peak_rb.getIntensity()) * 0.5;
    for (auto it = left_it; it != right_it; ++it)
    {
      if (abs(center_intensity - it->getIntensity()) < abs(center_intensity - nearest_it->getIntensity()))
      {
        nearest_it = it;
      }
    }
    return PeakIndex(spectrum_index, nearest_it - spectrum.begin());
  }

  std::unique_ptr<Painter1DBase> LayerData1DPeak::getPainter1D() const
  {
    return make_unique<Painter1DPeak>(this);
  }

  Annotation1DItem* LayerData1DPeak::addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color)
  {
    auto peak = getCurrentSpectrum()[peak_index.peak];
    auto* item = new Annotation1DPeakItem<decltype(peak)>(peak, text, color);
    item->setSelected(false);
    getCurrentAnnotations().push_front(item);
    return item;
  }

  void LayerData1DPeak::synchronizePeakAnnotations()
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->empty() || type != LayerDataBase::DT_PEAK)
    {
      return;
    }

    // get mutable access to the spectrum
    MSSpectrum& spectrum = getPeakDataMuteable()->getSpectrum(current_idx_);

    int ms_level = spectrum.getMSLevel();

    if (ms_level != 2)
      return;

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
      const Annotations1DContainer& las = getAnnotations(current_idx_);

      // no annotations so we don't need to synchronize
      bool has_peak_annotation(false);
      for (auto& a : las)
      {
        // only store peak annotations
        auto pa = dynamic_cast<Annotation1DPeakItem<Peak1D>*>(a);
        if (pa != nullptr)
        {
          has_peak_annotation = true;
          break;
        }
      }
      if (has_peak_annotation == false)
      {
        return;
      }
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

  void LayerData1DPeak::removePeakAnnotationsFromPeptideHit(const std::vector<Annotation1DItem*>& selected_annotations)
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->empty() || type != LayerDataBase::DT_PEAK)
    {
      return;
    }

    // no ID selected
    if (peptide_id_index == -1 || peptide_hit_index == -1)
    {
      return;
    }

    // get mutable access to the spectrum
    MSSpectrum& spectrum = getPeakDataMuteable()->getSpectrum(current_idx_);
    int ms_level = spectrum.getMSLevel();

    // wrong MS level
    if (ms_level < 2)
    {
      return;
    }

    // extract PeptideIdentification and PeptideHit if possible.
    // that this function returns prematurely is unlikely,
    // since we are deleting existing annotations,
    // that have to be somewhere, but better make sure
    vector<PeptideIdentification>& pep_ids = spectrum.getPeptideIdentifications();
    if (pep_ids.empty())
    {
      return;
    }
    vector<PeptideHit>& hits = pep_ids[peptide_id_index].getHits();
    if (hits.empty())
    {
      return;
    }
    PeptideHit& hit = hits[peptide_hit_index];
    vector<PeptideHit::PeakAnnotation> fas = hit.getPeakAnnotations();
    if (fas.empty())
    {
      return;
    }

    // all requirements fulfilled, PH in hit and annotations in selected_annotations
    vector<PeptideHit::PeakAnnotation> to_remove;
    // collect annotations, that have to be removed
    for (auto const& tmp_a : fas)
    {
      for (auto const& it : selected_annotations)
      {
        using ItemType = Peak1D;
        auto pa = dynamic_cast<Annotation1DPeakItem<ItemType>*>(it);
        // only search for peak annotations
        if (pa == nullptr)
        {
          continue;
        }

        if (fabs(tmp_a.mz - pa->getPeakPosition().getMZ()) < 1e-6)
        {
          if (String(pa->getText()).hasPrefix(tmp_a.annotation))
          {
            to_remove.push_back(tmp_a);
          }
        }
      }
    }
    // remove the collected annotations from the PeptideHit annotations
    for (auto const& tmp_a : to_remove)
    {
      fas.erase(std::remove(fas.begin(), fas.end(), tmp_a), fas.end());
    }
    if (!to_remove.empty())
    {
      hit.setPeakAnnotations(fas);
    }
  }

  void LayerData1DPeak::updatePeptideHitAnnotations_(PeptideHit& hit)
  {
    // copy user annotations to fragment annotation vector
    const Annotations1DContainer& las = getCurrentAnnotations();

    // initialize with an empty vector
    vector<PeptideHit::PeakAnnotation> fas;

    // do not change PeptideHit annotations, if there are no annotations on the spectrum
    bool annotations_changed(false);

    // for each annotation item on the canvas
    for (auto& a : las)
    {
      // only store peak annotations (skip general labels and distance annotations)
      auto pa = dynamic_cast<Annotation1DPeakItem<Peak1D>*>(a);
      if (pa == nullptr)
      {
        continue;
      }
      fas.push_back(pa->toPeakAnnotation());
      annotations_changed = true;
    }

    if (annotations_changed)
    {
      hit.setPeakAnnotations(fas);
    }
  }


}// namespace OpenMS
