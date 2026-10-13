// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerData1DPeak.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>
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

  std::optional<IdentificationData::MatchReference> LayerData1DPeak::getSelectedMatch() const
  {
    if (getPeakData() == nullptr || getPeakData()->getMSExperiment().empty() || peptide_id_index < 0 || peptide_hit_index < 0)
    {
      return std::nullopt;
    }
    const auto identifications = getPeakData()->getSpectrumIdentifications(getCurrentIndex());
    if (Size(peptide_id_index) >= identifications.size())
    {
      return std::nullopt;
    }
    const auto& identification = identifications[peptide_id_index];
    if (Size(peptide_hit_index) >= identification.matches.size())
    {
      return std::nullopt;
    }
    return IdentificationData::MatchReference {identification.run->getUuid(), identification.matches[peptide_hit_index]->getId()};
  }

  namespace
  {
    /// Replace the peak annotations of the match @p reference of @p data
    void setPeakAnnotations(IdentificationData& data, const IdentificationData::MatchReference& reference,
                            std::vector<PeptideHit::PeakAnnotation> annotations)
    {
      auto* run = data.findRunByUuid(reference.run_uuid);
      IdentificationData::MatchData match = run->getMatch(reference.match).getData();
      match.peak_annotations = std::move(annotations);
      run->replaceMatch(reference.match, match);
    }
  }

  void LayerData1DPeak::synchronizePeakAnnotations()
  {
    #ifdef DEBUG_IDENTIFICATION_VIEW    
    std::cout << "synchronizePeakAnnotations." << std::endl;
    #endif

    // Return if no valid peak layer attached
    if (getPeakData() == nullptr 
      || getPeakData()->getMSExperiment().empty() 
      || type != LayerDataBase::DT_PEAK)
    {
      return;
    }

    if (getCurrentSpectrum().getMSLevel() != 2)
    {
      return;
    }

    // no match selected
    const auto selected = getSelectedMatch();
    if (!selected)
    {
      return;
    }

    // copy user annotations to fragment annotation vector
    vector<PeptideHit::PeakAnnotation> fas;
    // do not change the match's annotations, if there are no annotations on the spectrum
    bool annotations_changed(false);
    // for each annotation item on the canvas
    for (auto& a : getCurrentAnnotations())
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
      setPeakAnnotations(getPeakDataMuteable()->getIdentificationData(), *selected, std::move(fas));
    }
  }

  void LayerData1DPeak::removePeakAnnotationsFromPeptideHit(const std::vector<Annotation1DItem*>& selected_annotations)
  {
    // Return if no valid peak layer attached
    if (getPeakData() == nullptr || getPeakData()->getMSExperiment().empty() || type != LayerDataBase::DT_PEAK)
    {
      return;
    }

    // wrong MS level
    if (getCurrentSpectrum().getMSLevel() < 2)
    {
      return;
    }

    // the selected match. That this function returns prematurely is unlikely, since we are deleting existing
    // annotations, that have to be somewhere, but better make sure
    const auto selected = getSelectedMatch();
    if (!selected)
    {
      return;
    }

    vector<PeptideHit::PeakAnnotation> fas =
      getPeakData()->getIdentificationData().findRunByUuid(selected->run_uuid)->getMatch(selected->match).peak_annotations;
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
          if (StringUtils::hasPrefix(fromQString(pa->getText()), tmp_a.annotation))
          {
            to_remove.push_back(tmp_a);
          }
        }
      }
    }
    // remove the collected annotations from the match's annotations
    for (auto const& tmp_a : to_remove)
    {
      fas.erase(std::remove(fas.begin(), fas.end(), tmp_a), fas.end());
    }
    if (!to_remove.empty())
    {
      setPeakAnnotations(getPeakDataMuteable()->getIdentificationData(), *selected, std::move(fas));
    }
  }


}// namespace OpenMS
