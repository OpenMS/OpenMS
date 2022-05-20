// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerDataPeak.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/VISUAL/Painter1DBase.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerVisibleData.h>

using namespace std;

namespace OpenMS
{
  LayerDataPeak::LayerDataPeak() : LayerDataBase(LayerDataBase::DT_PEAK)
  {
    flags.set(LayerDataBase::P_PRECURSORS);
  }

  std::unique_ptr<LayerVisibleData> LayerDataPeak::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = std::unique_ptr<LayerVisibleDataPeakMap>();
    ret->storeVisibleExperiment(*peak_map_.get(), visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerVisibleData> LayerDataPeak::storeFullData() const
  {
    auto ret = std::unique_ptr<LayerFullDataPeakMap>();
    ret->storeFullExperiment(*peak_map_.get());
    return ret;
  }

  LayerDataPeak::ProjectionData LayerDataPeak::getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const
  {
    ProjectionData result;

    // create projection data
    map<float, float> rt;
    map<int, float> mzint;
    map<int, float> mzsum;

    auto& peak_count = result.stats.number_of_datapoints;
    auto& intensity_max = result.stats.max_intensity;
    double total_intensity_sum = 0.0; 

    // divide visible range into 100 bins (much faster than using a constant, e.g. 0.05, leading to many peaks for large maps without more information)
    float mz_range = area.RangeMZ::getSpan();
    float mult = 100.0f / (std::isnan(mz_range) ? 1 : mz_range);

    MSSpectrum projection_mz;
    MSChromatogram projection_rt;

    for (auto i = getPeakData()->areaBeginConst(area.getMinRT(), area.getMaxRT(), area.getMinMZ(), area.getMaxMZ()); i != getPeakData()->areaEndConst(); ++i)
    {
      PeakIndex pi = i.getPeakIndex();
      if (filters.passes((*getPeakData())[pi.spectrum], pi.peak))
      {
        // summary stats
        ++peak_count;
        total_intensity_sum += i->getIntensity();
        intensity_max = max(intensity_max, i->getIntensity());

        // binning for m/z
        auto intensity = i->getIntensity();
        mzint[int(i->getMZ() * mult)] += intensity;
        // ... to later obtain an intensity weighted average m/z value
        mzsum[int(i->getMZ() * mult)] += i->getMZ() * intensity;

        // binning in RT (one value per scan)
        rt[i.getRT()] += i->getIntensity();
      }
    }

    // write to spectra/chrom
    projection_mz.resize(mzint.size() + 2);
    projection_mz[0].setMZ(area.getMinMZ());
    projection_mz[0].setIntensity(0.0);
    projection_mz.back().setMZ(area.getMaxMZ());
    projection_mz.back().setIntensity(0.0);

    projection_rt.resize(rt.size() + 2);
    projection_rt[0].setRT(area.getMinRT());
    projection_rt[0].setIntensity(0.0);
    projection_rt.back().setRT(area.getMaxRT());
    projection_rt.back().setIntensity(0.0);

    Size i = 1;
    map<int, float>::iterator intit = mzint.begin();

    for (auto it = mzsum.cbegin(); it != mzsum.cend(); ++it)
    {
      auto intensity = intit->second;
      projection_mz[i].setMZ(it->second / intensity);
      projection_mz[i].setIntensity(intensity);
      ++intit;
      ++i;
    }

    i = 1;
    for (auto it = rt.cbegin(); it != rt.cend(); ++it)
    {
      projection_rt[i].setMZ(it->first);
      projection_rt[i].setIntensity(it->second);
      ++i;
    }

    // create final datastructure

    // projection for m/z
    auto ptr_mz = make_unique<LayerData1DPeak>();
    MSExperiment exp_mz;
    exp_mz.addSpectrum(std::move(projection_mz));
    ptr_mz->setPeakData(ExperimentSharedPtrType(new ExperimentType(exp_mz)));

    // projection for RT
    auto ptr_rt = make_unique<LayerData1DChrom>();
    MSExperiment exp_rt;
    exp_mz.addChromatogram(std::move(projection_rt));
    ptr_rt->setPeakData(ExperimentSharedPtrType(new ExperimentType(exp_rt)));

    auto assign_axis = [&](auto unit, auto& layer) {
      switch (unit)
      {
        case DIM_UNIT::MZ:
          layer = std::move(ptr_mz);
          break;
        case DIM_UNIT::RT:
          layer = std::move(ptr_rt);
          break;
        default:
          // do nothing, leave projection empty
          break;
      }
    };
    assign_axis(unit_x, result.projection_ontoX);
    assign_axis(unit_y, result.projection_ontoY);

    return result;
  }

  LayerDataDefs::PointXYType LayerDataPeak::peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const
  {
    return mapper.map(getSpectrum(peak.spectrum)[peak.peak]);
  }

  std::unique_ptr<LayerStatistics> LayerDataPeak::getStats() const
  {
    return make_unique<LayerStatisticsPeakMap>(*peak_map_);
  }

  std::unique_ptr<LayerVisibleData> LayerData1DPeak::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = std::unique_ptr<LayerVisibleDataPeakMap>();
    ret->storeVisibleSpectrum(getCurrentSpectrum(), visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerVisibleData> LayerData1DPeak::storeFullData() const
  {
    return LayerDataPeak::storeFullData(); // just forward
  }

  QMenu* LayerData1DPeak::getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint)
  {
    QMenu* context_menu = new QMenu("Peak1D", nullptr);
    context_menu->addAction("Edit", [&]() {
      annot_item->editText();
      synchronizePeakAnnotations();
      need_repaint = true;
    });
    context_menu->addAction("Delete", [&]() {
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
    Peak1D peak_lt, peak_rb;
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
    const auto center_intensity = peak_lt.getIntensity() + peak_rb.getIntensity();
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
    PeakType peak = getCurrentSpectrum()[peak_index.peak];
    auto* item = new Annotation1DPeakItem<PeakType>(peak, text, color);
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
