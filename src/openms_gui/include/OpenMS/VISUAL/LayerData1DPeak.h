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

#pragma once

#include <OpenMS/VISUAL/LayerData1DBase.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>

namespace OpenMS
{
  
  class OPENMS_GUI_DLLAPI LayerData1DPeak : public LayerDataPeak, public LayerData1DBase
  {
  public:
    LayerData1DPeak()
      : LayerDataBase(DT_PEAK)
    {
    }

    LayerData1DPeak(const LayerDataPeak& base) : LayerDataBase(base), LayerDataPeak(base)
    {
    }


    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;
    std::unique_ptr<LayerStoreData> storeFullData() const override;

    std::unique_ptr<Painter1DBase> getPainter1D() const override;

    bool hasIndex(Size index) const override
    {
      return index < peak_map_->size();
    }

    RangeAllType getRangeForArea(const RangeAllType partial_range) const override
    {
      const auto& spec = getCurrentSpectrum();
      auto spec_filtered = SpectrumType();
      spec_filtered.insert(spec_filtered.begin(), spec.MZBegin(partial_range.getMinMZ()), spec.MZEnd(partial_range.getMaxMZ()));
      spec_filtered.updateRanges();
      return RangeAllType().assign(spec_filtered.getRange());
    }

    const ExperimentType::SpectrumType& getCurrentSpectrum() const
    {
      return LayerDataPeak::getSpectrum(current_idx_);
    }


    void updateRanges() override
    {
      LayerDataPeak::updateRanges();
    }

    RangeAllType getRange1D() const override
    {
      return RangeAllType().assign(getCurrentSpectrum().getRange());
    }

    // docu in base class
    QMenu* getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint) override;

    PeakIndex findClosestDataPoint(const RangeAllType& area) const override;

    // docu in base class
    Annotation1DItem* addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color) override;

    /// updates the PeakAnnotations in the current PeptideHit with manually changed annotations
    /// if no PeptideIdentification or PeptideHit for the spectrum exist, it is generated
    void synchronizePeakAnnotations();

    /// remove peak annotations in the given list from the currently active PeptideHit
    void removePeakAnnotationsFromPeptideHit(const std::vector<Annotation1DItem*>& selected_annotations);

    /// updates the PeakAnnotations in the current PeptideHit with manually changed annotations
    void updatePeptideHitAnnotations_(PeptideHit& hit);

  protected:
  };

}// namespace OpenMS
