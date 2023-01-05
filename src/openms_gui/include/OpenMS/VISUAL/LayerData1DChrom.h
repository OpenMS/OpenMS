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
#include <OpenMS/VISUAL/LayerDataChrom.h>

namespace OpenMS
{
  class OPENMS_GUI_DLLAPI LayerData1DChrom : public LayerDataChrom, public LayerData1DBase
  {
  public:
    LayerData1DChrom() : LayerDataBase(DT_CHROMATOGRAM)
    {
    }

    LayerData1DChrom(const LayerDataChrom& base) : LayerDataBase(base), LayerDataChrom(base)
    {
    }

    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;
    std::unique_ptr<LayerStoreData> storeFullData() const override;

    std::unique_ptr<Painter1DBase> getPainter1D() const override;

    bool hasIndex(Size index) const override
    {
      return index < chromatogram_map_->getNrChromatograms();
    }

    RangeAllType getRangeForArea(const RangeAllType partial_range) const override
    {
      const auto& chrom = getCurrentChrom();
      auto chrom_filtered = MSExperiment::ChromatogramType();
      chrom_filtered.insert(chrom_filtered.begin(), chrom.RTBegin(partial_range.getMinRT()), chrom.RTEnd(partial_range.getMaxRT()));
      chrom_filtered.updateRanges();
      return RangeAllType().assign(chrom_filtered.getRange());
    }

    RangeAllType getRange1D() const override
    {
      return RangeAllType().assign(getCurrentChrom().getRange());
    }

    const ExperimentType::ChromatogramType& getCurrentChrom() const
    {
      return getChromatogram(current_idx_);
    }

    void updateRanges() override
    {
      LayerDataChrom::updateRanges();
    }

    RangeAllType getRange() const override
    {
      // do NOT change the behaviour of getRange() for 1D, since we want the full RT range across all chroms
      // when scrolling in the list of chroms
      return LayerDataChrom::getRange();
    }

    // docu in base class
    QMenu* getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint) override;

    PeakIndex findClosestDataPoint(const RangeAllType& area) const override;

    // docu in base class
    Annotation1DItem* addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color) override;

  protected:
    /// Current cached spectrum
    //ExperimentType::SpectrumType cached_spectrum_;
  };

} // namespace OpenMS
