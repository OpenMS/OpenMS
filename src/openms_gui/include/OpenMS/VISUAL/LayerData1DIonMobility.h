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
#include <OpenMS/VISUAL/LayerDataIonMobility.h>

namespace OpenMS
{
  
  class OPENMS_GUI_DLLAPI LayerData1DIonMobility : public LayerDataIonMobility, public LayerData1DBase
  {
  public:
    LayerData1DIonMobility()
      : LayerDataBase(DT_PEAK)
    {
    }

    LayerData1DIonMobility(const LayerDataIonMobility& base) : LayerDataBase(base), LayerDataIonMobility(base)
    {
    }


    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;
    std::unique_ptr<LayerStoreData> storeFullData() const override;

    std::unique_ptr<Painter1DBase> getPainter1D() const override;

    bool hasIndex(Size index) const override
    {
      return index == 0;
    }

    RangeAllType getRangeForArea(const RangeAllType partial_range) const override
    {
      const auto& spec = getCurrentMobilogram();
      auto spec_filtered = Mobilogram();
      spec_filtered.insert(spec_filtered.begin(), spec.MBBegin(partial_range.getMinMobility()), spec.MBEnd(partial_range.getMaxMobility()));
      spec_filtered.updateRanges();
      return RangeAllType().assign(spec_filtered.getRange());
    }
    
    RangeAllType getRange1D() const override
    {
      return RangeAllType().assign(getCurrentMobilogram().getRange());
    }

    const Mobilogram& getCurrentMobilogram() const
    {
      return LayerDataIonMobility::getMobilogram(this->getCurrentIndex());
    }

    void updateRanges() override
    {
      LayerDataIonMobility::updateRanges();
    }

    RangeAllType getRange() const override
    {
      // do NOT change the behaviour of getRange() for 1D, since we want the full IM range across all mbs
      // when scrolling in the list of mbs
      return LayerDataIonMobility::getRange();
    }

    // docu in base class
    QMenu* getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint) override;

    PeakIndex findClosestDataPoint(const RangeAllType& area) const override;

    // docu in base class
    Annotation1DItem* addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color) override;

  protected:
  };

}// namespace OpenMS
