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

#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/LayerData1DBase.h>

namespace OpenMS
{
  class Annotation1DItem;

  /**
  @brief Class that stores the data for one layer of type PeakMap

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataPeak : public virtual LayerDataBase
  {
  public:

    using SpectrumType = ExperimentType::SpectrumType;
    using PeakType = SpectrumType::PeakType;

    /// Default constructor
    LayerDataPeak();
    /// no Copy-ctor (should not be needed)
    LayerDataPeak(const LayerDataPeak& ld) = delete;
    /// no assignment operator (should not be needed)
    LayerDataPeak& operator=(const LayerDataPeak& ld) = delete;
    /// move Ctor
    LayerDataPeak(LayerDataPeak&& ld) = default;
    /// move assignment
    LayerDataPeak& operator=(LayerDataPeak&& ld) = default;

    std::unique_ptr<LayerVisibleData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;
    std::unique_ptr<LayerVisibleData> storeFullData() const override;
    ProjectionData getProjection(const DIM_UNIT unit, const RangeAllType& area) const override;



    void updateRanges() override
    {
      peak_map_->updateRanges();
      // on_disc_peaks->updateRanges(); // note: this is not going to work since its on disk! We currently don't have a good way to access these ranges
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      r.assign(*peak_map_);
      return r;
    }

    const ExperimentType::SpectrumType getSpectrum(Size spectrum_idx) const
    {
      if ((*peak_map_)[spectrum_idx].size() > 0)
      {
        return (*peak_map_)[spectrum_idx];
      }
      else if (!on_disc_peaks->empty())
      {
        return on_disc_peaks->getSpectrum(spectrum_idx);
      }
      return (*peak_map_)[spectrum_idx];
    }

    PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const override;

    String getDataArrayDescription(const PeakIndex& peak_index) override
    {
      String status;
      const ExperimentType::SpectrumType& s = getSpectrum(peak_index.spectrum);
      for (Size m = 0; m < s.getFloatDataArrays().size(); ++m)
      {
        if (peak_index.peak < s.getFloatDataArrays()[m].size())
        {
          status += s.getFloatDataArrays()[m].getName() + ": " + s.getFloatDataArrays()[m][peak_index.peak] + " ";
        }
      }
      for (Size m = 0; m < s.getIntegerDataArrays().size(); ++m)
      {
        if (peak_index.peak < s.getIntegerDataArrays()[m].size())
        {
          status += s.getIntegerDataArrays()[m].getName() + ": " + s.getIntegerDataArrays()[m][peak_index.peak] + " ";
        }
      }
      for (Size m = 0; m < s.getStringDataArrays().size(); ++m)
      {
        if (peak_index.peak < s.getStringDataArrays()[m].size())
        {
          status += s.getStringDataArrays()[m].getName() + ": " + s.getStringDataArrays()[m][peak_index.peak] + " ";
        }
      }
      return status;
    }

    std::unique_ptr<LayerStatistics> getStats() const override;
  };


  
  class OPENMS_GUI_DLLAPI LayerData1DPeak : public LayerData1DBase, public LayerDataPeak
  {
  public:
    LayerData1DPeak()
      : LayerDataBase(DT_PEAK)
    {
    }

    std::unique_ptr<LayerVisibleData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;
    std::unique_ptr<LayerVisibleData> storeFullData() const override;

    std::unique_ptr<Painter1DBase> getPainter1D() const override;

    
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
      return cached_spectrum_;
    }

    void sortCurrentSpectrumByPosition()
    {
      cached_spectrum_.sortByPosition();
    }

    void updateRanges() override
    {
      LayerDataPeak::updateRanges();
      cached_spectrum_.updateRanges();
    }

    RangeAllType getRange() const override
    {
      return RangeAllType().assign(getCurrentSpectrum().getRange());
    }

    // docu in base class
    QMenu* getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint) override;

    PeakIndex findClosestDataPoint(const RangeAllType& area) const override;

    const ExperimentType::SpectrumType getSpectrum(Size spectrum_idx) const
    {
      if (spectrum_idx == current_idx_)
      {
        return cached_spectrum_;
      }
      return LayerDataPeak::getSpectrum(spectrum_idx);
    }

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
    /// Current cached spectrum
    ExperimentType::SpectrumType cached_spectrum_;


  };

}// namespace OpenMS
