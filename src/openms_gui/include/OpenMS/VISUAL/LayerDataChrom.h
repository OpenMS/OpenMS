// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

namespace OpenMS
{
  /// SharedPtr on OSWData
  typedef boost::shared_ptr<OSWData> OSWDataSharedPtrType;

  /**
  @brief Class that stores the data for one layer of type Chromatogram

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataChrom : public virtual LayerDataBase
  {
  public:
    /// Default constructor
    LayerDataChrom();
    /// Copy-ctor
    LayerDataChrom(const LayerDataChrom& ld) = default;
    /// no assignment operator (should not be needed)
    LayerDataChrom& operator=(const LayerDataChrom& ld) = delete;

    std::unique_ptr<Painter2DBase> getPainter2D() const override;

    std::unique_ptr<LayerData1DBase> to1DLayer() const override;

    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;

    std::unique_ptr<LayerStoreData> storeFullData() const override;

    ProjectionData getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const override;

    PeakIndex findHighestDataPoint(const RangeAllType& area) const override;

    void updateRanges() override
    {
      chromatogram_map_->updateRanges();
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      r.assign(*chromatogram_map_);
      return r;
    }

    std::unique_ptr<LayerStatistics> getStats() const override;

    PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const override;

    String getDataArrayDescription(const PeakIndex& peak_index) override;

    const ExperimentType::ChromatogramType& getChromatogram(Size idx) const
    {
      return chromatogram_map_->getChromatogram(idx);
    }

    
    /**
      @brief Set the current in-memory chrom data
    */
    void setChromData(ExperimentSharedPtrType p)
    {
      chromatogram_map_ = p;
    }

    /// Returns a mutable reference to the current chromatogram data
    const ExperimentSharedPtrType& getChromatogramData() const
    {
      return chromatogram_map_;
    }

    /// Returns a mutable reference to the current chromatogram data
    ExperimentSharedPtrType& getChromatogramData()
    {
      return chromatogram_map_;
    }

    /// Set the current on-disc data
    void setOnDiscPeakData(ODExperimentSharedPtrType p)
    {
      on_disc_peaks_ = p;
    }

    /// Returns a mutable reference to the on-disc data
    const ODExperimentSharedPtrType& getOnDiscPeakData() const
    {
      return on_disc_peaks_;
    }
    
    OSWDataSharedPtrType& getChromatogramAnnotation()
    {
      return chrom_annotation_;
    }

    const OSWDataSharedPtrType& getChromatogramAnnotation() const
    {
      return chrom_annotation_;
    }

    /// add annotation from an OSW sqlite file.
    void setChromatogramAnnotation(OSWData&& data);

  protected:
    /// chromatogram data
    ExperimentSharedPtrType chromatogram_map_ = ExperimentSharedPtrType(new ExperimentType());

    /// on disc chrom data
    ODExperimentSharedPtrType on_disc_peaks_ = ODExperimentSharedPtrType(new OnDiscMSExperiment());

    /// Chromatogram annotation data
    OSWDataSharedPtrType chrom_annotation_;
  };

} //namespace

