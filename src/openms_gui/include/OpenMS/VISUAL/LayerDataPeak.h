// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/LayerData1DBase.h>

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>

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
    /// Copy-ctor
    LayerDataPeak(const LayerDataPeak& ld) = default;
    /// no assignment operator (should not be needed)
    LayerDataPeak& operator=(const LayerDataPeak& ld) = delete;

    std::unique_ptr<Painter2DBase> getPainter2D() const override;

    std::unique_ptr<LayerData1DBase> to1DLayer() const override;

    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;

    std::unique_ptr<LayerStoreData> storeFullData() const override;

    ProjectionData getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const override;

    PeakIndex findHighestDataPoint(const RangeAllType& area) const override;

    void updateRanges() override
    {
      peak_map_->updateRanges();
      // on_disc_peaks_->updateRanges(); // note: this is not going to work since its on disk! We currently don't have a good way to access these ranges
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      r.assign(*peak_map_);
      return r;
    }

    PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const override;

    String getDataArrayDescription(const PeakIndex& peak_index) override;

    std::unique_ptr<LayerStatistics> getStats() const override;

    bool annotate(const std::vector<PeptideIdentification>& identifications, const std::vector<ProteinIdentification>& protein_identifications) override;

    const ExperimentType::SpectrumType& getSpectrum(Size spectrum_idx) const
    {
      if ((*peak_map_)[spectrum_idx].size() > 0)
      {
        return (*peak_map_)[spectrum_idx];
      }
      if (!on_disc_peaks_->empty())
      {
        static MSSpectrum local_spec;
        local_spec = on_disc_peaks_->getSpectrum(spectrum_idx);
        return local_spec;
      }
      return (*peak_map_)[spectrum_idx];
    }

    /**
    @brief Returns a const reference to the current in-memory peak data

    @note Depending on the caching strategy (on-disk or in-memory), all or some
    spectra may have zero size and contain only meta data since peak data is
    cached on disk.

    @note Do *not* use this function to access the current spectrum for the 1D view, use getCurrentSpectrum() instead.
    */
    const ConstExperimentSharedPtrType getPeakData() const;

    /**
    @brief Returns a mutable reference to the current in-memory peak data

    @note Depending on the caching strategy (on-disk or in-memory), all or some
    spectra may have zero size and contain only meta data since peak data is
    cached on disk.

    @note Do *not* use this function to access the current spectrum for the 1D view, use getCurrentSpectrum() instead.
    */
    const ExperimentSharedPtrType& getPeakDataMuteable()
    {
      return peak_map_;
    }

    /**
    @brief Set the current in-memory peak data
    */
    void setPeakData(ExperimentSharedPtrType p)
    {
      peak_map_ = p;
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


    
    /// Check whether the current layer should be represented as ion mobility
    bool isIonMobilityData() const
    {
      return this->getPeakData()->size() > 0 && this->getPeakData()->metaValueExists("is_ion_mobility") && this->getPeakData()->getMetaValue("is_ion_mobility").toBool();
    }

    void labelAsIonMobilityData() const
    {
      peak_map_->setMetaValue("is_ion_mobility", "true");
    }

    /// Check whether the current layer contains DIA (SWATH-MS) data
    bool isDIAData() const
    {
      return this->getPeakData()->size() > 0 && this->getPeakData()->metaValueExists("is_dia_data") && this->getPeakData()->getMetaValue("is_dia_data").toBool();
    }

    /// Label the current layer as DIA (SWATH-MS) data
    void labelAsDIAData()
    {
      peak_map_->setMetaValue("is_dia_data", "true");
    }

    /**
    @brief Check whether the current layer is a chromatogram

    This is needed because type will *not* distinguish properly between
    chromatogram and spectra data. This is due to the fact that we store
    chromatograms for display in 1D in a data layer using MSSpectrum and
    so the layer looks like PEAK data to tools.
    */
    bool chromatogram_flag_set() const
    {
      return this->getPeakData()->size() > 0 && this->getPeakData()->metaValueExists("is_chromatogram") && this->getPeakData()->getMetaValue("is_chromatogram").toBool();
    }

    /// set the chromatogram flag
    void set_chromatogram_flag()
    {
      peak_map_->setMetaValue("is_chromatogram", "true");
    }

    /// remove the chromatogram flag
    void remove_chromatogram_flag()
    {
      if (this->chromatogram_flag_set())
      {
        peak_map_->removeMetaValue("is_chromatogram");
      }
    }


  protected:
    /// peak data
    ExperimentSharedPtrType peak_map_ = ExperimentSharedPtrType(new ExperimentType());

    /// on disc peak data
    ODExperimentSharedPtrType on_disc_peaks_ = ODExperimentSharedPtrType(new OnDiscMSExperiment());
  };

}// namespace OpenMS
