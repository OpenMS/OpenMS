// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/LayerDataIonMobility.h>

#include <OpenMS/KERNEL/DimMapper.h>

#include <OpenMS/VISUAL/Painter2DBase.h>
#include <OpenMS/VISUAL/LayerData1DIonMobility.h>


using namespace std;

namespace OpenMS
{
  LayerDataIonMobility::LayerDataIonMobility() : LayerDataBase(LayerDataBase::DT_PEAK)
  {
  }

  LayerDataIonMobility::LayerDataIonMobility(const LayerDataIonMobility& ld)
    : LayerDataBase(static_cast<const LayerDataBase&>(ld))
  {
  }

  std::unique_ptr<Painter2DBase> LayerDataIonMobility::getPainter2D() const
  {
    return make_unique<Painter2DIonMobility>(this);
  }


  std::unique_ptr<LayerData1DBase> LayerDataIonMobility::to1DLayer() const
  {
    return make_unique<LayerData1DIonMobility>(*this);
  }

  std::unique_ptr<LayerStoreData> LayerDataIonMobility::storeVisibleData(const RangeAllType& /*visible_range*/, const DataFilters& /*layer_filters*/) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    // does not exist yet...
    /*auto ret = make_unique<LayerStoreDataMobilogramVisible>();
    ret->storeVisibleMobilogram(single_mobilogram_, visible_range, layer_filters);
    return ret;*/
  }

  std::unique_ptr<LayerStoreData> LayerDataIonMobility::storeFullData() const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    // does not exist yet...
    /*auto ret = make_unique<LayerStoreDataMobilogramAll>();
    ret->storeFullMobilograms(*peak_map_.get());
    return ret;*/
  }

  LayerDataIonMobility::ProjectionData LayerDataIonMobility::getProjection(const DIM_UNIT /*unit_x*/, const DIM_UNIT /*unit_y*/, const RangeAllType& /*area*/) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    /*ProjectionData result;
    return result;*/
  }

  PointXYType LayerDataIonMobility::peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const
  {
    if (peak.spectrum != 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Currently only one mobilogram is supported!", String(peak.spectrum));
    }
    return mapper.map(single_mobilogram_, peak.peak);
  }

  String LayerDataIonMobility::getDataArrayDescription(const PeakIndex& peak_index)
  {
    if (peak_index.spectrum != 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Currently only one mobilogram is supported!", String(peak_index.spectrum));
    }
    // no array data exists for mobilograms (yet)
    String status;
    return status;
  }

  std::unique_ptr<LayerStatistics> LayerDataIonMobility::getStats() const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    // does not exist yet... 
    //return make_unique<LayerStatisticsIonMobilityMap>(single_mobilogram_);
  }

} // namespace OpenMS