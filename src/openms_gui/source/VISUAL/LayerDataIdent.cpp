// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerDataIdent.h>

#include <OpenMS/KERNEL/DimMapper.h>

#include <OpenMS/VISUAL/Painter2DBase.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

using namespace std;

namespace OpenMS
{
  std::unique_ptr<Painter2DBase> LayerDataIdent::getPainter2D() const
  {
    return make_unique<Painter2DIdent>(this);
  }

  std::unique_ptr<LayerStoreData> LayerDataIdent::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = make_unique<LayerStoreDataIdentVisible>();
    ret->storeVisibleIdent(peptides_, visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerStoreData> LayerDataIdent::storeFullData() const
  {
    auto ret = make_unique<LayerStoreDataIdentAll>();
    ret->storeFullIdent(peptides_);
    return ret;
  }

  LayerDataDefs::ProjectionData LayerDataIdent::getProjection(const DIM_UNIT /*unit_x*/, const DIM_UNIT /*unit_y*/, const RangeAllType& /*area*/) const
  { // currently only a stub
    ProjectionData proj;
    return proj;
  }

  PointXYType LayerDataIdent::peakIndexToXY(const PeakIndex& /*peak*/, const DimMapper<2>& /*mapper*/) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  std::unique_ptr<LayerStatistics> LayerDataIdent::getStats() const
  {
    return make_unique<LayerStatisticsIdent>(peptides_);
  }
} // namespace OpenMS
