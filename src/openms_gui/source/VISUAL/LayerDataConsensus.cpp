// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerDataConsensus.h> 

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/KERNEL/DimMapper.h>

#include <OpenMS/VISUAL/Painter2DBase.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

using namespace std;

namespace OpenMS
{
  /// Default constructor
  LayerDataConsensus::LayerDataConsensus(ConsensusMapSharedPtrType& map) : LayerDataBase(LayerDataBase::DT_CONSENSUS)
  {
    consensus_map_ = map;
  }

  std::unique_ptr<Painter2DBase> LayerDataConsensus::getPainter2D() const
  {
    return make_unique<Painter2DConsensus>(this);
  }

  std::unique_ptr<LayerStoreData> LayerDataConsensus::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = make_unique<LayerStoreDataConsensusMapVisible>();
    ret->storeVisibleCM(*consensus_map_.get(), visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerStoreData> LayerDataConsensus::storeFullData() const
  {
    auto ret = make_unique<LayerStoreDataConsensusMapAll>();
    ret->storeFullCM(*consensus_map_.get());
    return ret;
  }

  PeakIndex LayerDataConsensus::findHighestDataPoint(const RangeAllType& area) const
  {
    using IntType = MSExperiment::ConstAreaIterator::PeakType::IntensityType;
    auto max_int = numeric_limits<IntType>::lowest();
    PeakIndex max_pi;
    for (ConsensusMapType::ConstIterator i = getConsensusMap()->begin(); i != getConsensusMap()->end(); ++i)
    {
      // consensus feature in visible area?
      if (area.containsRT(i->getRT()) && area.containsMZ(i->getMZ()) && filters.passes(*i))
      {
        if (i->getIntensity() > max_int)
        {
          max_int = i->getIntensity();
          max_pi = PeakIndex(i - getConsensusMap()->begin());
        }
      }
    }
    return max_pi;
  }

  PointXYType LayerDataConsensus::peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const
  {
    return mapper.map(peak.getFeature(*getConsensusMap()));
  }


  /*std::unique_ptr<Painter1DBase> LayerDataConsensus::getPainter1D() const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  } */
  
  std::unique_ptr<LayerStatistics> LayerDataConsensus::getStats() const
  {
    return make_unique<LayerStatisticsConsensusMap>(*consensus_map_);
  }

  bool LayerDataConsensus::annotate(const vector<PeptideIdentification>& identifications, const vector<ProteinIdentification>& protein_identifications)
  {
    IDMapper mapper;
    mapper.annotate(*getConsensusMap(), identifications, protein_identifications);

    return true;
  }

}// namespace OpenMS
