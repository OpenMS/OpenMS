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

#include <vector>

namespace OpenMS
{

  /**
  @brief Class that stores the data for one layer of type ConsensusMap

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataConsensus : public virtual LayerDataBase
  {
  public:
    /// Default constructor
    LayerDataConsensus(ConsensusMapSharedPtrType& map);
    /// no Copy-ctor (should not be needed)
    LayerDataConsensus(const LayerDataConsensus& ld) = delete;
    /// no assignment operator (should not be needed)
    LayerDataConsensus& operator=(const LayerDataConsensus& ld) = delete;

    std::unique_ptr<Painter2DBase> getPainter2D() const override;

    std::unique_ptr<LayerData1DBase> to1DLayer() const override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;

    std::unique_ptr<LayerStoreData> storeFullData() const override;

    ProjectionData getProjection(const DIM_UNIT /*unit_x*/, const DIM_UNIT /*unit_y*/, const RangeAllType& /*area*/) const override
    { // currently only a stub
      ProjectionData proj;
      return proj;
    }

    PeakIndex findHighestDataPoint(const RangeAllType& area) const override;

    PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const override;

    void updateRanges() override
    {
      consensus_map_->updateRanges();
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      r.assign(*getConsensusMap());
      return r;
    }

    std::unique_ptr<LayerStatistics> getStats() const override;

    bool annotate(const std::vector<PeptideIdentification>& identifications, const std::vector<ProteinIdentification>& protein_identifications) override;
  
    /// Returns a const reference to the consensus feature data
    const ConsensusMapSharedPtrType& getConsensusMap() const
    {
      return consensus_map_;
    }

    /// Returns current consensus map (mutable)
    ConsensusMapSharedPtrType& getConsensusMap()
    {
      return consensus_map_;
    }
  protected:
    /// consensus feature data
    ConsensusMapSharedPtrType consensus_map_ = ConsensusMapSharedPtrType(new ConsensusMapType());
  };

}// namespace OpenMS
