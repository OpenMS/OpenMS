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