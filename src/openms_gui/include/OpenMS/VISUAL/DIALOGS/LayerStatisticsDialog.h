// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/LayerData.h>

#include <QtWidgets/QDialog>
#include <QtWidgets/QPushButton>

#include <utility>
#include <map>

namespace Ui
{
  class LayerStatisticsDialogTemplate;
}

namespace OpenMS
{
  class SpectrumWidget;
  class SpectrumCanvas;

  /**
      @brief Dialog showing statistics about the data of the current layer

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI LayerStatisticsDialog :
    public QDialog
  {
    Q_OBJECT

public:

    /// Constructor
    LayerStatisticsDialog(SpectrumWidget * parent);

    ~LayerStatisticsDialog();

protected slots:

    /// Shows the distribution according to the clicked button
    void showDistribution_();

protected:

    /**
        @brief Struct representing the statistics about one meta information
    */
    struct MetaStatsValue_
    {
      MetaStatsValue_(unsigned long c = 0, double mi = 0, double ma = 0, double a = 0)
      {
        count = c;
        min = mi;
        max = ma;
        avg = a;
      }

      unsigned long count;
      double min, max, avg;
    };

    /// Iterates over peaks of a spectrum
    typedef LayerData::ExperimentType::SpectrumType::ConstIterator PeakIterator_;
    /// Iterates over features of a feature map
    typedef LayerData::FeatureMapType::ConstIterator FeatureIterator_;
    /// Iterates over features of a feature map
    typedef LayerData::ConsensusMapType::ConstIterator ConsensusIterator_;
    /// Iterates over the meta_stats map
    typedef std::map<UInt, MetaStatsValue_>::iterator MetaIterator_;

    /// Computes the statistics of a peak layer
    void computePeakStats_();
    /// Computes the statistics of a feature layer
    void computeFeatureStats_();
    /// Computes the statistics of a consensus feature layer
    void computeConsensusStats_();
    /// Computes the statistics of all meta data contained in the FloatDataArray or IntegerDataArray of an MSSpectrum
    template <typename MetaDataIterator>
    void computeMetaDataArrayStats_(MetaDataIterator begin, MetaDataIterator end);
    /// Brings the meta values of one @p meta_interface (a peak or feature) into the statistics
    void bringInMetaStats_(const MetaInfoInterface & meta_interface);
    /// Computes the averages of all meta values stored in meta_stats and meta_array_stats
    void computeMetaAverages_();

    /// Map containing the statistics about all meta information of the peaks/features in the layer
    std::map<UInt, MetaStatsValue_> meta_stats_;
    /// Map containing the statistics about the FloatDataArrays of all spectra in this layer
    std::map<String, MetaStatsValue_> meta_array_stats_;
    /// The canvas of the layer
    SpectrumCanvas * canvas_;
    /// The LayerData object we compute statistics about
    const LayerData& layer_data_;
    /// Minimum intensity value
    double min_intensity_;
    /// Maximum intensity value
    double max_intensity_;
    /// Average intensity value
    double avg_intensity_;
    /// Minimum charge value
    double min_charge_;
    /// Maximum charge value
    double max_charge_;
    /// Average charge value
    double avg_charge_;
    /// Minimum quality value
    double min_quality_;
    /// Maximum quality value
    double max_quality_;
    /// Average quality value
    double avg_quality_;
    /// Minimum number of elements (for consensus features only)
    double min_elements_;
    /// Maximum number of elements (for consensus features only)
    double max_elements_;
    /// Average number of elements (for consensus features only)
    double avg_elements_;

private:
    ///Not implemented
    LayerStatisticsDialog();

    Ui::LayerStatisticsDialogTemplate* ui_;

  };

  template <typename MetaDataIterator>
  void LayerStatisticsDialog::computeMetaDataArrayStats_(MetaDataIterator begin, MetaDataIterator end)
  {
    for (MetaDataIterator meta_array_it = begin; meta_array_it != end; meta_array_it++)
    {
      String meta_name = meta_array_it->getName();
      MetaStatsValue_ meta_stats_value;
      std::map<String, MetaStatsValue_>::iterator it = meta_array_stats_.find(meta_name);
      if (it != meta_array_stats_.end())       // stats about this meta name already exist -> bring this value in
      {
        meta_stats_value = it->second;
        for (typename MetaDataIterator::value_type::const_iterator value_it = meta_array_it->begin(); value_it != meta_array_it->end(); value_it++)
        {
          float value = *value_it;
          meta_stats_value.count++;
          if (value < meta_stats_value.min)
          {
            meta_stats_value.min = value;
          }
          else if (value > meta_stats_value.max)
          {
            meta_stats_value.max = value;
          }
          meta_stats_value.avg += value;
        }
        it->second = meta_stats_value;
      }
      else if (meta_array_it->size() > 0)    // meta name has not occurred before, create new stats for it:
      {
        float init_value = *(meta_array_it->begin());
        meta_stats_value = MetaStatsValue_(0, init_value, init_value, 0);
        for (typename MetaDataIterator::value_type::const_iterator value_it = meta_array_it->begin(); value_it != meta_array_it->end(); value_it++)
        {
          float value = *value_it;
          meta_stats_value.count++;
          if (value < meta_stats_value.min)
          {
            meta_stats_value.min = value;
          }
          else if (value > meta_stats_value.max)
          {
            meta_stats_value.max = value;
          }
          meta_stats_value.avg += value;
        }
        meta_array_stats_.insert(make_pair(meta_name, meta_stats_value));
      }
    }
  }

}
