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

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/LayerStatisticsDialog.h>
#include <ui_LayerStatisticsDialog.h>

#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/METADATA/MetaInfo.h>

using namespace std;

namespace OpenMS
{

  LayerStatisticsDialog::LayerStatisticsDialog(SpectrumWidget * parent) :
    QDialog(parent),
    canvas_(parent->canvas()),
    layer_data_(canvas_->getCurrentLayer()),
    ui_(new Ui::LayerStatisticsDialogTemplate)
  {
    ui_->setupUi(this);
    
    if (layer_data_.type == LayerData::DT_PEAK)
    {
      computePeakStats_();
    }
    else if (layer_data_.type == LayerData::DT_FEATURE)
    {
      computeFeatureStats_();

      // add two rows for charge and quality
      ui_->table_->setRowCount(ui_->table_->rowCount() + 2);
      QTableWidgetItem * item = new QTableWidgetItem();
      item->setText(QString("Charge"));
      ui_->table_->setVerticalHeaderItem(1, item);
      item = new QTableWidgetItem();
      item->setText(QString("Quality"));
      ui_->table_->setVerticalHeaderItem(2, item);

      // add computed charge and quality stats to the table
      item = new QTableWidgetItem();
      item->setText("-");
      ui_->table_->setItem(1, 0, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(min_charge_, 'f', 2));
      ui_->table_->setItem(1, 1, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(max_charge_, 'f', 2));
      ui_->table_->setItem(1, 2, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(avg_charge_, 'f', 2));
      ui_->table_->setItem(1, 3, item);

      item = new QTableWidgetItem();
      item->setText("-");
      ui_->table_->setItem(2, 0, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(min_quality_, 'f', 2));
      ui_->table_->setItem(2, 1, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(max_quality_, 'f', 2));
      ui_->table_->setItem(2, 2, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(avg_quality_, 'f', 2));
      ui_->table_->setItem(2, 3, item);

    }
    else if (layer_data_.type == LayerData::DT_CONSENSUS)
    {
      computeConsensusStats_();

      // add three rows: charge, quality and elements
      ui_->table_->setRowCount(ui_->table_->rowCount() + 3);
      QTableWidgetItem * item = new QTableWidgetItem();
      item->setText(QString("Charge"));
      ui_->table_->setVerticalHeaderItem(1, item);
      item = new QTableWidgetItem();
      item->setText(QString("Quality"));
      ui_->table_->setVerticalHeaderItem(2, item);
      item = new QTableWidgetItem();
      item->setText(QString("Elements"));
      ui_->table_->setVerticalHeaderItem(3, item);

      // add computed charge and quality stats to the table
      item = new QTableWidgetItem();
      item->setText("-");
      ui_->table_->setItem(1, 0, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(min_charge_, 'f', 2));
      ui_->table_->setItem(1, 1, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(max_charge_, 'f', 2));
      ui_->table_->setItem(1, 2, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(avg_charge_, 'f', 2));
      ui_->table_->setItem(1, 3, item);

      item = new QTableWidgetItem();
      item->setText("-");
      ui_->table_->setItem(2, 0, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(min_quality_, 'f', 2));
      ui_->table_->setItem(2, 1, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(max_quality_, 'f', 2));
      ui_->table_->setItem(2, 2, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(avg_quality_, 'f', 2));
      ui_->table_->setItem(2, 3, item);

      item = new QTableWidgetItem();
      item->setText("-");
      ui_->table_->setItem(3, 0, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(min_elements_, 'f', 2));
      ui_->table_->setItem(3, 1, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(max_elements_, 'f', 2));
      ui_->table_->setItem(3, 2, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(avg_elements_, 'f', 2));
      ui_->table_->setItem(3, 3, item);

    }
    else if (layer_data_.type == LayerData::DT_CHROMATOGRAM)
    {
      //TODO CHROM
    }
    // add computed intensity stats to the table
    QTableWidgetItem * item = new QTableWidgetItem();
    item->setText("-");
    ui_->table_->setItem(0, 0, item);
    item = new QTableWidgetItem();
    item->setText(QString::number(min_intensity_, 'f', 2));
    ui_->table_->setItem(0, 1, item);
    item = new QTableWidgetItem();
    item->setText(QString::number(max_intensity_, 'f', 2));
    ui_->table_->setItem(0, 2, item);
    item = new QTableWidgetItem();
    item->setText(QString::number(avg_intensity_, 'f', 2));
    ui_->table_->setItem(0, 3, item);
    QPushButton * button = new QPushButton("intensity", ui_->table_);
    ui_->table_->setCellWidget(0, 4, button);
    connect(button, SIGNAL(clicked()), this, SLOT(showDistribution_()));


    // add computed stats about meta infos in the FloatDataArrays of the spectra to the table
    for (std::map<String, MetaStatsValue_>::const_iterator it = meta_array_stats_.begin(); it != meta_array_stats_.end(); ++it)
    {
      ui_->table_->setRowCount(ui_->table_->rowCount() + 1);
      String name = it->first;

      item = new QTableWidgetItem();
      item->setText(name.toQString());
      ui_->table_->setVerticalHeaderItem(ui_->table_->rowCount() - 1, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(it->second.count));
      ui_->table_->setItem(ui_->table_->rowCount() - 1, 0, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(it->second.min, 'f', 2));
      ui_->table_->setItem(ui_->table_->rowCount() - 1, 1, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(it->second.max, 'f', 2));
      ui_->table_->setItem(ui_->table_->rowCount() - 1, 2, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(it->second.avg, 'f', 2));
      ui_->table_->setItem(ui_->table_->rowCount() - 1, 3, item);

      if (it->second.count >= 2 && it->second.min < it->second.max)
      {
        button = new QPushButton(name.toQString(), ui_->table_);
        ui_->table_->setCellWidget(ui_->table_->rowCount() - 1, 4, button);
        connect(button, SIGNAL(clicked()), this, SLOT(showDistribution_()));
      }
    }

    // add peak/featurewise collected meta stats to the table
    String name;
    for (MetaIterator_ it = meta_stats_.begin(); it != meta_stats_.end(); it++)
    {
      ui_->table_->setRowCount(ui_->table_->rowCount() + 1);
      name = MetaInfo::registry().getName(it->first);

      item = new QTableWidgetItem();
      item->setText(name.toQString());
      ui_->table_->setVerticalHeaderItem(ui_->table_->rowCount() - 1, item);

      item = new QTableWidgetItem();
      item->setText(QString::number(it->second.count));
      ui_->table_->setItem(ui_->table_->rowCount() - 1, 0, item);

      if (it->second.min <= it->second.max)      // if (min <= max) --> value numerical
      {
        item = new QTableWidgetItem();
        item->setText(QString::number(it->second.min, 'f', 2));
        ui_->table_->setItem(ui_->table_->rowCount() - 1, 1, item);

        item = new QTableWidgetItem();
        item->setText(QString::number(it->second.max, 'f', 2));
        ui_->table_->setItem(ui_->table_->rowCount() - 1, 2, item);

        item = new QTableWidgetItem();
        item->setText(QString::number(it->second.avg, 'f', 2));
        ui_->table_->setItem(ui_->table_->rowCount() - 1, 3, item);

        if (it->second.count >= 2 && it->second.min < it->second.max)
        {
          button = new QPushButton(name.toQString(), ui_->table_);
          ui_->table_->setCellWidget(ui_->table_->rowCount() - 1, 4, button);
          connect(button, SIGNAL(clicked()), this, SLOT(showDistribution_()));
        }
      }
      else       // min > max --> meta value was not numerical --> statistics only about the count
      {
        item = new QTableWidgetItem();
        item->setText("-");
        ui_->table_->setItem(ui_->table_->rowCount() - 1, 1, item);
        item = new QTableWidgetItem();
        item->setText("-");
        ui_->table_->setItem(ui_->table_->rowCount() - 1, 2, item);
        item = new QTableWidgetItem();
        item->setText("-");
        ui_->table_->setItem(ui_->table_->rowCount() - 1, 3, item);
      }
    }
  }
  
  LayerStatisticsDialog::~LayerStatisticsDialog()
  {
    delete ui_;
  }


  void LayerStatisticsDialog::computePeakStats_()
  {
    min_intensity_ = canvas_->getCurrentMinIntensity();
    max_intensity_ = canvas_->getCurrentMaxIntensity();
    avg_intensity_ = 0;
    unsigned long divisor = 0;
    for (LayerData::ExperimentType::ConstIterator it_rt = layer_data_.getPeakData()->begin(); it_rt != layer_data_.getPeakData()->end(); it_rt++)
    {
      for (PeakIterator_ it_peak = it_rt->begin(); it_peak != it_rt->end(); it_peak++)
      {
        avg_intensity_ += it_peak->getIntensity();
        divisor++;
      }
      // collect stats about the meta data arrays of this spectrum
      computeMetaDataArrayStats_(it_rt->getFloatDataArrays().begin(), it_rt->getFloatDataArrays().end());
      computeMetaDataArrayStats_(it_rt->getIntegerDataArrays().begin(), it_rt->getIntegerDataArrays().end());
    }
    if (divisor != 0)
      avg_intensity_ /= (double)divisor;
    computeMetaAverages_();
  }

  void LayerStatisticsDialog::computeFeatureStats_()
  {
    min_intensity_ = canvas_->getCurrentMinIntensity();
    max_intensity_ = canvas_->getCurrentMaxIntensity();
    avg_intensity_ = 0;
    if (!layer_data_.getFeatureMap()->empty())
    {
      min_charge_ = layer_data_.getFeatureMap()->begin()->getCharge();
      max_charge_ = layer_data_.getFeatureMap()->begin()->getCharge();
      avg_charge_ = 0;

      min_quality_ = layer_data_.getFeatureMap()->begin()->getOverallQuality();
      max_quality_ = layer_data_.getFeatureMap()->begin()->getOverallQuality();
      avg_quality_ = 0;
    }

    unsigned long divisor = 0;
    for (FeatureIterator_ it = layer_data_.getFeatureMap()->begin(); it != layer_data_.getFeatureMap()->end(); it++)
    {
      if (it->getCharge() < min_charge_)
        min_charge_ = it->getCharge();
      if (it->getCharge() > max_charge_)
        max_charge_ = it->getCharge();
      if (it->getOverallQuality() < min_quality_)
        min_quality_ = it->getOverallQuality();
      if (it->getOverallQuality() > max_quality_)
        max_quality_ = it->getOverallQuality();
      avg_intensity_ += it->getIntensity();
      avg_charge_ += it->getCharge();
      avg_quality_ += it->getOverallQuality();
      divisor++;
      const MetaInfoInterface & mii = static_cast<MetaInfoInterface>(*it);
      bringInMetaStats_(mii);
    }
    if (divisor != 0)
    {
      avg_intensity_ /= (double)divisor;
      avg_charge_ /= (double)divisor;
      avg_quality_ /= (double)divisor;
    }
    computeMetaAverages_();
  }

  void LayerStatisticsDialog::computeConsensusStats_()
  {
    min_intensity_ = canvas_->getCurrentMinIntensity();
    max_intensity_ = canvas_->getCurrentMaxIntensity();
    avg_intensity_ = 0;
    if (!layer_data_.getConsensusMap()->empty())
    {
      min_charge_ = layer_data_.getConsensusMap()->begin()->getCharge();
      max_charge_ = layer_data_.getConsensusMap()->begin()->getCharge();
      avg_charge_ = 0;

      min_quality_ = layer_data_.getConsensusMap()->begin()->getQuality();
      max_quality_ = layer_data_.getConsensusMap()->begin()->getQuality();
      avg_quality_ = 0;

      min_elements_ = layer_data_.getConsensusMap()->begin()->size();
      max_elements_ = layer_data_.getConsensusMap()->begin()->size();
      avg_elements_ = 0;
    }

    unsigned long divisor = 0;
    for (ConsensusIterator_ it = layer_data_.getConsensusMap()->begin(); it != layer_data_.getConsensusMap()->end(); it++)
    {
      if (it->getCharge() < min_charge_)
        min_charge_ = it->getCharge();
      if (it->getCharge() > max_charge_)
        max_charge_ = it->getCharge();
      if (it->getQuality() < min_quality_)
        min_quality_ = it->getQuality();
      if (it->getQuality() > max_quality_)
        max_quality_ = it->getQuality();
      if (it->size() < min_elements_)
        min_elements_ = it->size();
      if (it->size() > max_elements_)
        max_elements_ = it->size();
      avg_intensity_ += it->getIntensity();
      avg_charge_ += it->getCharge();
      avg_quality_ += it->getQuality();
      avg_elements_ += it->size();
      divisor++;
    }
    if (divisor != 0)
    {
      avg_intensity_ /= (double)divisor;
      avg_charge_ /= (double)divisor;
      avg_quality_ /= (double)divisor;
      avg_elements_ /= (double)divisor;
    }
  }

  void LayerStatisticsDialog::bringInMetaStats_(const MetaInfoInterface & meta_interface)
  {
    vector<UInt> new_meta_keys;
    meta_interface.getKeys(new_meta_keys);
    for (vector<UInt>::iterator it_meta_index = new_meta_keys.begin(); it_meta_index != new_meta_keys.end(); ++it_meta_index)
    {
      const DataValue & next_value = meta_interface.getMetaValue(*it_meta_index);
      MetaIterator_ it = meta_stats_.find(*it_meta_index);
      if (it != meta_stats_.end())      // stats about this meta index already exist -> bring this value in
      {
        it->second.count++;
        if (next_value.valueType() == DataValue::INT_VALUE || next_value.valueType() == DataValue::DOUBLE_VALUE)
        {
          double val = (double)next_value;
          if (val < it->second.min)
            it->second.min = val;
          if (val > it->second.max)
            it->second.max = val;
          it->second.avg += val;
        }
      }
      else       // meta index has not occurred before, create new stats for it:
      {
        MetaStatsValue_ meta_stats_value;
        if (next_value.valueType() == DataValue::INT_VALUE || next_value.valueType() == DataValue::DOUBLE_VALUE)
        {
          double val = (double)next_value;
          meta_stats_value = MetaStatsValue_(1, val, val, val);
        }
        else
        {
          meta_stats_value = MetaStatsValue_(1, 1, 0, 0);        // min=1 > max=0 (illegal) indicates that value is not numerical
        }
        meta_stats_.insert(make_pair(*it_meta_index, meta_stats_value));
      }
    }
  }

  void LayerStatisticsDialog::computeMetaAverages_()
  {
    for (MetaIterator_ it = meta_stats_.begin(); it != meta_stats_.end(); it++)
    {
      if (it->second.count != 0)
      {
        it->second.avg /= (double)it->second.count;
      }
    }
    for (std::map<String, MetaStatsValue_>::iterator it = meta_array_stats_.begin(); it != meta_array_stats_.end(); ++it)
    {
      if (it->second.count != 0)
      {
        it->second.avg /= (float)it->second.count;
      }
    }
  }

  void LayerStatisticsDialog::showDistribution_()
  {
    QPushButton * button = qobject_cast<QPushButton *>(sender());
    QString text = button->text();

    if (text == "intensity")
    {
      qobject_cast<SpectrumWidget *>(parent())->showIntensityDistribution();
    }
    else
    {
      qobject_cast<SpectrumWidget *>(parent())->showMetaDistribution(String(text));
    }
  }

} // namespace
