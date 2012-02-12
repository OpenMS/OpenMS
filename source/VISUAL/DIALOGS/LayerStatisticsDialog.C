// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/LayerStatisticsDialog.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>

using namespace std;

namespace OpenMS
{

	LayerStatisticsDialog::LayerStatisticsDialog(SpectrumWidget* parent)
		: QDialog(parent)
	{
		setupUi(this);
				
		canvas_ = parent->canvas();
		layer_data_ = canvas_->getCurrentLayer();
		
		if(layer_data_.type == LayerData::DT_PEAK)
		{			
			computePeakStats_();
		}
		else if(layer_data_.type == LayerData::DT_FEATURE)
		{
			computeFeatureStats_();
			
			// add two rows for charge and quality
			table_->setRowCount(table_->rowCount() + 2);
			QTableWidgetItem* item = new QTableWidgetItem();
			item->setText(QString("Charge"));
			table_->setVerticalHeaderItem(1, item);
			item = new QTableWidgetItem();
			item->setText(QString("Quality"));
			table_->setVerticalHeaderItem(2, item);
			
			// add computed charge and quality stats to the table
			item = new QTableWidgetItem();
			item->setText("-");
			table_->setItem(1,0,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(min_charge_,'f',2));
			table_->setItem(1,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_charge_,'f',2));
			table_->setItem(1,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_charge_,'f',2));
			table_->setItem(1,3,item);
			
			item = new QTableWidgetItem();
			item->setText("-");
			table_->setItem(2,0,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(min_quality_,'f',2));
			table_->setItem(2,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_quality_,'f',2));
			table_->setItem(2,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_quality_,'f',2));
			table_->setItem(2,3,item);
			
		}
		else if (layer_data_.type == LayerData::DT_CONSENSUS)
		{
			computeConsensusStats_();

			// add thres rows for charge, quality and elements
			table_->setRowCount(table_->rowCount() + 3);
			QTableWidgetItem* item = new QTableWidgetItem();
			item->setText(QString("Charge"));
			table_->setVerticalHeaderItem(1, item);
			item = new QTableWidgetItem();
			item->setText(QString("Quality"));
			table_->setVerticalHeaderItem(2, item);
			item = new QTableWidgetItem();
			item->setText(QString("Elements"));
			table_->setVerticalHeaderItem(3, item);
			
			// add computed charge and quality stats to the table
			item = new QTableWidgetItem();
			item->setText("-");
			table_->setItem(1,0,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(min_charge_,'f',2));
			table_->setItem(1,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_charge_,'f',2));
			table_->setItem(1,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_charge_,'f',2));
			table_->setItem(1,3,item);
			
			item = new QTableWidgetItem();
			item->setText("-");
			table_->setItem(2,0,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(min_quality_,'f',2));
			table_->setItem(2,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_quality_,'f',2));
			table_->setItem(2,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_quality_,'f',2));
			table_->setItem(2,3,item);

			item = new QTableWidgetItem();
			item->setText("-");
			table_->setItem(3,0,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(min_elements_,'f',2));
			table_->setItem(3,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_elements_,'f',2));
			table_->setItem(3,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_elements_,'f',2));
			table_->setItem(3,3,item);

		}
		else if (layer_data_.type == LayerData::DT_CHROMATOGRAM)
		{
			//TODO CHROM
		}
		// add computed intensity stats to the table
		QTableWidgetItem* item = new QTableWidgetItem();
		item->setText("-");
		table_->setItem(0,0,item);
		item = new QTableWidgetItem();
		item->setText(QString::number(min_intensity_,'f',2));
		table_->setItem(0,1,item);
		item = new QTableWidgetItem();
		item->setText(QString::number(max_intensity_,'f',2));
		table_->setItem(0,2,item);
		item = new QTableWidgetItem();
		item->setText(QString::number(avg_intensity_,'f',2));
		table_->setItem(0,3,item);
		QPushButton* button = new QPushButton("intensity", table_);
		table_->setCellWidget(0,4,button);
		connect(button, SIGNAL(clicked()), this, SLOT(showDistribution_()));
		
		
		// add computed stats about meta infos in the FloatDataArrays of the spectra to the table
    for(std::map<String, MetaStatsValue_>::const_iterator it = meta_array_stats_.begin(); it != meta_array_stats_.end(); ++it)
		{
			table_->setRowCount(table_->rowCount()+1);
			String name = it->first;
			
			item = new QTableWidgetItem();
			item->setText(name.toQString());
			table_->setVerticalHeaderItem(table_->rowCount()-1, item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(it->second.count));
			table_->setItem(table_->rowCount()-1, 0, item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(it->second.min,'f',2));
			table_->setItem(table_->rowCount()-1, 1, item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(it->second.max,'f',2));
			table_->setItem(table_->rowCount()-1, 2, item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(it->second.avg,'f',2));
			table_->setItem(table_->rowCount()-1, 3, item);
			
			if (it->second.count>=2 && it->second.min<it->second.max)
			{
				button = new QPushButton(name.toQString(), table_);
				table_->setCellWidget(table_->rowCount()-1,4,button);
				connect(button, SIGNAL(clicked()), this, SLOT(showDistribution_()));
			}
		}
		
		// add peak/featurewise collected meta stats to the table
		String name;
		for(MetaIterator_ it = meta_stats_.begin(); it != meta_stats_.end(); it++)
		{
			table_->setRowCount(table_->rowCount()+1);
			name = MetaInfo::registry().getName(it->first);
			
			item = new QTableWidgetItem();
			item->setText(name.toQString());
			table_->setVerticalHeaderItem(table_->rowCount()-1, item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(it->second.count));
			table_->setItem(table_->rowCount()-1, 0, item);
			
			if(it->second.min <= it->second.max) // if (min <= max) --> value numerical
			{
				item = new QTableWidgetItem();
				item->setText(QString::number(it->second.min,'f',2));
				table_->setItem(table_->rowCount()-1, 1, item);
				
				item = new QTableWidgetItem();
				item->setText(QString::number(it->second.max,'f',2));
				table_->setItem(table_->rowCount()-1, 2, item);
				
				item = new QTableWidgetItem();
				item->setText(QString::number(it->second.avg,'f',2));
				table_->setItem(table_->rowCount()-1, 3, item);
				
				if (it->second.count>=2 && it->second.min<it->second.max)
				{
					button = new QPushButton(name.toQString(), table_);
					table_->setCellWidget(table_->rowCount()-1,4,button);
					connect(button, SIGNAL(clicked()), this, SLOT(showDistribution_()));
				}
			}
			else // min > max --> meta value was not numerical --> statistics only about the count
			{
				item = new QTableWidgetItem();
				item->setText("-");
				table_->setItem(table_->rowCount()-1, 1, item);
				item = new QTableWidgetItem();
				item->setText("-");
				table_->setItem(table_->rowCount()-1, 2, item);
				item = new QTableWidgetItem();
				item->setText("-");
				table_->setItem(table_->rowCount()-1, 3, item);
			}
		}
	}
	
	void LayerStatisticsDialog::computePeakStats_()
	{
		min_intensity_ = canvas_->getCurrentMinIntensity();
		max_intensity_ = canvas_->getCurrentMaxIntensity();
		avg_intensity_ = 0;
		unsigned long divisor = 0;
                for(LayerData::ExperimentType::ConstIterator it_rt = layer_data_.getPeakData()->begin(); it_rt != layer_data_.getPeakData()->end(); it_rt++)
		{
			for(PeakIterator_ it_peak = it_rt->begin(); it_peak != it_rt->end(); it_peak++)
			{
				avg_intensity_ += it_peak->getIntensity();
				divisor++;
			}
			// collect stats about the meta data arrays of this spectrum
			computeMetaDataArrayStats_(it_rt->getFloatDataArrays().begin(),it_rt->getFloatDataArrays().end());
			computeMetaDataArrayStats_(it_rt->getIntegerDataArrays().begin(),it_rt->getIntegerDataArrays().end());
		}
		if (divisor != 0) avg_intensity_ /= (DoubleReal)divisor;
		computeMetaAverages_();
	}
	
	void LayerStatisticsDialog::computeFeatureStats_()
	{
		min_intensity_ = canvas_->getCurrentMinIntensity();
		max_intensity_ = canvas_->getCurrentMaxIntensity();
		avg_intensity_ = 0;
    if(!layer_data_.getFeatureMap()->empty())
		{
      min_charge_ = layer_data_.getFeatureMap()->begin()->getCharge();
      max_charge_ = layer_data_.getFeatureMap()->begin()->getCharge();
			avg_charge_ = 0;
		
      min_quality_ = layer_data_.getFeatureMap()->begin()->getOverallQuality();
      max_quality_ = layer_data_.getFeatureMap()->begin()->getOverallQuality();
			avg_quality_ = 0;
		}
		
		unsigned long divisor = 0;
    for(FeatureIterator_ it = layer_data_.getFeatureMap()->begin(); it != layer_data_.getFeatureMap()->end(); it++)
		{
			if(it->getCharge() < min_charge_) min_charge_ = it->getCharge();
			if(it->getCharge() > max_charge_) max_charge_ = it->getCharge();
			if(it->getOverallQuality() < min_quality_) min_quality_ = it->getOverallQuality();
			if(it->getOverallQuality() > max_quality_) max_quality_ = it->getOverallQuality();
			avg_intensity_ += it->getIntensity();
			avg_charge_ += it->getCharge();
			avg_quality_ += it->getOverallQuality();
			divisor++;
			const MetaInfoInterface& mii = static_cast<MetaInfoInterface>(*it);
			bringInMetaStats_(mii);
		}
		if(divisor != 0)
		{
			avg_intensity_ /= (DoubleReal)divisor;
			avg_charge_ /= (DoubleReal)divisor;
			avg_quality_ /= (DoubleReal)divisor;
		}
		computeMetaAverages_();
	}

	void LayerStatisticsDialog::computeConsensusStats_()
	{
		min_intensity_ = canvas_->getCurrentMinIntensity();
		max_intensity_ = canvas_->getCurrentMaxIntensity();
		avg_intensity_ = 0;
    if(!layer_data_.getConsensusMap()->empty())
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
    for(ConsensusIterator_ it = layer_data_.getConsensusMap()->begin(); it != layer_data_.getConsensusMap()->end(); it++)
		{
			if(it->getCharge() < min_charge_) min_charge_ = it->getCharge();
			if(it->getCharge() > max_charge_) max_charge_ = it->getCharge();
			if(it->getQuality() < min_quality_) min_quality_ = it->getQuality();
			if(it->getQuality() > max_quality_) max_quality_ = it->getQuality();
			if(it->size() < min_elements_) min_elements_ = it->size();
			if(it->size() > max_elements_) max_elements_ = it->size();
			avg_intensity_ += it->getIntensity();
			avg_charge_ += it->getCharge();
			avg_quality_ += it->getQuality();
			avg_elements_ += it->size();
			divisor++;
		}
		if(divisor != 0)
		{
			avg_intensity_ /= (DoubleReal)divisor;
			avg_charge_ /= (DoubleReal)divisor;
			avg_quality_ /= (DoubleReal)divisor;
			avg_elements_ /= (DoubleReal)divisor;
		}
	}

	void LayerStatisticsDialog::bringInMetaStats_(const MetaInfoInterface& meta_interface)
	{
		vector<UInt> new_meta_keys;
		meta_interface.getKeys(new_meta_keys);
    for(vector<UInt>::iterator it_meta_index = new_meta_keys.begin(); it_meta_index != new_meta_keys.end(); ++it_meta_index)
		{
			const DataValue& next_value = meta_interface.getMetaValue(*it_meta_index);
			MetaIterator_ it = meta_stats_.find(*it_meta_index);
			if(it != meta_stats_.end()) // stats about this meta index already exist -> bring this value in
			{
				it->second.count++;
				if(next_value.valueType() == DataValue::INT_VALUE || next_value.valueType() == DataValue::DOUBLE_VALUE)
				{
					DoubleReal val = (DoubleReal)next_value;
					if( val< it->second.min) it->second.min = val;
					if( val > it->second.max) it->second.max = val;
					it->second.avg += val;
				}
			}
			else // meta index has not occurred before, create new stats for it:
			{
				MetaStatsValue_ meta_stats_value;
				if(next_value.valueType() == DataValue::INT_VALUE || next_value.valueType() == DataValue::DOUBLE_VALUE)
				{
					DoubleReal val = (DoubleReal)next_value;
					meta_stats_value = MetaStatsValue_(1,val,val,val);
				}
				else
				{
					meta_stats_value = MetaStatsValue_(1,1,0,0); // min=1 > max=0 (illegal) indicates that value is not numerical
				}
				meta_stats_.insert(make_pair(*it_meta_index, meta_stats_value));
			}
		}
	}
	
	
	void LayerStatisticsDialog::computeMetaAverages_()
	{
		for(MetaIterator_ it = meta_stats_.begin(); it != meta_stats_.end(); it++)
		{
			if(it->second.count != 0)
			{
				it->second.avg /= (DoubleReal)it->second.count;
			}
		}
    for(std::map<String, MetaStatsValue_>::iterator it = meta_array_stats_.begin(); it != meta_array_stats_.end(); ++it)
		{
			if(it->second.count != 0)
			{
				it->second.avg /= (Real)it->second.count;
			}
		}
	}
	
	void LayerStatisticsDialog::showDistribution_()
	{
		QPushButton* button = qobject_cast<QPushButton*>(sender());
		QString text = button->text();
		
		if (text=="intensity")
		{
			qobject_cast<SpectrumWidget*>(parent())->showIntensityDistribution();
		}
		else
		{
			qobject_cast<SpectrumWidget*>(parent())->showMetaDistribution(String(text));
		}
	}
	
} // namespace
