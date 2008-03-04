// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/LayerStatisticsDialog.h>

using namespace std;

namespace OpenMS
{

	LayerStatisticsDialog::LayerStatisticsDialog(SpectrumWidget* parent)
		: QDialog(parent)
	{
		setupUi(this);
		
		QPushButton* button = new QPushButton("Show", tableWidget);
		tableWidget->setCellWidget(0,4,button);
		connect(button, SIGNAL(clicked()), parent, SLOT(showIntensityDistribution()));
		
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
			tableWidget->setRowCount(tableWidget->rowCount() + 2);
			QTableWidgetItem* item = new QTableWidgetItem();
			item->setText(QString("Charge"));
			tableWidget->setVerticalHeaderItem(1, item);
			item = new QTableWidgetItem();
			item->setText(QString("Quality"));
			tableWidget->setVerticalHeaderItem(2, item);
			
			// add computed charge and quality stats to the table
			item = new QTableWidgetItem();
			item->setText(QString::number(min_charge_,'f',3));
			tableWidget->setItem(1,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_charge_,'f',3));
			tableWidget->setItem(1,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_charge_,'f',3));
			tableWidget->setItem(1,3,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(min_quality_,'f',3));
			tableWidget->setItem(2,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_quality_,'f',3));
			tableWidget->setItem(2,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_quality_,'f',3));
			tableWidget->setItem(2,3,item);
			
		}
		
		// add computed intensity stats to the table
		QTableWidgetItem* item = new QTableWidgetItem();
		item->setText(QString::number(min_intensity_,'f',3));
		tableWidget->setItem(0,1,item);
		item = new QTableWidgetItem();
		item->setText(QString::number(max_intensity_,'f',3));
		tableWidget->setItem(0,2,item);
		item = new QTableWidgetItem();
		item->setText(QString::number(avg_intensity_,'f',3));
		tableWidget->setItem(0,3,item);
		
		// add all computed meta stats to the table
		String name;
		MetaStatsValue_* meta_stats_value;
		for(MetaIterator_ it = meta_stats_.begin(); it != meta_stats_.end(); it++)
		{
			tableWidget->setRowCount(tableWidget->rowCount()+1);
			name = MetaInfo::registry().getName(it->first);
			meta_stats_value = it->second;
			
			item = new QTableWidgetItem();
			item->setText(name.toQString());
			tableWidget->setVerticalHeaderItem(tableWidget->rowCount()-1, item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(meta_stats_value->count));
			tableWidget->setItem(tableWidget->rowCount()-1, 0, item);
			
			if(meta_stats_value->min <= meta_stats_value->max) // if (min <= max) --> value numerical
			{
				item = new QTableWidgetItem();
				item->setText(QString::number(meta_stats_value->min,'f',3));
				tableWidget->setItem(tableWidget->rowCount()-1, 1, item);
				
				item = new QTableWidgetItem();
				item->setText(QString::number(meta_stats_value->max,'f',3));
				tableWidget->setItem(tableWidget->rowCount()-1, 2, item);
				
				item = new QTableWidgetItem();
				item->setText(QString::number(meta_stats_value->avg,'f',3));
				tableWidget->setItem(tableWidget->rowCount()-1, 3, item);
			}
			else // min > max --> meta value was not numerical --> statistics only about the count
			{
				item = new QTableWidgetItem();
				item->setText("-");
				tableWidget->setItem(tableWidget->rowCount()-1, 1, item);
				tableWidget->setItem(tableWidget->rowCount()-1, 2, item);
				tableWidget->setItem(tableWidget->rowCount()-1, 3, item);
			}
		}
	}
	
	void LayerStatisticsDialog::computePeakStats_()
	{
		min_intensity_ = canvas_->getCurrentMinIntensity();
		max_intensity_ = canvas_->getCurrentMaxIntensity();
		avg_intensity_ = 0;
		unsigned long divisor = 0;
		for(RTIterator_ it_rt = layer_data_.peaks.begin(); it_rt != layer_data_.peaks.end(); it_rt++)
		{
			for(PeakIterator_ it_peak = it_rt->begin(); it_peak != it_rt->end(); it_peak++)
			{
				avg_intensity_ += it_peak->getIntensity();
				divisor++;
				const MetaInfoInterface& mii = static_cast<MetaInfoInterface>(*it_peak);
				bringInMetaStats_(mii);
			}
		}
		if (divisor != 0) avg_intensity_ /= (DoubleReal)divisor;
		computeMetaAverages_();
	}
	
	
	void LayerStatisticsDialog::computeFeatureStats_()
	{
		min_intensity_ = canvas_->getCurrentMinIntensity();
		max_intensity_ = canvas_->getCurrentMaxIntensity();
		avg_intensity_ = 0;
		if(!layer_data_.features.empty())
		{
			min_charge_ = layer_data_.features.begin()->getCharge();
			max_charge_ = layer_data_.features.begin()->getCharge();
			avg_charge_ = 0;
		
			min_quality_ = layer_data_.features.begin()->getOverallQuality();
			max_quality_ = layer_data_.features.begin()->getOverallQuality();
			avg_quality_ = 0;
		}
		
		unsigned long divisor = 0;
		for(FeatureIterator_ it_feature = layer_data_.features.begin(); it_feature != layer_data_.features.end(); it_feature++)
		{
			if(it_feature->getCharge() < min_charge_) min_charge_ = it_feature->getCharge();
			if(it_feature->getCharge() > max_charge_) max_charge_ = it_feature->getCharge();
			if(it_feature->getOverallQuality() < min_quality_) min_quality_ = it_feature->getOverallQuality();
			if(it_feature->getOverallQuality() > max_quality_) max_quality_ = it_feature->getOverallQuality();
			avg_intensity_ += it_feature->getIntensity();
			avg_charge_ += it_feature->getCharge();
			avg_quality_ += it_feature->getOverallQuality();
			divisor++;
			const MetaInfoInterface& mii = static_cast<MetaInfoInterface>(*it_feature);
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
	
	
	void LayerStatisticsDialog::bringInMetaStats_(const MetaInfoInterface& meta_interface)
	{
		vector<UInt> new_meta_keys;
		DataValue next_value;
		meta_interface.getKeys(new_meta_keys);
		for(vector<UInt>::iterator it_meta_index = new_meta_keys.begin(); it_meta_index != new_meta_keys.end(); it_meta_index++)
		{
			next_value = meta_interface.getMetaValue(*it_meta_index);
			MetaIterator_ it = meta_stats_.find(*it_meta_index);
			if(it != meta_stats_.end()) // stats about this meta index already exist -> bring this value in
			{
				MetaStatsValue_* meta_stats_value = it->second;
				meta_stats_value->count++;
				if(next_value.valueType() == DataValue::INT_VALUE || next_value.valueType() == DataValue::DOUBLE_VALUE)
				{
					if((DoubleReal)next_value < meta_stats_value->min) meta_stats_value->min = (DoubleReal)next_value;
					if((DoubleReal)next_value > meta_stats_value->max) meta_stats_value->max = (DoubleReal)next_value;
					meta_stats_value->avg += (DoubleReal)next_value;
				}
			}
			else // meta index has not occurred before, create new stats for it:
			{
				MetaStatsValue_* meta_stats_value;
				if(next_value.valueType() == DataValue::INT_VALUE || next_value.valueType() == DataValue::DOUBLE_VALUE)
				{
					DoubleReal val = (DoubleReal)next_value;
					meta_stats_value = new MetaStatsValue_(1,val,val,val);
				}
				else meta_stats_value = new MetaStatsValue_(1,1,0,0); // min=1 > max=0 (illegal) indicates that value is not numerical
				pair<UInt, MetaStatsValue_*> p = make_pair(*it_meta_index, meta_stats_value);
				meta_stats_.insert(p);
			}
		}
	}
	
	
	void LayerStatisticsDialog::computeMetaAverages_()
	{
		for(MetaIterator_ it = meta_stats_.begin(); it != meta_stats_.end(); it++)
		{
			MetaStatsValue_* meta_stats_value = it->second;
			if(meta_stats_value->count != 0)
			{
				meta_stats_value->avg /= (DoubleReal)meta_stats_value->count;
			}
		}
	}
	
} // namespace
