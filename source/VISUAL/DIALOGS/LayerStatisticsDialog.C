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

#include <iostream>

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
		
		LayerData layer_data = parent->canvas()->getCurrentLayer();
		
		map<String, MetaStatsValue*> meta_stats;
		vector<String> new_meta_keys;
		DataValue new_value;
		QTableWidgetItem* item;
		DoubleReal min_intensity = parent->canvas()->getCurrentMinIntensity();
		DoubleReal max_intensity = parent->canvas()->getCurrentMaxIntensity();
		DoubleReal avg_intensity = 0;
		
		if(layer_data.type == LayerData::DT_PEAK)
		{			
			unsigned long divisor = 0;
			for(RTIterator it_rt = layer_data.peaks.begin(); it_rt != layer_data.peaks.end(); it_rt++)
			{
				for(PeakIterator it_peak = it_rt->begin(); it_peak != it_rt->end(); it_peak++)
				{
					avg_intensity += it_peak->getIntensity();
					divisor++;
					it_peak->getKeys(new_meta_keys);
					for(vector<String>::iterator it_meta = new_meta_keys.begin(); it_meta != new_meta_keys.end(); it_meta++)
					{
						new_value = it_peak->getMetaValue(*it_meta);
						MetaIterator it = meta_stats.find(*it_meta);
						if(it != meta_stats.end())
						{
							MetaStatsValue* meta_stats_value = it->second;
							meta_stats_value->count++;
							if(new_value.valueType() == DataValue::INT_VALUE || new_value.valueType() == DataValue::DOUBLE_VALUE)
							{
								if((DoubleReal)new_value < meta_stats_value->min) meta_stats_value->min = (DoubleReal)new_value;
								if((DoubleReal)new_value > meta_stats_value->max) meta_stats_value->max = (DoubleReal)new_value;
								meta_stats_value->avg += (DoubleReal)new_value;
							}
						}
						else // new meta info has not occurred before, create stats for it:
						{
							MetaStatsValue* m;
							if(new_value.valueType() == DataValue::INT_VALUE || new_value.valueType() == DataValue::DOUBLE_VALUE)
							{
								DoubleReal val = (DoubleReal)new_value;
								m = new MetaStatsValue(1,val,val,val);
							}
							else m = new MetaStatsValue(1,1,0,0); // min > max (illegal) indicates value not numerical
							pair<String, MetaStatsValue*> p = make_pair(*it_meta, m);
							meta_stats.insert(p);
						}
					}
				}
			}
			if (divisor != 0) avg_intensity /= (DoubleReal)divisor;
			for(MetaIterator it = meta_stats.begin(); it != meta_stats.end(); it++)
			{
				MetaStatsValue* meta_stats_value = it->second;
				if(meta_stats_value->count != 0)
				{
					meta_stats_value->avg /= (DoubleReal)meta_stats_value->count;
				}
			}
		}
		
		else if(layer_data.type == LayerData::DT_FEATURE)
		{
			tableWidget->setRowCount(tableWidget->rowCount() + 2);
			item = new QTableWidgetItem();
			item->setText(QString("Charge"));
			tableWidget->setVerticalHeaderItem(1, item);
			item = new QTableWidgetItem();
			item->setText(QString("Quality"));
			tableWidget->setVerticalHeaderItem(2, item);
			
			DoubleReal min_charge = 0, max_charge = 0, avg_charge = 0, min_quality = 0, max_quality = 0, avg_quality = 0;
			
			if(!layer_data.features.empty())
			{
				min_charge = layer_data.features.begin()->getCharge();
				max_charge = layer_data.features.begin()->getCharge();
				avg_charge = 0;
			
				min_quality = layer_data.features.begin()->getOverallQuality();
				max_quality = layer_data.features.begin()->getOverallQuality();
				avg_quality = 0;
			}
			
			unsigned long divisor = 0;
			for(FeatureIterator it_feature = layer_data.features.begin(); it_feature != layer_data.features.end(); it_feature++)
			{
				if(it_feature->getCharge() < min_charge) min_charge = it_feature->getCharge();
				if(it_feature->getCharge() > max_charge) max_charge = it_feature->getCharge();
				if(it_feature->getOverallQuality() < min_quality) min_quality = it_feature->getOverallQuality();
				if(it_feature->getOverallQuality() > max_quality) max_quality = it_feature->getOverallQuality();
				avg_intensity += it_feature->getIntensity();
				avg_charge += it_feature->getCharge();
				avg_quality += it_feature->getOverallQuality();
				divisor++;
				it_feature->getKeys(new_meta_keys);
				for(vector<String>::iterator it_meta = new_meta_keys.begin(); it_meta != new_meta_keys.end(); it_meta++)
				{
					new_value = it_feature->getMetaValue(*it_meta);
					MetaIterator it = meta_stats.find(*it_meta);
					if(it != meta_stats.end())
					{
						MetaStatsValue* meta_stats_value = it->second;
						meta_stats_value->count++;
						if(new_value.valueType() == DataValue::INT_VALUE || new_value.valueType() == DataValue::DOUBLE_VALUE)
						{
							if((DoubleReal)new_value < meta_stats_value->min) meta_stats_value->min = (DoubleReal)new_value;
							if((DoubleReal)new_value > meta_stats_value->max) meta_stats_value->max = (DoubleReal)new_value;
							meta_stats_value->avg += (DoubleReal)new_value;
						}
					}
					else // new meta info has not occurred before, create stats for it:
					{
						MetaStatsValue* m;
						if(new_value.valueType() == DataValue::INT_VALUE || new_value.valueType() == DataValue::DOUBLE_VALUE)
						{
							DoubleReal val = (DoubleReal)new_value;
							m = new MetaStatsValue(1,val,val,val);
						}
						else m = new MetaStatsValue(1,1,0,0); // min > max (illegal) indicates value not numerical
						pair<String, MetaStatsValue*> p = make_pair(*it_meta, m);
						meta_stats.insert(p);
					}
				}
			}
			if(divisor != 0)
			{
				avg_intensity /= (DoubleReal)divisor;
				avg_charge /= (DoubleReal)divisor;
				avg_quality /= (DoubleReal)divisor;
			}
			MetaStatsValue* meta_stats_value;
			for(MetaIterator it = meta_stats.begin(); it != meta_stats.end(); it++)
			{
				meta_stats_value = it->second;
				if(meta_stats_value->count != 0)
				{
					meta_stats_value->avg /= (DoubleReal)meta_stats_value->count;
				}
			}
			
			item = new QTableWidgetItem();
			item->setText(QString::number(min_charge,'f',3));
			tableWidget->setItem(1,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_charge,'f',3));
			tableWidget->setItem(1,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_charge,'f',3));
			tableWidget->setItem(1,3,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(min_quality,'f',3));
			tableWidget->setItem(2,1,item);
					
			item = new QTableWidgetItem();
			item->setText(QString::number(max_quality,'f',3));
			tableWidget->setItem(2,2,item);
			
			item = new QTableWidgetItem();
			item->setText(QString::number(avg_quality,'f',3));
			tableWidget->setItem(2,3,item);
			
		}
		
		item = new QTableWidgetItem();
		item->setText(QString::number(min_intensity,'f',3));
		tableWidget->setItem(0,1,item);
		
		item = new QTableWidgetItem();
		item->setText(QString::number(max_intensity,'f',3));
		tableWidget->setItem(0,2,item);
		
		item = new QTableWidgetItem();
		item->setText(QString::number(avg_intensity,'f',3));
		tableWidget->setItem(0,3,item);
		
		String name;
		MetaStatsValue* meta_stats_value;
		for(MetaIterator it = meta_stats.begin(); it != meta_stats.end(); it++)
		{
			tableWidget->setRowCount(tableWidget->rowCount()+1);
			name = it->first;
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
			else // min > max --> meta value was not numerical
			{
				item = new QTableWidgetItem();
				item->setText("-");
				tableWidget->setItem(tableWidget->rowCount()-1, 1, item);
				tableWidget->setItem(tableWidget->rowCount()-1, 2, item);
				tableWidget->setItem(tableWidget->rowCount()-1, 3, item);
			}
		}
	}
} // namespace
