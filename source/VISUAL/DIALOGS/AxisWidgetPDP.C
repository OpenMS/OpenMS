// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: AxisWidgetPDP.C,v 1.3 2006/03/28 12:53:15 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/VISUAL/DIALOGS/AxisWidgetPDP.h>
#include <OpenMS/VISUAL/AxisWidget.h>

#include <qcheckbox.h>
#include <qlayout.h>
#include <qlistbox.h>
#include <qgroupbox.h>
#include <qlabel.h>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{
		
		AxisWidgetPDP::AxisWidgetPDP( AxisWidget* manager, QWidget* parent, const char* name, WFlags f)
			:PreferencesDialogPage(manager,parent,name,f)
		{
			help_ = "This is the preferences dialog of the axis widget!"
							"<br>";
			QVBoxLayout* layout;
		  QVBoxLayout* box_layout;
		  QGroupBox* legend_box;
		 	QWidget* background = new QWidget(parent);
		                        
			layout = new QVBoxLayout( background, 11, 6, "AxisWidgetPDP layout"); 
			layout->setResizeMode(QLayout::Fixed);
		
			legend_box = new QGroupBox(0,Qt::Vertical,"Legend",background);
			legend_box->layout()->setSpacing( 6 );
			legend_box->layout()->setMargin( 11 );
			box_layout = new QVBoxLayout( legend_box->layout() );
			box_layout->setAlignment( Qt::AlignTop );
		
		  unit_box_ = 0;
		  const std::vector<std::string>& units = dynamic_cast<AxisWidget*>(manager_)->getExtendedLegend();
			if (units.size()>0){   // There is more then one unit available as legend
				QLabel* legend_label = new QLabel( legend_box, "legendLabel" );
				box_layout->addWidget( legend_label );
				legend_label->setText( tr( "Units:" ) );
		    unit_box_ = new QListBox( legend_box, "unitBox" );
				box_layout->addWidget( unit_box_ );
		    unit_box_->clear();
				for (std::vector<std::string>::const_iterator it = units.begin(); it!=units.end(); it++)
					unit_box_->insertItem((*it).c_str());
		  }
		  
			show_legend_ = new QCheckBox( legend_box, "showLegend" );
			box_layout->addWidget( show_legend_ );
			layout->addWidget( legend_box );
		
			show_legend_->setText( tr( "Legend visible on axis" ) );
			load();
		}
					
		AxisWidgetPDP::~AxisWidgetPDP()
		{
		}
		
		void AxisWidgetPDP::load()
		{
		  if (unit_box_ != 0){   // There is more then one unit available as legend
		 		unit_box_->setSelected(dynamic_cast<AxisWidget*>(manager_)->getUnit(),true);
			}
		  show_legend_->setChecked(dynamic_cast<AxisWidget*>(manager_)->isLegendShown());
		}
		
		void AxisWidgetPDP::save()
		{
			if (unit_box_ != 0){   // There is more then one unit available as legend
				dynamic_cast<AxisWidget*>(manager_)->setUnit(unit_box_->currentItem());
			}
			dynamic_cast<AxisWidget*>(manager_)->showLegend(show_legend_->isChecked());
		}

	} // namespace Internal

} //namespace


