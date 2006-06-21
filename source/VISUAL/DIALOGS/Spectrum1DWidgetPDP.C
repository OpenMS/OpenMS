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
// $Id: Spectrum1DWidgetPDP.C,v 1.3 2006/03/28 12:53:15 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/VISUAL/DIALOGS/Spectrum1DWidgetPDP.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

// Qt
#include <qlayout.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qcheckbox.h>
#include <qcombobox.h>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{
		
		Spectrum1DWidgetPDP::Spectrum1DWidgetPDP( Spectrum1DWidget* manager, QWidget* parent, const char* name, WFlags f)
			: PreferencesDialogPage(manager,parent,name,f)
		{
			help_ = "This is the preferences dialog of 1D spectrum!"
							"<br>";
		
			QGridLayout* grid;
		
			//1D View Tab
			grid = new QGridLayout(this,3,2);
		
			grid->setMargin(6);
			grid->setSpacing(4);	
			
			colors_ = manager->client("Canvas", this);
			
			grid->addWidget(colors_,0,0);
		
			QGroupBox* box = new QGroupBox(2,Qt::Horizontal,"Mapping",this);
			QLabel* label = new QLabel("Map m/z to: ",box);
			axis_mapping_ = new QComboBox(false, box, "read-only combobox");
			axis_mapping_->insertItem("X-Axis");
			axis_mapping_->insertItem("Y-Axis");  
			grid->addWidget(box,1,0);
		
			box = new QGroupBox(2,Qt::Horizontal,"Axis orientation (from lower left corner)",this);
			label = new QLabel("Select x axis orientation: ",box);
			x_axis_orientation_ = new QComboBox(false, box, "read-only combobox");
			x_axis_orientation_->insertItem("Ascending");
			x_axis_orientation_->insertItem("Descending");
		
			label = new QLabel("Select y axis orientation: ",box);
			y_axis_orientation_ = new QComboBox(false, box, "read-only combobox");
			y_axis_orientation_->insertItem("Ascending");
			y_axis_orientation_->insertItem("Descending");
			grid->addWidget(box,1,1);
		
			box = new QGroupBox( 2,Qt::Vertical,"Intensity scale",this);
			log_check_box_ = new QCheckBox("Logarithmic", box, "log CheckBox" );
			rel_check_box_ = new QCheckBox("Relative(%)",box, "log CheckBox" );
			grid->addWidget(box,0,1);
		
			box = new QGroupBox( 2,Qt::Vertical, "X-Axis",this);
			x_page_= manager->client("X-Axis",box);
			grid->addWidget(box,2,0);
		
			box = new QGroupBox( 2,Qt::Vertical,"Y-Axis",this);
			y_page_ = manager->client("Y-Axis",box);
			grid->addWidget(box,2,1);
		
			load();
		}
		
		Spectrum1DWidgetPDP::~Spectrum1DWidgetPDP()
		{
			
		}
		
		void Spectrum1DWidgetPDP::load()
		{
			Spectrum1DWidget* w = dynamic_cast<Spectrum1DWidget*>(manager_);
			
		  (w->canvas()->getMappingInfo().isMzToXAxis()) ? axis_mapping_->setCurrentText("X-Axis") : axis_mapping_->setCurrentText("Y-Axis");
		  (w->canvas()->getMappingInfo().isXAxisAsc())? x_axis_orientation_->setCurrentText("Ascending"): x_axis_orientation_->setCurrentText("Descending");
		  (w->canvas()->getMappingInfo().isYAxisAsc())? y_axis_orientation_->setCurrentText("Ascending"): y_axis_orientation_->setCurrentText("Descending");
		
			log_check_box_->setChecked(w->isLogIntensity());
			rel_check_box_->setChecked(!w->canvas()->isAbsoluteIntensity());
		}
		
		void Spectrum1DWidgetPDP::save()
		{
			Spectrum1DWidget* w = dynamic_cast<Spectrum1DWidget*>(manager_);
			
			(axis_mapping_->currentText()=="X-Axis") ? w->switchAxis(false) : w->switchAxis(true);
			(x_axis_orientation_->currentText()=="Ascending") ? w->setMirroredXAxis(false) : w->setMirroredXAxis(true);
			(y_axis_orientation_->currentText()=="Ascending") ? w->setMirroredYAxis(false) : w->setMirroredYAxis(true);
					
			if (log_check_box_->isChecked())
			{
				w->setIntensityModificationLog();
			}
			else
			{
				w->setIntensityModificationNone();
			}
			if (rel_check_box_->isChecked())
			{
				w->intensityAxisRelative();
			}
			else
			{
				w->intensityAxisAbsolute();
			}
			x_page_->save();
			y_page_->save();
		
			colors_->save();
		}

	} // namespace Internal

} //namespace


