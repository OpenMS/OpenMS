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
// $Id: Spectrum2DWidgetPDP.C,v 1.4 2006/05/30 13:48:19 cfriedle Exp $
// $Author: cfriedle $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

//OpenMS
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DWidgetPDP.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>

// Qt
#include <qlayout.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qcombobox.h>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{
	
		Spectrum2DWidgetPDP::Spectrum2DWidgetPDP( Spectrum2DWidget* manager, QWidget* parent, const char* name, WFlags f)
			: PreferencesDialogPage(manager,parent,name,f)
		{
			help_ = "This is the preferences dialog of 2D spectrum !"
							"<br>";
		
			QGridLayout* grid;
		
			grid = new QGridLayout(this,1,1);
			
			canvas_ = manager->client("Canvas", this);
			grid->addMultiCellWidget(canvas_, 0,1,0,1);
			
			QGroupBox* box = new QGroupBox(2,Qt::Horizontal,"Mapping",this);
			QLabel* label;
			label = new QLabel("Map m/z to: ",box);
			axis_mapping_ = new QComboBox(false, box, "read-only combobox");
			axis_mapping_->insertItem("X-Axis");
			axis_mapping_->insertItem("Y-Axis");  
			grid->addWidget(box,2,0);
		
			box = new QGroupBox(2,Qt::Horizontal,"Axis orientation (from lower left corner)",this);
			label = new QLabel("Select x axis orientation: ",box);
			x_axis_orientation_ = new QComboBox(false, box, "read-only combobox");
			x_axis_orientation_->insertItem("Ascending");
			x_axis_orientation_->insertItem("Descending");
		
			label = new QLabel("Select y axis orientation: ",box);
			y_axis_orientation_ = new QComboBox(false, box, "read-only combobox");
			y_axis_orientation_->insertItem("Ascending");
			y_axis_orientation_->insertItem("Descending");
			grid->addWidget(box,2,1);
		
			load();
		}
		
		Spectrum2DWidgetPDP::~Spectrum2DWidgetPDP()
		{
			
		}
		
		void Spectrum2DWidgetPDP::load()
		{
			Spectrum2DWidget* w = dynamic_cast<Spectrum2DWidget*>(manager_);
		  (w->canvas()->getMappingInfo().isMzToXAxis()) ? axis_mapping_->setCurrentText("X-Axis") : axis_mapping_->setCurrentText("Y-Axis");
		  (w->canvas()->getMappingInfo().isXAxisAsc())? x_axis_orientation_->setCurrentText("Ascending"): x_axis_orientation_->setCurrentText("Descending");
		  (w->canvas()->getMappingInfo().isYAxisAsc())? y_axis_orientation_->setCurrentText("Ascending"): y_axis_orientation_->setCurrentText("Descending");
		}
		
		void Spectrum2DWidgetPDP::save()
		{
			Spectrum2DWidget* w = dynamic_cast<Spectrum2DWidget*>(manager_);
			(axis_mapping_->currentText()=="X-Axis") ? w->switchAxis(false) : w->switchAxis(true);
			(x_axis_orientation_->currentText()=="Ascending") ? w->setMirroredXAxis(false) : w->setMirroredXAxis(true);
			(y_axis_orientation_->currentText()=="Ascending") ? w->setMirroredYAxis(false) : w->setMirroredYAxis(true);
		
			canvas_->save();
		}

	} // namespace Internal

} //namespace


