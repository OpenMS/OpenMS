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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

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
			grid->setMargin(6);
			grid->setSpacing(4);	
			
			canvas_ = manager->client("Canvas", this);
			grid->addMultiCellWidget(canvas_, 0,1,0,1);
			
			QGroupBox* box = new QGroupBox(2,Qt::Horizontal,"Mapping",this);
			QLabel* label;
			label = new QLabel("Map m/z to: ",box);
			axis_mapping_ = new QComboBox(false, box, "read-only combobox");
			axis_mapping_->insertItem("X-Axis");
			axis_mapping_->insertItem("Y-Axis");  
			grid->addWidget(box,2,0);

			load();
		}
		
		Spectrum2DWidgetPDP::~Spectrum2DWidgetPDP()
		{
			
		}
		
		void Spectrum2DWidgetPDP::load()
		{
			Spectrum2DWidget* w = dynamic_cast<Spectrum2DWidget*>(manager_);
		  (w->canvas()->isMzToXAxis()) ? axis_mapping_->setCurrentText("X-Axis") : axis_mapping_->setCurrentText("Y-Axis");
		}
		
		void Spectrum2DWidgetPDP::save()
		{
			Spectrum2DWidget* w = dynamic_cast<Spectrum2DWidget*>(manager_);
			(axis_mapping_->currentText()=="X-Axis") ? w->mzToXAxis(true) : w->mzToXAxis(false);
		
			canvas_->save();
		}

	} // namespace Internal

} //namespace


