// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <QtGui/QLayout>
#include <QtGui/QLabel>
#include <QtGui/QGroupBox>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{
	
		Spectrum2DWidgetPDP::Spectrum2DWidgetPDP( Spectrum2DWidget* manager, QWidget* parent)
			: PreferencesDialogPage(manager,parent)
		{
			help_ = "This is the preferences dialog of a 2D view of a map!";
		
			QGridLayout* grid = new QGridLayout(this);
			
			//Spectrum2DCanvas
			canvas_ = manager->client("Canvas", this);
			grid->addWidget(canvas_, 0,0,1,2);
			
			//mapping
			QGroupBox* box = addBox(grid,1,0,"Mapping");
			axis_mapping_ = new QComboBox( box);
			axis_mapping_->insertItem(0,"X-Axis");
			axis_mapping_->insertItem(1,"Y-Axis");
			addWidget(box->layout(),0,"Map m/z to:",axis_mapping_);
			finish(box->layout());			
			
			load();
		}
		
		Spectrum2DWidgetPDP::~Spectrum2DWidgetPDP()
		{
			
		}
		
		void Spectrum2DWidgetPDP::load()
		{
			Spectrum2DWidget* w = dynamic_cast<Spectrum2DWidget*>(manager_);
			if (w->canvas()->isMzToXAxis())
			{
				axis_mapping_->setCurrentIndex(0);
			}
			else
			{
				axis_mapping_->setCurrentIndex(1);
			}
		}
		
		void Spectrum2DWidgetPDP::save()
		{
			Spectrum2DWidget* w = dynamic_cast<Spectrum2DWidget*>(manager_);
			(axis_mapping_->currentText()=="X-Axis") ? w->mzToXAxis(true) : w->mzToXAxis(false);
		
			canvas_->save();
		}

	} // namespace Internal

} //namespace


