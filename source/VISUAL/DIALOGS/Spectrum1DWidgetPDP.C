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

#include <OpenMS/VISUAL/DIALOGS/Spectrum1DWidgetPDP.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>

// Qt
#include <QtGui/QLayout>
#include <QtGui/QLabel>
#include <QtGui/QGroupBox>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{
		
		Spectrum1DWidgetPDP::Spectrum1DWidgetPDP( Spectrum1DWidget* manager, QWidget* parent)
			: PreferencesDialogPage(manager,parent)
		{
			help_ = "This is the preferences dialog of a displayed spectrum!";
		
			QGridLayout* grid = new QGridLayout(this);
			
			//Specetrum1DCanvas settings
			colors_ = manager->client("Canvas", this);
			grid->addWidget(colors_,0,0);
			
			//mapping
			QGroupBox* box = addBox(grid,1,0,"Mapping");
			axis_mapping_ = new QComboBox( box);
			axis_mapping_->insertItem(0,"X-Axis");
			axis_mapping_->insertItem(1,"Y-Axis");  
			addWidget(box->layout(),0,"Map m/z to:",axis_mapping_);
			
			finish(grid);
			
			load();
		}
		
		Spectrum1DWidgetPDP::~Spectrum1DWidgetPDP()
		{
			
		}
		
		void Spectrum1DWidgetPDP::load()
		{
			Spectrum1DWidget* w = dynamic_cast<Spectrum1DWidget*>(manager_);
			
		  if (w->canvas()->isMzToXAxis()) 
		  {
		  	axis_mapping_->setCurrentIndex(0);
			}
			else
			{
				axis_mapping_->setCurrentIndex(1);
			}
		}
		
		void Spectrum1DWidgetPDP::save()
		{
			Spectrum1DWidget* w = dynamic_cast<Spectrum1DWidget*>(manager_);
			
			(axis_mapping_->currentText()=="X-Axis") ? w->mzToXAxis(true) : w->mzToXAxis(false);

			colors_->save();
		}

	} // namespace Internal

} //namespace


