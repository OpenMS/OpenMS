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

#include <OpenMS/VISUAL/DIALOGS/Spectrum1DCanvasPDP.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/ColorSelector.h>

// Qt
#include <QtGui/QLayout>
#include <QtGui/QLabel>
#include <QtGui/QGroupBox>
#include <QtGui/QGridLayout>


using namespace std;

namespace OpenMS
{

	namespace Internal
	{
	
		Spectrum1DCanvasPDP::Spectrum1DCanvasPDP(Spectrum1DCanvas* manager, QWidget* parent)
			: PreferencesDialogPage(manager,parent)
		{
			help_ = "This is the preferences dialog of a displayed spectrum!";
		
			QGridLayout* grid = new QGridLayout(this);
			
			QGroupBox* box = addBox(grid,0,0,"Colors:");

			peak_color_ = new ColorSelector(box);
			addWidget(box->layout(),0,"Peak color:", peak_color_);

			icon_color_ = new ColorSelector(box);
			addWidget(box->layout(),1,"Icon color:", icon_color_);
			
			high_color_ = new ColorSelector(box);
			addWidget(box->layout(),2,"Highlighted peak color:", high_color_);
			
			back_color_ = new ColorSelector(box);
			addWidget(box->layout(),3,"Background color:", back_color_);
		
			load();
		}
		
		Spectrum1DCanvasPDP::~Spectrum1DCanvasPDP()
		{
			
		}
		
		void Spectrum1DCanvasPDP::load()
		{
			Spectrum1DCanvas* w = dynamic_cast<Spectrum1DCanvas*>(manager_);
			
			peak_color_->setColor(QColor(w->getPrefAsString("Preferences:1D:PeakColor").c_str()));
			icon_color_->setColor(QColor(w->getPrefAsString("Preferences:1D:IconColor").c_str()));
			high_color_->setColor(QColor(w->getPrefAsString("Preferences:1D:HighColor").c_str()));
			back_color_->setColor(QColor(w->getPrefAsString("Preferences:1D:BackgroundColor").c_str()));
		}
		
		void Spectrum1DCanvasPDP::save()
		{
			Spectrum1DCanvas* w = dynamic_cast<Spectrum1DCanvas*>(manager_);
			
			w->setPref("Preferences:1D:PeakColor",peak_color_->getColor().name().toAscii().data());
			w->setPref("Preferences:1D:HighColor",high_color_->getColor().name().toAscii().data());
			w->setPref("Preferences:1D:IconColor",icon_color_->getColor().name().toAscii().data());
			w->setPref("Preferences:1D:BackgroundColor",back_color_->getColor().name().toAscii().data());
		
			w->repaintAll();
		}

	} // namespace Internal

} //namespace


