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
#include <qlayout.h>
#include <qlabel.h>
#include <qgroupbox.h>


using namespace std;

namespace OpenMS
{

	namespace Internal
	{
	
		Spectrum1DCanvasPDP::Spectrum1DCanvasPDP(Spectrum1DCanvas* manager, QWidget* parent, const char* name, WFlags f)
			: PreferencesDialogPage(manager,parent,name,f)
		{
			help_ = "This is the preferences dialog of 1D spectrum!"
							"<br>";
		
			QGridLayout* grid;
		
			//1D View Tab
			grid = new QGridLayout(this, 1, 1);
			
			QGroupBox* box = new QGroupBox(2, Qt::Horizontal,"Colors",this);
			QLabel* label;
			label = new QLabel("Peak color: ",box);
			peak_color_ = new ColorSelector(box);
			label = new QLabel("Icon color: ",box);
			icon_color_ = new ColorSelector(box);
			label = new QLabel("Highlighted peak color: ",box);
			high_color_ = new ColorSelector(box);
			label = new QLabel("Background color: ",box);
			back_color_ = new ColorSelector(box);
			
			grid->addWidget(box, 0, 0);
		
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
			
			w->setPref("Preferences:1D:PeakColor",peak_color_->getColor().name().ascii());
			w->setPref("Preferences:1D:HighColor",high_color_->getColor().name().ascii());
			w->setPref("Preferences:1D:IconColor",icon_color_->getColor().name().ascii());
			w->setPref("Preferences:1D:BackgroundColor",back_color_->getColor().name().ascii());
		
			w->repaintAll();
		}

	} // namespace Internal

} //namespace


