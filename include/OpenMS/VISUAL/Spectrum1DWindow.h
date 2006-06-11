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
// $Id: Spectrum1DWindow.h,v 1.15 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM1DWINDOW_H
#define OPENMS_VISUAL_SPECTRUM1DWINDOW_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/SpectrumWindow.h>

namespace OpenMS
{
	class Spectrum1DWidget;

	/**
		@brief Window for visualization of spectrum
		
		
		
		@ingroup spectrum_widgets
	*/
	class Spectrum1DWindow : public SpectrumWindow
	{
		Q_OBJECT
		public:
			Spectrum1DWindow(QWidget* parent=0, const char* name="Spectrum1DWindow", WFlags f=0);
			~Spectrum1DWindow();
			
			Spectrum1DWidget* widget();
			
			virtual void setDrawMode(QAction* a);
			virtual int getDrawMode();
			
			void switchAxis(bool b);
			void setMirroredXAxis(bool b);
			void setMirroredYAxis(bool b);
			///hands the preferences of the parent object down to the child
			virtual void setMainPreferences(const Param& prefs);

			///PreferencesManager
			virtual PreferencesDialogPage* createPreferences(QWidget* parent);

			bool getSnapToMax();
			void setSnapToMax(bool b);
			
		public slots:
			void setLoXHiXNoEmit(double,double);
      virtual void showGoToDialog();				
		protected:
			virtual void createContextMenu_();
		signals:
			void loXHiXChanged(double,double);	
			void changeLoXHiX(double,double);
	};
}
#endif

