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
// $Id: Spectrum3DWindow.h,v 1.8 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUM3DWINDOW_H
#define OPENMS_VISUAL_SPECTRUM3DWINDOW_H

#ifndef OPENMS_CONFIG_H
#include <OpenMS/config.h>
#endif
#include <OpenMS/VISUAL/SpectrumWindow.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>

namespace OpenMS
{
	class Spectrum3DWidget;
	
	/**
		@brief Window for 3D-visualization of map data
		
		
		
		@ingroup spectrum_widgets
	*/
	class Spectrum3DWindow : public SpectrumWindow
	{
		Q_OBJECT
		
		public:
			/// Constructor
			Spectrum3DWindow(QWidget* parent=0, const char* name="Spectrum3DWindow", WFlags f=0);
			/// Destructor
			virtual ~Spectrum3DWindow();
			///
			Spectrum3DWidget* widget();
		
			///
			virtual PreferencesDialogPage* createPreferences(QWidget* parent);

			///
		  virtual void showGoToDialog();

	protected:
		void createContextMenu_();
	};
} //namespace
	
#endif
