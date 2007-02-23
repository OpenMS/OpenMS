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

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM1DCANVASPDP_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM1DCANVASPDP_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>

namespace OpenMS
{
	class ColorSelector;
	class Spectrum1DCanvas;
	
	namespace Internal
	{
	
		///Preferences dialog page of a Spectrum1DCanvas widget (internal use only)	
		class Spectrum1DCanvasPDP: public PreferencesDialogPage
		{
			Q_OBJECT
			
			public:
				/// Constructor
				Spectrum1DCanvasPDP( Spectrum1DCanvas* manager, QWidget* parent = 0);
				/// Destructor
				virtual ~Spectrum1DCanvasPDP();
				// Docu in base class
				virtual void load();
				// Docu in base class
				virtual void save();
			
			protected:
				ColorSelector* peak_color_;
				ColorSelector* high_color_; // color of highlighted peak
				ColorSelector* icon_color_;
				ColorSelector* back_color_; // color of background
		};
		
	} //namespace Internal
	
} //namespace OpenMS

#endif //OPENMS_VISUAL_DIALOGS_SPECTRUM1DCANVASPDP_H

