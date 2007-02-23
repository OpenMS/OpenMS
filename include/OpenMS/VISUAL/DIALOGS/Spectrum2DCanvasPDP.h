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

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM2DCANVASPDP_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM2DCANVASPDP_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>

class QRadioButton;
class QSpinBox;

namespace OpenMS
{
	class MultiGradientSelector;
	class Spectrum2DCanvas;
	class ColorSelector;

	namespace Internal
	{
	
		///Preferences dialog page of a Spectrum2DCanvas (internal use only)	
		class Spectrum2DCanvasPDP: public PreferencesDialogPage
		{
			Q_OBJECT
			
			public:
				/// Constructor
				Spectrum2DCanvasPDP( Spectrum2DCanvas* manager, QWidget* parent = 0);
				///  Destructor
				virtual ~Spectrum2DCanvasPDP();
				// Docu in base class
				virtual void load();
				// Docu in base class
				virtual void save();
	
			protected:
				QRadioButton* dot_mode_black_;
				QRadioButton* dot_mode_gradient_;
				MultiGradientSelector* dot_gradient_;
				MultiGradientSelector* surface_gradient_;
			  ColorSelector* background_color_;
			  QSpinBox* marching_squares_steps_;
			  QSpinBox* interpolation_steps_;
			  QSpinBox* contour_steps_;
		};
	
	} //namespace Internal
	
} //namespace OpenMS

#endif //OPENMS_VISUAL_DIALOGS_SPECTRUM2DCANVASPDP_H


