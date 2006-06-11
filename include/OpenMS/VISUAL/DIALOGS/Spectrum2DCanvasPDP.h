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
// $Id: Spectrum2DCanvasPDP.h,v 1.5 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
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
				Spectrum2DCanvasPDP( Spectrum2DCanvas* manager, QWidget* parent = 0, const char* name = "Spectrum2DCanvasPDP", WFlags f = 0);
				virtual ~Spectrum2DCanvasPDP();
				virtual void load();
				virtual void save();
	
			protected:
				QRadioButton* dot_mode_black_;
				QRadioButton* dot_mode_gradient_;
				MultiGradientSelector* dot_gradient_;
				MultiGradientSelector* surface_gradient_;
			  ColorSelector* background_color_;
			  QSpinBox* marching_squares_steps_;
			  QSpinBox* dot_interpolation_steps_;
			  QSpinBox* surface_interpolation_steps_;
		};
	
	} //namespace Internal
	
} //namespace OpenMS

#endif //OPENMS_VISUAL_DIALOGS_SPECTRUM2DCANVASPDP_H


