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
// $Id: Spectrum2DWidgetPDP.h,v 1.4 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM2DWIDGETPDP_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM2DWIDGETPDP_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>

class QComboBox;

namespace OpenMS
{
	class Spectrum2DWidget;

	namespace Internal
	{

		///Preferences dialog page of a Spectrum2DWidget (internal use only)	
		class Spectrum2DWidgetPDP: public PreferencesDialogPage
		{
			Q_OBJECT
			
			public:
				Spectrum2DWidgetPDP( Spectrum2DWidget* manager, QWidget* parent = 0, const char* name = "Spectrum2DWidgetPDP", WFlags f = 0);
				virtual ~Spectrum2DWidgetPDP();
				virtual void load();
				virtual void save();
			protected:
				PreferencesDialogPage* canvas_;
			  QComboBox* axis_mapping_;
				QComboBox* x_axis_orientation_;
	  		QComboBox* y_axis_orientation_;
		};
	
	} //namespace Internal
	
} //namespace OpenMS

#endif //OPENMS_VISUAL_DIALOGS_SPECTRUM2DWIDGETPDP_H

