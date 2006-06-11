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
// $Id: Spectrum1DWidgetPDP.h,v 1.4 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM1DWIDGETPDP_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM1DWIDGETPDP_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>

class QCheckBox;
class QComboBox;

namespace OpenMS
{
	class Spectrum1DWidget;

	namespace Internal
	{

		///Preferences dialog page of a Spectrum1DWidget (internal use only)	
		class Spectrum1DWidgetPDP: public PreferencesDialogPage
		{
			Q_OBJECT
	
			public:
				Spectrum1DWidgetPDP( Spectrum1DWidget* manager, QWidget* parent = 0, const char* name = "Spectrum1DWidgetPDP", WFlags f = 0);
				virtual ~Spectrum1DWidgetPDP();
				virtual void load();
				virtual void save();
			
			protected:            
				PreferencesDialogPage* colors_;
				QComboBox* axis_mapping_;
				QComboBox* x_axis_orientation_;
				QComboBox* y_axis_orientation_;
				QCheckBox* log_check_box_;
				QCheckBox* rel_check_box_;
				PreferencesDialogPage* x_page_;
				PreferencesDialogPage* y_page_;
		};
	
	} //namespace Internal
	
} //namespace OpenMS

#endif //OPENMS_VISUAL_DIALOGS_SPECTRUM1DWIDGETPDP_H

