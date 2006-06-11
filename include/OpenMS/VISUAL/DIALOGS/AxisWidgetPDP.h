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
// $Id: AxisWidgetPDP.h,v 1.4 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_AXISWIDGETPDP_H
#define OPENMS_VISUAL_DIALOGS_AXISWIDGETPDP_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>

class QCheckBox;
class QListBox;

namespace OpenMS
{
	class AxisWidget;
	
	namespace Internal
	{
		///Preferences dialog page of an AxisWidget (internal use only)
		class AxisWidgetPDP: public PreferencesDialogPage
		{
			Q_OBJECT
			
			public:
				AxisWidgetPDP( AxisWidget* manager, QWidget* parent = 0, const char* name = "AxisWidtetPDP_", WFlags f = 0);
				virtual ~AxisWidgetPDP();
				virtual void load();
				virtual void save();
			protected:
				QCheckBox* show_legend_;
				QListBox* unit_box_;
		};
		
	} //namespace Internal

} //namespace OpenMS

#endif //OPENMS_VISUAL_DIALOGS_AXISWIDGETPDP_H

