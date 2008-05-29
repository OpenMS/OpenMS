// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/TOPPViewOpenDialogTemplate.h>

namespace OpenMS 
{
	class Param;
	class String;
	/**
		@brief Dataset opening options for TOPPView
		
		@ingroup TOPPView_elements
	*/
	class TOPPViewOpenDialog
		: public QDialog,
  		public Ui::TOPPViewOpenDialogTemplate
	{
		Q_OBJECT
		
		public:			
			/// Constructor
			TOPPViewOpenDialog(const String& data_name, Param& preferences, QWidget* parent = 0 );
			/// Destructor
			virtual ~TOPPViewOpenDialog();
			
			/// Returns true, if 2D mode is to be used for maps
			bool viewMapAs2D() const;
			/// Returns of the low intensity peaks should be hidden
			bool isCutoffEnabled() const;
			/// Returns true, if the data should be opened in a new window
			bool openAsNewWindow() const;
			
		protected:
			/// Preferences
			Param& prefs_;
	};

}
#endif // OPENMS_VISUAL_DIALOGS_TOPPVIEWOPENDIALOG_H

