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


#ifndef OPENMS_VISUAL_DIALOGS_PREFERENCESDIALOG_H
#define OPENMS_VISUAL_DIALOGS_PREFERENCESDIALOG_H

// QT
#include <QtGui/QDialog>

// STL
#include <string>
#include <vector>

namespace OpenMS
{
	class PreferencesDialogPage;
	class PreferencesManager;
	class ListStack;
	
	/**
		@brief Main dialog for the PreferencesManager classes.	
		
		@ingroup Dialogs
	*/
	class PreferencesDialog: public QDialog
	{
		Q_OBJECT

		public:
			///constructor
			PreferencesDialog();
			///destructor
			virtual ~PreferencesDialog();


			/// add a new PreferencesDailogPage
			void addPage(std::string name, PreferencesDialogPage* page, PreferencesManager* creator, bool highlight, PreferencesManager* parent = 0);

		protected slots:
			/// ok button pressed
			void ok_();
			/// cancel button pressed
			void cancel_();
			/// apply button pressed
			void apply_();
			/// help button pressed
			void help_();
		
		protected:
			/// the ListStack in which the PreferencesDialogPages are displayed
			ListStack* stack_;
			/// a vector with pointers to all pages
			std::vector<PreferencesDialogPage*> pages_;
	};
} // namespace OpenMS


#endif

