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
// $Id: PreferencesDialogPage.h,v 1.4 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_PREFERENCESDIALOGPAGE_H
#define OPENMS_VISUAL_DIALOGS_PREFERENCESDIALOGPAGE_H

// QT
#include <qwidget.h>

// STL
#include <string>

// OpenMS
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	class PreferencesManager;
	
	/**
		@brief Base class for all PreferencesManager dialog pages.
		
		
		
		@ingroup Dialogs
	*/
	class PreferencesDialogPage: public QWidget
	{
		Q_OBJECT

		public:
			///constructor
			PreferencesDialogPage( PreferencesManager* manager, QWidget* parent = 0, const char* name = "PreferencesDialogPage", WFlags f = 0);
			///destructor
			virtual ~PreferencesDialogPage();

			/// returns the help text
			const std::string& getHelpText() const;
			///sets the help text
			void setHelpText(const std::string& text);

			/// load values from the PreferencesManager's Param
			virtual void load();
			/// write changes in Dialog to the PreferencesManager's Param
			virtual void save();



		protected:
			///Pointer to the PreferencesManager
			PreferencesManager* manager_;
			///stores the help text
			std::string help_;

	};
} // namespace OpenMS


#endif

