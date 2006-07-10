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
// $Id: PreferencesManager.h,v 1.7 2006/03/03 17:55:35 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_PREFERENCESMANAGER_H
#define OPENMS_VISUAL_PREFERENCESMANAGER_H

// STL
#include <string>
#include <vector>

// OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialogPage.h>


namespace OpenMS
{
	class PreferencesDialog;
	class String;
	/**
		@brief Base class for all classes that use a Param to store preferences and want to display
		       dialog pages that allow editing the preferences.
		
		
		
		@ingroup Visual
	*/
	class PreferencesManager
	{
		public:
			///constructor
			PreferencesManager();
			///destructor
			virtual ~PreferencesManager();

			/// returns the parent PreferencesManager
			PreferencesManager* getParent();
			/// sets the parent PreferencesManager
			void setParent(PreferencesManager* parent);

			/// adds a child PreferencesManager and associated name (sets this PreferencesManger as the child's parent)
			void addClient(PreferencesManager* client, const std::string& name, bool isIncluded=false);
			/// removes a child PreferencesManager
			void removeClient(PreferencesManager* client, bool isIncluded=false);

			/// changes the name of a child PreferencesManager
			void setClientName(PreferencesManager* client, const std::string& name, bool isIncluded=false);

			/**
				@brief Creates and shows the preferences dialog.
				
				For return values see QDialog::exec() return values.
			*/
			SignedInt showPreferencesDialog();

			/**
				@brief Creates a PreferencesDialogPage with the preferences settings on it.
				
			 @param parent parent widget to be used for the background widget
			*/
			virtual PreferencesDialogPage* createPreferences(QWidget* parent) = 0;

			/**
				@brief Return the preferences dialog page of a client.
				
				It looks up a client with name @p name form incl_clients_. 
				Then creates its PreferencesDialogPage page with @p parent as parent.
				
				@return a pointer to the PreferencesDialogPage, or 0 if there is no client with this name.
			*/
			PreferencesDialogPage* client(std::string name,QWidget* parent);
			
			///if this is true, this page wants to be shown first in the main dialog
			bool isActive();
			///sets if the page is shown first in the main dialog
			void setActive(bool);

			/**
			  @brief Sets the preference @p name to @p value.
				
				This method is mainly used to set preferences from the preferences dialog.
			*/
			void setPref(const String& name, const String& value);
			/**
				@brief Sets the preference @p name to @p value.
				
				This method is mainly used to set preferences from the preferences dialog.
			*/
			void setPref(const String& name, SignedInt value);
			///returns the preference entry @p name
			const DataValue& getPref(const String& name) const;
			///returns the preference entry @p name as SignedInt
			SignedInt getPrefAsInt(const String& name) const;
			///returns the preference entry @p name as String
			String getPrefAsString(const String& name) const;
			///removes the preference entry @p name
			void removePref(const String& name);
		protected:

			/// calls createPreferences(QWidget* parent) and adds the widget into the stack with name <i>name</i>
			virtual void createPreferences_(PreferencesDialog* dialog, const std::string& name);

			/// parent preferences manager (that manages the preferences dialog)
			PreferencesManager* parent_;
			/// list of all client preferences managers (that are to be managed)
			std::map<PreferencesManager*,std::string> clients_;
			/// list of all client preferences managers (that are included into this page)
			std::map<PreferencesManager*,std::string> incl_clients_;
			/// here the actual preferences are stored
			Param prefs_;
			/// if this is true, this page wants to be shown first in the main dialog
		  bool is_active_;

	};
} // namespace OpenMS


#endif

