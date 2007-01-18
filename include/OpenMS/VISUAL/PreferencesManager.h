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

	  GUI classes that store preferences and create a preferences
	  dialog for them, are derived from PreferencesManager. 
	  The preferences have to be stored
	  in the private member PreferencesManager::prefs_, an object of type Param (see
	  FORMAT/Param.h). PreferencesManager classes are arranged in a tree structure,
	  where one PreferencesManager (the parent) can have several clients. Clients are
	  added to a parent along with a name, which is used later on in
	  the visualization of the preferences dialog. See getParent(),
	  setParent(parent), addClient(client, name), removeClient(client)
	  and setClientName(client, name). When the showPreferencesDialog()
	  method of any object in the tree is called, the root object
	  creates a PreferencesDialog and shows it. On the left of the
	  dialog a tree view that corresponds to the object tree is shown
	  and on the right the preferences of the selected entry of the
	  tree view are shown (see figure). All preferences pages on the
	  right side are objects of classes derived from
	  PreferencesDialogPage (abbreviated PDP). 
		
		\image html PreferencesManager.png
		
	  <b>Using the default mechanism for adding new Widgets</b>

	  In order to store Preferences of a GUI class and a create
	  preferences dialog for them, derive your class from PreferencesManager. If you
	  want to use the default preferences mechanism, just add the
	  object as client to a parent. Now your object will have its own
	  dialog page in the preferences dialog. However this page is still
	  empty! Next, you have to derive a class from PDP in order to fill
	  your object's preferences page with life (This new class may be a
	  private class of your PreferencesManager object, which is defined in the same
	  header/source). In the constructor of the new PDP class you
	  create all the widgets you want to show on the preferences page.
	  In oder to see the new PDP in the preferences dialog, you have to
	  reimplement the createPreferences(parent) method of you PreferencesManager class.
	  By default this method creates a PDP object and returns a pointer
	  to it. Now you have to create and return a pointer to a object of
	  your class derived from PDP. Now you can see the preferences page
	  but it does not interact with your object at all... The PDP
	  object has a pointer to its creator(or manager) which is an
	  argument to its constructor. Through this pointer, the load()
	  method reads the current settings from the manager and save()
	  writes the preferences back to the manager. Finally, you have to
	  reimplement the load() and save() methods of the PDP. Important:
	  This mechanism works only when there are public set/get methods
	  for all preferences, but they should be there anyway.

	  <b>Parents managing clients themselves</b>

	  In cases where you do not want a client to be shown as a child in
	  the tree structure, you can just integrate the client's
	  preferences page into the preferences page of the parent. To do
	  that you have to modify the instructions above in the following
	  manner:
	
	  <ul>
	    <li>Set the the parent property of the client, but do not add
	    the client as child to the parent = &gt; The preference page is
	    not shown in the dialog anymore.</li>
	
	    <li>Integrate the child PDP into the client PDP by calling
	    createPreferences() of the child in the createPreferences()
	    method of the client.</li>
	
	    <li>Call the load()/save() method of the client in the
	    load()/save() method of the parent</li>
	  </ul>

	  @note Preferences that are used very often or are time-consuming to
				  construct from prefs_ can be stored in a second member for faster
				  access. But they should nevertheless be stored in prefs_, as
				  prefs_ is used to persistently store the current state of
				  TOPPView.			
		
		
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

