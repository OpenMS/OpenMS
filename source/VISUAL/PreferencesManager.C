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


#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/VISUAL/DIALOGS/PreferencesDialog.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <sstream>
#include <iostream>

#include <qdialog.h>
#include <qpushbutton.h>
#include <qlayout.h>

using namespace std;

namespace OpenMS
{

	PreferencesManager::PreferencesManager()
	:parent_(0), 
	 clients_(),
	 is_active_(0)
	{  
	}
	
	PreferencesManager::~PreferencesManager()
	{
		if (parent_!=0)
		{
			parent_->removeClient(this);
		}
	}
	
	PreferencesManager* PreferencesManager::getParent()
	{
			return parent_;
	}
	
	void PreferencesManager::setParent(PreferencesManager* parent)
	{	
		parent_ = parent;
	}
	
	void PreferencesManager::addClient(PreferencesManager* client, const string& name, bool isIncluded)
	{ 
			if (isIncluded){
			incl_clients_[client] = name;
			}else{
			clients_[client] = name;
			}
		client->setParent(this);
	}
	
	
	void PreferencesManager::removeClient(PreferencesManager* client, bool isIncluded)
	{ 	
		if (isIncluded){
			incl_clients_.erase(client);
		}else{
			clients_.erase(client);
		}
	}
	
	void PreferencesManager::setClientName(PreferencesManager* client, const std::string& name, bool isIncluded)
	{
		if (isIncluded){
			incl_clients_[client] = name;
		}else{
			clients_[client] = name;
		} 
	}
	
	SignedInt PreferencesManager::showPreferencesDialog()
	{
		if (parent_!=0)
		{
			return parent_->showPreferencesDialog();
		}
	
		//dialog
		PreferencesDialog dialog;
		
		createPreferences_(&dialog,"Default Preferences");
		
		
		if (dialog.exec())
		{
			return 1;
		}
		else
		{
			return 0;
		}
		
	}
	
	
	void PreferencesManager::createPreferences_(PreferencesDialog* dialog, const string& name)
	{
		dialog->addPage(name,createPreferences(dialog),this,getParent());
		//client preferences
		for (map<PreferencesManager*,string>::iterator it = clients_.begin();it!=clients_.end();++it)
		{
			it->first->createPreferences_(dialog,it->second);
		}
	}
	
	PreferencesDialogPage* PreferencesManager::client(std::string name, QWidget* parent)
	{
		for (map<PreferencesManager*,string>::iterator it = incl_clients_.begin(); it!=incl_clients_.end();++it)
		{
			if (it->second == name)
				return it->first->createPreferences(parent);
		}
		return 0;
	}
	
	bool PreferencesManager::isActive()
	{
		return is_active_;
	}
	
	void PreferencesManager::setActive(bool b)
	{
		is_active_ = b;
	}
	
	void PreferencesManager::setPref(const String& name,const String& value)
	{
		prefs_.setValue(name, value);
	}
	
	void PreferencesManager::setPref(const String& name,SignedInt value)
	{
		prefs_.setValue(name, value);
	}
	
	const DataValue& PreferencesManager::getPref(const String& name) const
	{
		return prefs_.getValue(name);
	}
	
	void PreferencesManager::removePref(const String& name)
	{
		return prefs_.remove(name);
	}
	
	///returns the preference entry @p name as SignedInt
	SignedInt PreferencesManager::getPrefAsInt(const String& name) const
	{
		//cout << name << ": " << (string)prefs_.getValue(name) << endl;
		return (SignedInt)(prefs_.getValue(name));
	}
	
	///returns the preference entry @p name as String
	String PreferencesManager::getPrefAsString(const String& name) const
	{
		return prefs_.getValue(name).toString();
	}

} //namespace OpenMS
