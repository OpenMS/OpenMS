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

// OpenMS includes
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/DATASTRUCTURES/Param.h>


// QT includes
#include <QtGui/QButtonGroup>

// STL includes
#include <iostream>

using namespace std;

namespace OpenMS
{

	TOPPViewOpenDialog::TOPPViewOpenDialog(const String& data_name, Param& preferences, QWidget * parent)
		: QDialog(parent),
			prefs_(preferences)
	{
		setupUi(this);
		
		//init map view
		QButtonGroup* button_group = new QButtonGroup(this);
		button_group->addButton(d2_);
		button_group->addButton(d3_);
		if ((String)(prefs_.getValue("preferences:default_map_view"))=="3d")
		{
			d3_->setChecked(true);
		}
		else
		{
			d2_->setChecked(true);
		}

		//init intensity cutoff
		button_group = new QButtonGroup(this);
		button_group->addButton(cutoff_);
		button_group->addButton(nocutoff_);
		if ((String)(prefs_.getValue("preferences:intensity_cutoff"))=="off")
		{
			nocutoff_->setChecked(true);
		}
		else
		{
			cutoff_->setChecked(true);
		}
		
		//init open as
		button_group = new QButtonGroup(this);
		button_group->addButton(window_);
		button_group->addButton(layer_);
		window_->setChecked(true);
			
		//do file/DB specific stuff
		setWindowTitle((String("Open data options for ") + data_name).toQString());
	}
	
	TOPPViewOpenDialog::~TOPPViewOpenDialog()
	{
	}

	bool TOPPViewOpenDialog::viewMapAs2D() const
	{
		if (d2_->isChecked()) return true;
		return false;
	}
	
	bool TOPPViewOpenDialog::isCutoffEnabled() const
	{
		if (cutoff_->isChecked()) return true;
		return false;
	}
	
	bool TOPPViewOpenDialog::openAsNewWindow() const
	{
		if (window_->isChecked()) return true;
		return false;	
	}
}



