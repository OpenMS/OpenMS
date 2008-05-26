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
#include <OpenMS/FORMAT/DB/DBAdapter.h>
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPViewOpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/DBOpenDialog.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>


// QT includes
#include <QtGui/QMessageBox>
#include <QtGui/QRadioButton>
#include <QtGui/QFileDialog>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QInputDialog>
#include <QtGui/QComboBox>
#include <QtGui/QButtonGroup>

// STL includes
#include <iostream>

using namespace std;

namespace OpenMS
{

	TOPPViewOpenDialog::TOPPViewOpenDialog(const String& data, bool file, Param& preferences, QWidget * parent)
		: QDialog(parent),
			data_(data),
			is_file_(file),
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

		//init force file type
		FileHandler fh;
		force_->insertItem(0,"Detect automatically",0);
		for (int i=1; i< FileHandler::SIZE_OF_TYPE; ++i)
		{
			FileHandler::Type type = (FileHandler::Type)i;
			if (type!=FileHandler::PARAM && type!=FileHandler::IDXML && type!=FileHandler::CONSENSUSXML)
			{
				force_->insertItem(force_->count(),fh.typeToName(type).c_str(),i);
			}
		}
			
		//connect meta data browsing
		connect(metadata_,SIGNAL(pressed()),this,SLOT(showMetaData_()));

		//do file/DB specific stuff
		if (is_file_)
		{
			setWindowTitle((String("Open file: ") + File::basename(data_)).toQString());
		}
		else
		{
			setWindowTitle((String("Open database entry: ") + data_.toQString()).toQString());
			force_->setEnabled(false);
		}
	}
	
	TOPPViewOpenDialog::~TOPPViewOpenDialog()
	{
	}

	void TOPPViewOpenDialog::showMetaData_()
	{
		MSExperiment<> exp;
		
		//try to open file or database entry
		try
		{
			if (is_file_)
			{
				FileHandler fh;
				fh.getOptions().setMetadataOnly(true);
				fh.loadExperiment(data_,exp,forcedFileType());
			}
			else
			{
				DBConnection con;
				con.connect(prefs_.getValue("preferences:db:name"), prefs_.getValue("preferences:db:login"),prefs_.getValue("DBPassword"),prefs_.getValue("preferences:db:host"),(Int)prefs_.getValue("preferences:db:port"));
				DBAdapter db(con);
				db.getOptions().setMetadataOnly(true);
				db.loadExperiment(data_.toInt(), exp);
			}
		}
		catch (Exception::Base& e)
		{
			QMessageBox::critical(this,"Error",(String("Error while reading data: ")+e.what()).c_str());
      return;
		}
		MSMetaDataExplorer dlg(false, this);
		dlg.setWindowTitle("Meta data");			
		dlg.visualize(exp);
 	 	dlg.exec();
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

  FileHandler::Type TOPPViewOpenDialog::forcedFileType() const
  {
  	return (FileHandler::Type) (force_->currentIndex());
  }

}



