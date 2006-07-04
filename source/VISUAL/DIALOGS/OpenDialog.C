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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/OpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/DBSpectrumSelectorDialog.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/DBConnection.h>

// QT includes
#include <qmessagebox.h>
#include <qbuttongroup.h>
#include <qradiobutton.h>
#include <qfiledialog.h>
#include <qlineedit.h>
#include <qpushbutton.h>
#include <qinputdialog.h>
#include <qcombobox.h>

// STL includes
#include <iostream>

using namespace std;
namespace OpenMS
{

	OpenDialog::OpenDialog( Param& preferences, QWidget * parent, const char * name, WFlags fl)
		: OpenDialogTemplate(parent,name,fl),
			prefs_(preferences)
	{
		if (!(bool(getPrefAsInt_("Preferences:DefaultMapView"))))
		{
			d3_radio->setChecked(true);
		}
		if ((String)(getPref_("Preferences:MapIntensityCutoff"))!="None")
		{
			mower->setCurrentText("Noise Estimator");
		}

		FileHandler fh;
		filetypes_->insertItem("Detect automatically",0);
		for (int i=1; i< FileHandler::SIZE_OF_TYPE; ++i)
		{
			filetypes_->insertItem(fh.typeToName((FileHandler::Type)(i)),i);
		}
	}
	
	OpenDialog::~OpenDialog()
	{
		
	}
	
	void OpenDialog::browse()
	{
		//browse files
		if (source->selected()==dynamic_cast<QButton*>(file_radio))
		{
			String filter_all = "all files (*.dta;*.dta2d;*.DTA;*.DTA2D";
			String filter_single = "dta files (*.dta;*.DTA);;dta2d files (*.dta2d;*.DTA2D)";
	#ifdef ANDIMS_DEF
			filter_single += ";;ANDI/MS files (*.cdf;*.CDF)";
			filter_all +=";*.cdf;*.CDF";
	#endif
			filter_single +=";;mzXML files (*.mzXML;*.mzxml;*.MZXML);;mzData files (*.mzData;*.mzdata;*.MZDATA);;feature map (*.feat;*.FEAT);;all files (*.*)";
			filter_all += ";*.mzXML;*.mzxml;*.MZXML;*.mzData;*.mzdata;*.MZDATA;*.feat;*.FEAT);;" + filter_single;
		
		 	QStringList files = QFileDialog::getOpenFileNames(filter_all.c_str(), prefs_.getValue("Preferences:DefaultPath").toString().c_str(), this,"open file dialog", "Select file(s) to open");
			//check if the dialog was canceled
			if (files.size()!=0)
			{
				ok_button->setEnabled(true);
				names_.clear();
				name_label->setText( "" );
				for(QStringList::iterator it=files.begin();it!=files.end();it++)
				{
					QFileInfo path_info( *it );
					// set path in the first file
					if (it==files.begin())
					{
						path_label->setText( path_info.dirPath() );
					}
					names_.push_back(path_info.filePath().ascii());
					name_label->setText( name_label->text() + path_info.fileName() + " " );
				}
			}
		}
		else 
		//browse DB
		{
	#ifdef DB_DEF
			DBConnection db;
			try
			{
				if (prefs_.getValue("DBPassword").isEmpty())
				{
					stringstream ss;
					ss << "Enter password for user '" << getPref_("Preferences:DB:Login") << "' at '"<< getPref_("Preferences:DB:Host")<<":"<<getPref_("Preferences:DB:Port")<<"' : ";
					bool ok;
					QString text = QInputDialog::getText("TOPPView Database Password", ss.str().c_str(), QLineEdit::Password,QString::null, &ok, this);
					if ( ok )
						{
							prefs_.setValue("DBPassword",text.ascii());
						}
				}
		
				if (!(prefs_.getValue("DBPassword").isEmpty()))
				{
					//cout <<"Name:'"<< getPref_("Preferences:DB:Name") <<"' Login:'"<<getPref_("Preferences:DB:Login")<<"' PW:'"<<prefs_.getValue("DBPassword")<<"' Host:'"<<getPref_("Preferences:DB:Host")<<"' Port:'"<<getPref_("Preferences:DB:Port")<<"'"<<endl;
					db.connect(getPref_("Preferences:DB:Name"), getPref_("Preferences:DB:Login"),getPref_("DBPassword"),getPref_("Preferences:DB:Host"),getPrefAsInt_("Preferences:DB:Port"));
					vector<UnsignedInt> result;
		
					DBSpectrumSelectorDialog dialog(db,result,this);
					if (dialog.exec() && result.size()!=0)
					{
						ok_button->setEnabled(true);
						names_.clear();
						name_label->setText( "" );
						QString convert;
						for (vector<UnsignedInt>::iterator it = result.begin();it!=result.end();++it)
						{
							names_.push_back(String(*it));
							
							name_label->setText( name_label->text() + " " + convert.number(*it) );
						}
					}
				}
			}
			catch (DBConnection::InvalidQuery er)
			{
				prefs_.remove("DBPassword");
				stringstream ss;
				ss << "Unable to log in to the database server.\nCheck the login data in preferences!\n\nDatabase error message:\n"<<er.what();
				QMessageBox::warning ( this, "Connection problem", ss.str().c_str(), QMessageBox::Ok , QMessageBox::NoButton );
		
			}
	#endif
		}
	}
	
	void OpenDialog::showMetadata()
	{
		//TODO implement!
		QMessageBox::information( this, "Show metadata","This feature is not implemented yet!" );
	}
	
	const vector<String>& OpenDialog::getNames() const
	{
		return names_;
	}
	
	OpenDialog::DataSource OpenDialog::getSource() const
	{
		// files
		if (source->selected()==dynamic_cast<QButton*>(file_radio))
		{
			return FILE;
		}
		return DB;
	}
	
	const DataValue& OpenDialog::getPref_(const String& name) const
	{
		return prefs_.getValue(name);
	}
	
	SignedInt OpenDialog::getPrefAsInt_(const String& name) const
	{
		return (SignedInt)(prefs_.getValue(name));
	}
	
	bool OpenDialog::isViewMaps2D() const
	{
		if (open_maps->selected()==dynamic_cast<QButton*>(d2_radio))
		{
			return true;
		}
		return false;
	}
	
	OpenDialog::Mower OpenDialog::getMower() const
	{
		if (mower->currentItem()==0)
		{
			return NO_MOWER;
		}
		return NOISE_ESTIMATOR;
	}
	
	bool OpenDialog::isOpenAsNewTab() const
	{
		if (open_in->selected()==dynamic_cast<QButton*>(newtab_radio))
		{
			return true;
		}
		return false;	
	}

  FileHandler::Type OpenDialog::forcedFileType() const
  {
  	return (FileHandler::Type) (filetypes_->currentItem());
  }

}



