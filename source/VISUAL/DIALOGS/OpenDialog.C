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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/FORMAT/DB/DBAdapter.h>
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/OpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/DBSpectrumSelectorDialog.h>
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

// STL includes
#include <iostream>

using namespace std;

namespace OpenMS
{

	OpenDialog::OpenDialog( Param& preferences, QWidget * parent)
		: QDialog(parent),
			prefs_(preferences)
	{
		setupUi(this);
		
		if ((String)(getPref_("preferences:default_map_view"))=="3d")
		{
			d3_radio->setChecked(true);
		}
		if ((String)(getPref_("preferences:intensity_cutoff"))!="none")
		{
			mower->setCurrentIndex(mower->findText("noise_estimator"));
		}

		FileHandler fh;
		filetypes_->insertItem(0,"Detect automatically",0);
		for (int i=1; i< FileHandler::SIZE_OF_TYPE; ++i)
		{
			filetypes_->insertItem(filetypes_->count(),fh.typeToName((FileHandler::Type)(i)).c_str(),i);
		}
	}
	
	OpenDialog::~OpenDialog()
	{
		
	}
	
	void OpenDialog::browse()
	{
		//browse files
		if (file_radio->isChecked())
		{
			String filter_all = "readable files (*.dta *.dta2d";
			String filter_single = "dta files (*.dta);;dta2d files (*.dta2d)";
	#ifdef ANDIMS_DEF
			filter_all +=" *.cdf";
			filter_single += ";;ANDI/MS files (*.cdf)";
	#endif
			filter_all += " *.mzXML *.mzData *.featureXML *.featurePairsXML);;" ;
			filter_single +=";;mzXML files (*.mzXML);;mzData files (*.mzData);;feature map (*.featureXML);;feature pairs (*.pairs);;all files (*.*)";
		
		 	QStringList files = QFileDialog::getOpenFileNames(this, "Open file(s)", prefs_.getValue("preferences:default_path").toString().c_str(), (filter_all+ filter_single).c_str());
			//check if the dialog was canceled
			if (files.size()!=0)
			{
				ok_button->setEnabled(true);
				names_.clear();
				name_label->setText( "" );
				for(QStringList::iterator it=files.begin();it!=files.end();it++)
				{
					// set path in the first file
					if (it==files.begin())
					{
						path_label->setText( File::path((*it).toAscii().data()).c_str() );
					}
					names_.push_back((*it).toAscii().data());
					name_label->setText( (String(name_label->text().toAscii().data()) + File::basename((*it).toAscii().data()) + " ").c_str() );
				}
			}
		}
		else 
		//browse DB
		{
			DBConnection db;
			try
			{
				if (prefs_.getValue("DBPassword").isEmpty())
				{
					stringstream ss;
					ss << "Enter password for user '" << getPref_("preferences:db:login") << "' at '"<< getPref_("preferences:db:host")<<":"<<getPref_("preferences:db:port")<<"' : ";
					bool ok;
					QString text = QInputDialog::getText(this, "TOPPView Database Password", ss.str().c_str(), QLineEdit::Password,QString::null, &ok);
					if ( ok )
					{
						prefs_.setValue("DBPassword",text.toAscii().data());
					}
				}
		
				if (!(prefs_.getValue("DBPassword").isEmpty()))
				{
					//cout <<"Name:'"<< getPref_("preferences:db:name") <<"' Login:'"<<getPref_("preferences:db:login")<<"' PW:'"<<prefs_.getValue("DBPassword")<<"' Host:'"<<getPref_("preferences:db:host")<<"' Port:'"<<getPref_("preferences:db:port")<<"'"<<endl;
					db.connect(getPref_("preferences:db:name"), getPref_("preferences:db:login"),getPref_("DBPassword"),getPref_("preferences:db:host"),getPrefAsInt_("preferences:db:port"));
					vector<UInt> result;
		
					DBSpectrumSelectorDialog dialog(db,result,this);
					if (dialog.exec() && result.size()!=0)
					{
						ok_button->setEnabled(true);
						names_.clear();
						name_label->setText( "" );
						QString convert;
						for (vector<UInt>::iterator it = result.begin();it!=result.end();++it)
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
				QMessageBox::warning ( this, "Connection problem", ss.str().c_str(), QMessageBox::Ok , Qt::NoButton );
		
			}
		}
	}
	
	void OpenDialog::showMetadata()
	{
		if (names_.size()!=0)
		{
			MSExperiment<> exp;

			//browse files
			if (file_radio->isChecked())
			{
				FileHandler fh;
				fh.getOptions().setMetadataOnly(true);
				fh.loadExperiment(names_[0],exp);
			}
			else
			{
				DBConnection con;
				con.connect(getPref_("preferences:db:name"), getPref_("preferences:db:login"),getPref_("DBPassword"),getPref_("preferences:db:host"),getPrefAsInt_("preferences:db:port"));
				DBAdapter db(con);
				db.getOptions().setMetadataOnly(true);
				db.loadExperiment(names_[0].toInt(), exp);
			}
			
			MSMetaDataExplorer dlg(false, this);
			dlg.setWindowTitle("Meta data");			
			dlg.visualize(exp);
   	 	dlg.exec();
		}
	}
	
	const vector<String>& OpenDialog::getNames() const
	{
		return names_;
	}
	
	OpenDialog::DataSource OpenDialog::getSource() const
	{
		// files
		if (file_radio->isChecked())
		{
			return FILE;
		}
		return DB;
	}
	
	const DataValue& OpenDialog::getPref_(const String& name) const
	{
		return prefs_.getValue(name);
	}
	
	Int OpenDialog::getPrefAsInt_(const String& name) const
	{
		return (Int)(prefs_.getValue(name));
	}
	
	bool OpenDialog::isViewMaps2D() const
	{
		if (d2_radio->isChecked())
		{
			return true;
		}
		return false;
	}
	
	OpenDialog::Mower OpenDialog::getMower() const
	{
		if (mower->currentIndex()==0)
		{
			return NO_MOWER;
		}
		return NOISE_ESTIMATOR;
	}
	
	bool OpenDialog::isOpenAsNewTab() const
	{
		if (newtab_radio->isChecked())
		{
			return true;
		}
		return false;	
	}

  FileHandler::Type OpenDialog::forcedFileType() const
  {
  	return (FileHandler::Type) (filetypes_->currentIndex());
  }

}



