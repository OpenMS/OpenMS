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
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/OpenDialog.h>
#include <OpenMS/VISUAL/DIALOGS/DBSpectrumSelectorDialog.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>


//OpenMS meta data include
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/METADATA/Gradient.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/ProcessingMethod.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/Acquisition.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/SpectrumSettings.h>

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
		
		if (!(bool(getPrefAsInt_("Preferences:DefaultMapView"))))
		{
			d3_radio->setChecked(true);
		}
		if ((String)(getPref_("Preferences:MapIntensityCutoff"))!="None")
		{
			mower->setCurrentIndex(mower->findText("Noise Estimator"));
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
			filter_all += " *.mzXML *.mzData *.feat *.pairs);;" ;
			filter_single +=";;mzXML files (*.mzXML);;mzData files (*.mzData);;feature map (*.feat);;feature pairs (*.pairs);;all files (*.*)";
		
		 	QStringList files = QFileDialog::getOpenFileNames(this, "Open file(s)", prefs_.getValue("Preferences:DefaultPath").toString().c_str(), (filter_all+ filter_single).c_str());
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
					ss << "Enter password for user '" << getPref_("Preferences:DB:Login") << "' at '"<< getPref_("Preferences:DB:Host")<<":"<<getPref_("Preferences:DB:Port")<<"' : ";
					bool ok;
					QString text = QInputDialog::getText(this, "TOPPView Database Password", ss.str().c_str(), QLineEdit::Password,QString::null, &ok);
					if ( ok )
					{
						prefs_.setValue("DBPassword",text.toAscii().data());
					}
				}
		
				if (!(prefs_.getValue("DBPassword").isEmpty()))
				{
					//cout <<"Name:'"<< getPref_("Preferences:DB:Name") <<"' Login:'"<<getPref_("Preferences:DB:Login")<<"' PW:'"<<prefs_.getValue("DBPassword")<<"' Host:'"<<getPref_("Preferences:DB:Host")<<"' Port:'"<<getPref_("Preferences:DB:Port")<<"'"<<endl;
					db.connect(getPref_("Preferences:DB:Name"), getPref_("Preferences:DB:Login"),getPref_("DBPassword"),getPref_("Preferences:DB:Host"),getPrefAsInt_("Preferences:DB:Port"));
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
			FileHandler().loadExperiment(names_[0],exp);
			
			//------------------------------------------------------------------------------
			//	             Set meta data object
			//------------------------------------------------------------------------------
			DateTime date;
			date.now();
			Date today;
			today.today();
			
			//------------------------------------------------------------------------------
			//	             Set ProteinIdentification object
			//------------------------------------------------------------------------------
			ProteinIdentification pid;
			pid.setProteinSignificanceThreshold(3.2);
			pid.setDateTime(date);
			
			std::vector< ProteinHit > prot_hits;  
			ProteinHit *phit1 = new ProteinHit();
			ProteinHit *phit2 = new ProteinHit();
			ProteinHit *phit3 = new ProteinHit();
			ProteinHit *phit4 = new ProteinHit();
			Real score1 =4.4;
			Real score2 =6.4;
			Real score3 =2.4;
			Real score4 =9.4;
			phit1->setScore(score1);
			phit2->setScore(score2);
			phit3->setScore(score3);
			phit4->setScore(score4);
			phit1->setAccession("1BNL");
			phit2->setAccession("1BY0");
			phit3->setAccession("ACC392");
			phit4->setAccession("ACD392");
			prot_hits.push_back(*phit1);
			prot_hits.push_back(*phit2);
			prot_hits.push_back(*phit3);
			prot_hits.push_back(*phit4);
			
			pid.setProteinHits(prot_hits);
			
			std::vector<ProteinIdentification> protIDs;
			protIDs.push_back(pid);
			
			
			
			//------------------------------------------------------------------------------
			//	             Set PeptideHit object
			//------------------------------------------------------------------------------
			PeptideHit pep;
			PeptideHit pep1;
			PeptideHit pep2;
			PeptideHit pep3;
			pep.setScore(1.2);
			pep.setScoreType("XCorr");
			pep.setSequence("ARRAY");
			pep.addProteinIndex(date,"1BNL");
			pep.addProteinIndex(date,"1BY0");
			pep.addProteinIndex(date,"ACC392");
			//pep.setRank(2);
			
			pep1.setScore(2.2);
			pep1.setScoreType("XCorr");
			pep1.setSequence("TTRRAY");
			pep1.addProteinIndex(date,"1BY0");
			pep1.addProteinIndex(date,"ACC392");
			//pep1.setRank(1);
			
			pep2.setScore(6.2);
			pep2.setScoreType("XCorr");
			pep2.setSequence("TTYYRAY");
			pep2.addProteinIndex(date,"1BNL");
			pep2.addProteinIndex(date,"ACC392");
			pep2.addProteinIndex(date,"ACD392");
			//pep2.setRank(3);
			
			pep3.setScore(8.4);
			pep3.setScoreType("XCorr");
			pep3.setSequence("AARRTTYYRAY");
			
			std::vector< PeptideHit > pep_hits; 
			pep_hits.push_back(pep); 
			pep_hits.push_back(pep1); 
			pep_hits.push_back(pep2); 
			pep_hits.push_back(pep3); 
			
			//------------------------------------------------------------------------------
			//	             Set Identification object
			//------------------------------------------------------------------------------
			Identification ident;
			ident.setProteinSignificanceThreshold(3.2);
			ident.setPeptideSignificanceThreshold(4.8);
			ident.setDateTime(date);
			
			
			ident.setPeptideAndProteinHits(pep_hits, prot_hits);
			
			std::vector<Identification> Idents;
			Idents.push_back(ident);
		
			SpectrumSettings ss;
			ss.setType(SpectrumSettings::PEAKS);
			ss.setComment("A lot of space for comments here ...");	
			ss.setIdentifications(Idents);
	
			//-------------------------------------------------------------------------
			//		end of meta data testing
			//-------------------------------------------------------------------------
						
			//MSMetaDataExplorer dlg(false, this);
			MSMetaDataExplorer dlg(true, this);
			dlg.setWindowTitle("Meta data");
			dlg.add(&ss);
			
			//dlg.add(&exp);
			
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



