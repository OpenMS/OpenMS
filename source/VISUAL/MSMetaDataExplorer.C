// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
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


//OpenMS
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/VISUAL/VISUALIZER/SampleVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/DigestionVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ModificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/TaggingVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/HPLCVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/GradientVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/SoftwareVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/SourceFileVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ContactPersonVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/InstrumentVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/IonSourceVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/IonDetectorVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/MassAnalyzerVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ProcessingMethodVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ProteinIdentificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ProteinHitVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/PeptideHitVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/ExperimentalSettingsVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionInfoVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoDescriptionVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/PrecursorVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/InstrumentSettingsVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/PeptideIdentificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/SpectrumSettingsVisualizer.h>

#include <QtGui/QLabel>
#include <QtGui/QStackedWidget>
#include <QtGui/QHBoxLayout>
#include <QtGui/QPushButton>
#include <QtGui/QMessageBox>
#include <QtGui/QSplitter>

class QLayoutItem;

using namespace std;


//class QtGui/QLayoutItem;

namespace OpenMS
{
	
	MSMetaDataExplorer::MSMetaDataExplorer(bool editable, QWidget *parent, bool modal)
	: QDialog(parent), 
		editable_(editable)
	{    
	  setModal(modal);
	  
	  //splitter to hold treewidget and the right part (on a dummy widget)
	  QGridLayout* grid = new QGridLayout(this);
	  QSplitter* splitter = new QSplitter(Qt::Horizontal,this);
	  grid->addWidget(splitter,0,0);
		
		//Create the tree for exploring data 
		treeview_ = new QTreeWidget(this);
		treeview_->setColumnCount(3);
		treeview_->setHeaderLabel("Browse in Metadata tree");
		treeview_->setRootIsDecorated(true);
		treeview_->setColumnHidden(1,true);
		treeview_->setColumnHidden(2,true);
	 	splitter->addWidget(treeview_);
		
		//create a dummy widget to hold a grid layout and the item stack
		QWidget* dummy = new QWidget(splitter);
		splitter->addWidget(dummy);
		grid = new QGridLayout(dummy);
		grid->setColumnStretch(0,1);
		
		//Create WidgetStack for managing visible metadata
		ws_ = new QStackedWidget(dummy);    
		grid->addWidget(ws_,0,0,1,3);
				
		if(isEditable())
		{
			saveallbutton_ = new QPushButton("OK", dummy);
		 	cancelbutton_ = new QPushButton("Cancel", dummy);
			grid->addWidget(saveallbutton_, 1, 1);
			grid->addWidget(cancelbutton_,1 ,2);
			connect(saveallbutton_, SIGNAL(clicked()), this, SLOT(saveAll_())  );
			connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject())  );
		}
		else 
		{
			closebutton_ = new QPushButton("Close", dummy);
			grid->addWidget(closebutton_,1 ,2);
			connect(closebutton_, SIGNAL(clicked()), this, SLOT(reject())  );
		}
		
	  connect(treeview_, SIGNAL(itemClicked(QTreeWidgetItem*,int)), this, SLOT(showDetails_(QTreeWidgetItem*,int))  );
	   
		status_list_="";
		
	}//end of constructor
	
	
	bool MSMetaDataExplorer::isEditable()
	{
			return editable_;
	}
	
	void MSMetaDataExplorer::setStatus(std::string status)
	{
		status_list_ = status_list_+ "\n"+ status;    
	}
	
	
	void MSMetaDataExplorer::showDetails_(QTreeWidgetItem *item,int /*column*/)
	{
	  ws_->setCurrentIndex(item->text(1).toInt());
	}
	
	void MSMetaDataExplorer::connectVisualizer_(BaseVisualizer* ptr)
	{
		connect(ptr, SIGNAL(sendStatus(std::string)), this, SLOT(setStatus(std::string))  );			
	}
		
	//Save all changes
	void MSMetaDataExplorer::saveAll_()
	{
		try
		{
			//call internal store function of all active visualizer objects
			for (int i = 0; i < ws_->count(); ++i) 
			{
			  dynamic_cast<BaseVisualizer*>(ws_->widget(i))->store_();
			}
			if(status_list_.length() != 0)
			{
				status_list_ = status_list_ + "\n"+ "\n"+ "Invalid modifications will not be saved.";
				QMessageBox::warning(this,tr("Save warning"),status_list_.c_str());
			}
			else
			{
				QMessageBox::information(this,tr("Save information"),tr("Your modifications have been saved."));
			}
		}
		catch(exception& e)
		{
			cout<<"Exception while trying to save modifications."<<endl<<e.what()<<endl;
		}
		
		//close dialog
		accept();	
	}
	
	//-------------------------------------------------------------------------------
	//	overloaded visualize functions to call the corresponding data visualizer
	//	(in alphabetical oeder)
	//-------------------------------------------------------------------------------

	//Visualizing AcquisitionInfo object
	void MSMetaDataExplorer::visualize_(AcquisitionInfo& meta, QTreeWidgetItem* parent)
	{
		AcquisitionInfoVisualizer *visualizer = new AcquisitionInfoVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "Acquisition Info" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
			
		//Get Aquisition objects
		for(UInt i=0; i< meta.size(); ++i)
		{ 
			visualize_(meta[i], item);
		}
		
		connectVisualizer_(visualizer);
	}
		
	//Visualizing Acquisition object
	void MSMetaDataExplorer::visualize_(Acquisition& meta, QTreeWidgetItem* parent)
	{
		AcquisitionVisualizer *visualizer = new AcquisitionVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    String name = String("Acquisition Nr. ") + meta.getNumber();
    labels << name.c_str() << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing ContactPerson object
	void MSMetaDataExplorer::visualize_(ContactPerson& meta, QTreeWidgetItem* parent)
	{
		ContactPersonVisualizer *visualizer = new ContactPersonVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "ContactPerson" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing Digestion object
	void MSMetaDataExplorer::visualize_(Digestion& meta, QTreeWidgetItem* parent)
	{
		DigestionVisualizer *visualizer = new DigestionVisualizer(isEditable(),this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "Digestion" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item); 
		connectVisualizer_(visualizer);
	}
	
	//Visualizing ExperimentalSettings object
	void MSMetaDataExplorer::visualize_(ExperimentalSettings& meta, QTreeWidgetItem* parent)
	{
		ExperimentalSettingsVisualizer *visualizer = new ExperimentalSettingsVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "ExperimentalSettings" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		
		
		//check for Sample
		visualize_(meta.getSample(), item);
		
		//check for ProteinIdentification
		for(UInt i=0; i<meta.getProteinIdentifications().size(); ++i)    
		{
			visualize_(meta.getProteinIdentifications()[i], item);
		}
		
		//check for ProcessingMethod
		visualize_(meta.getProcessingMethod(), item);
		
		//check for Instrument
		visualize_(meta.getInstrument(), item);
		
		//check for SourceFile
		visualize_(meta.getSourceFile(), item);
		
		//check for ContactPersons
		for(UInt i=0; i<meta.getContacts().size(); ++i)    
		{
			visualize_(meta.getContacts()[i], item);
		}
		
		//check for Software
		visualize_(meta.getSoftware(), item);
		
		//check for HPLC
		visualize_(meta.getHPLC(), item);
					
		connectVisualizer_(visualizer);
	}
	
	//Visualizing Gradient object
	void MSMetaDataExplorer::visualize_(Gradient& meta, QTreeWidgetItem* parent)
	{
		GradientVisualizer *visualizer = new GradientVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "Gradient" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing HPLC object
	void MSMetaDataExplorer::visualize_(HPLC& meta, QTreeWidgetItem* parent)
	{
		HPLCVisualizer *visualizer = new HPLCVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "HPLC" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
		
		
		visualize_(meta.getGradient(), item);
		
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing PeptideIdentification object
	void MSMetaDataExplorer::visualize_(PeptideIdentification& meta, QTreeWidgetItem* parent)
	{
		PeptideIdentificationVisualizer *visualizer = new PeptideIdentificationVisualizer(isEditable(), this, this); 
		
    QStringList labels;
    int id = ws_->addWidget(visualizer);
    labels << QString("PeptideIdentification %1").arg(meta.getScoreType().c_str()) << QString::number(id);

		visualizer->load(meta,id);  

    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
			
		//check for proteins and peptides hits
		meta.assignRanks();
			
		//list all peptides hits in the tree
		for(UInt i=0; i<meta.getHits().size(); ++i)    
		{
			visualize_(const_cast<PeptideHit&>(meta.getHits()[i]), item);
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);

		connectVisualizer_(visualizer);
	}

	//Visualizing ProteinIdentification object
	void MSMetaDataExplorer::visualize_(ProteinIdentification& meta, QTreeWidgetItem* parent)
	{
		ProteinIdentificationVisualizer *visualizer = new ProteinIdentificationVisualizer(isEditable(), this, this); 
		
    QStringList labels;
    int id = ws_->addWidget(visualizer);
    labels << QString("ProteinIdentification %1").arg(meta.getSearchEngine().c_str()) << QString::number(id);
    
    visualizer->load(meta,id);  
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		//check for proteinhits objects
		meta.assignRanks();

		for(UInt i=0; i<meta.getHits().size(); ++i)    
		{
			visualize_(const_cast<ProteinHit&>(meta.getHits()[i]), item);
		}
		
		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing InstrumentSettings object
	void MSMetaDataExplorer::visualize_(InstrumentSettings& meta, QTreeWidgetItem* parent)
	{
		InstrumentSettingsVisualizer *visualizer = new InstrumentSettingsVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "InstrumentSettings" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing Instrument object
	void MSMetaDataExplorer::visualize_(Instrument& meta, QTreeWidgetItem* parent)
	{
		InstrumentVisualizer* visualizer = new InstrumentVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "Instrument" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);

		//visualize IonSource object
		visualize_(meta.getIonSource(), item);
		
		//visualize IonDetector object
		visualize_(meta.getIonDetector(), item);
			
		//Check for MassAnalyzers
		vector<MassAnalyzer>& v= meta.getMassAnalyzers();  
		if(v.size() != 0)
		{
			for(UInt i=0; i<v.size(); ++i)    
			{
				visualize_(v[i], item);
			}
		}
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing IonDetector object
	void MSMetaDataExplorer::visualize_(IonDetector& meta, QTreeWidgetItem* parent)
	{
		IonDetectorVisualizer *visualizer = new IonDetectorVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "IonDetector" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing IonSource object
	void MSMetaDataExplorer::visualize_(IonSource& meta, QTreeWidgetItem* parent)
	{
		IonSourceVisualizer *visualizer = new IonSourceVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "IonSource" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);		
		connectVisualizer_(visualizer);	
	}
	
	
	//Visualizing MassAnalyzer object
	void MSMetaDataExplorer::visualize_(MassAnalyzer& meta, QTreeWidgetItem* parent)
	{
		MassAnalyzerVisualizer *visualizer = new MassAnalyzerVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "MassAnalyzer" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing MetaInfoDescription object
	void MSMetaDataExplorer::visualize_(MetaInfoDescription& meta,  QTreeWidgetItem* parent, const String& key)
	{
		MetaInfoDescriptionVisualizer *visualizer = new MetaInfoDescriptionVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    String name = String ("MetaInfoDescription ") + key;
    labels << name.c_str() << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
		
		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
			
		//Check for source file
		visualize_(meta.getSourceFile(), item);
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing MetaInfoInterface object
	void MSMetaDataExplorer::visualize_(MetaInfoInterface& meta, QTreeWidgetItem* parent)
	{
		MetaInfoVisualizer *visualizer = new MetaInfoVisualizer(isEditable(),this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "MetaInfo" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing modification object
	void MSMetaDataExplorer::visualize_(Modification& meta, QTreeWidgetItem* parent)
	{
		ModificationVisualizer *visualizer = new ModificationVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "Modification" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item); 
		connectVisualizer_(visualizer); 
	}
	
	
	//Visualizing PeptideHit object
	void MSMetaDataExplorer::visualize_(PeptideHit& meta, QTreeWidgetItem* parent)
	{
		PeptideHitVisualizer *visualizer = new PeptideHitVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
		String name = String("Pep ") + meta.getSequence() + " (" + meta.getScore() + ')';
		QString qs_name( name.c_str() );
    
		QStringList labels;
    labels << qs_name << QString::number(ws_->addWidget(visualizer)) << QString::number(meta.getScore());
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		} 

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);

		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing Precursor object
	void MSMetaDataExplorer::visualize_(Precursor& meta, QTreeWidgetItem* parent)
	{
		PrecursorVisualizer* visualizer = new PrecursorVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "Precursor" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing ProcessingMethod object
	void MSMetaDataExplorer::visualize_(ProcessingMethod& meta, QTreeWidgetItem* parent)
	{
		ProcessingMethodVisualizer *visualizer = new ProcessingMethodVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "ProcessingMethod" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}

		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing ProteinHit object
	void MSMetaDataExplorer::visualize_(ProteinHit& meta, QTreeWidgetItem* parent)
	{
		ProteinHitVisualizer *visualizer = new ProteinHitVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
		String name = String("Prot ") + meta.getAccession() + " (" + meta.getScore() + ')';
		QString qs_name( name.c_str() );
		
    QStringList labels;
    labels << qs_name << QString::number(ws_->addWidget(visualizer)) << QString::number(meta.getScore());
		
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		} 
		
		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		
		connectVisualizer_(visualizer);
	}
	
		
	//Visualizing sample object
	void MSMetaDataExplorer::visualize_(Sample& meta, QTreeWidgetItem* parent )
	{	
		SampleVisualizer *visualizer = new SampleVisualizer(isEditable(), this); 
    visualizer->load(meta);
    
    QStringList labels;
    labels << (String("Sample ") + meta.getName()).c_str() << QString::number(ws_->addWidget(visualizer));
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
				
		//check for treatments
		if(meta.countTreatments() != 0)
		{	
			for(Int i=0; i<meta.countTreatments(); ++i)
			{
				if(meta.getTreatment(i).getType()=="Digestion")
				{
					visualize_((const_cast<Digestion&>(dynamic_cast<const Digestion&>(meta.getTreatment(i)) ) ), item );
				}
				else if(meta.getTreatment(i).getType()== "Modification")
				{
					//Cast SampleTreatment reference to a const modification reference
					visualize_((const_cast<Modification&>(dynamic_cast<const Modification&>(meta.getTreatment(i))) ), item );
										
				}
				else if(meta.getTreatment(i).getType()=="Tagging")
				{
					visualize_((const_cast<Tagging&>(dynamic_cast<const Tagging&>(meta.getTreatment(i)) )), item );
				}
			}
		}
		
		//Check for subsamples
		vector<Sample>& v = meta.getSubsamples();  
		
		if(v.size() != 0)
		{
			for(UInt i=0; i<v.size(); ++i)    
			{
				visualize_(v[i], item);
			}
		}
		
		//check for metainfo objects
		visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		connectVisualizer_(visualizer);
	}
	
	
	//Visualizing Software object
	void MSMetaDataExplorer::visualize_(Software& meta, QTreeWidgetItem* parent)
	{
		SoftwareVisualizer *visualizer = new SoftwareVisualizer(isEditable(), this);
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "Software" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
		connectVisualizer_(visualizer);
	}

	
	//Visualizing SourceFile object
	void MSMetaDataExplorer::visualize_(SourceFile& meta, QTreeWidgetItem* parent)
	{
		SourceFileVisualizer *visualizer = new SourceFileVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "SourceFile" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		} 
		connectVisualizer_(visualizer);
	}
	
	
	
	//Visualizing SpectrumSettings object
	void MSMetaDataExplorer::visualize_(SpectrumSettings& meta, QTreeWidgetItem* parent)
	{
		SpectrumSettingsVisualizer *visualizer = new SpectrumSettingsVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "SpectrumSettings" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
		
		
		//check for InstrumentSettings
		visualize_(meta.getInstrumentSettings(), item);
		
		//check for PeptideIdentification
		for(UInt i=0; i<meta.getPeptideIdentifications().size(); ++i)    
		{
			visualize_(meta.getPeptideIdentifications()[i], item);
		}
		
		//check for Precursor
		visualize_(meta.getPrecursor(), item);
				
		//check for MetaInfoDescription
		for(map<String, MetaInfoDescription>::iterator iter = meta.getMetaInfoDescriptions().begin(); iter != meta.getMetaInfoDescriptions().end(); iter++ ) 
		{
			visualize_(iter->second, item, iter->first );
		}			
		
		//check for AcquisitionInfo
		visualize_(meta.getAcquisitionInfo(), item);
			
		connectVisualizer_(visualizer);
	}
	
	
	
	//Visualizing tagging object
	void MSMetaDataExplorer::visualize_(Tagging& meta, QTreeWidgetItem* parent)
	{
		TaggingVisualizer *visualizer = new TaggingVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
    QStringList labels;
    labels << "Tagging" << QString::number(ws_->addWidget(visualizer));
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		}
		connectVisualizer_(visualizer);
	}

	void MSMetaDataExplorer::filterHits_(DoubleReal threshold, bool higher_better, int tree_item_id)
	{		
		// find item in tree belonging to PeptideIdentification object
		QTreeWidgetItem *item = treeview_->findItems(QString::number(tree_item_id),Qt::MatchExactly | Qt::MatchRecursive, 1).first();
	
		//set the items visible or not visible depending to their score and the current threshold
		for(int i=0; i<item->childCount(); ++i)
		{
			QTreeWidgetItem* child = item->child(i);
	
			if( (higher_better && child->text(2).toFloat() <= threshold) || (!higher_better && child->text(2).toFloat() >= threshold) )
			{
				child->setHidden(true);
			}
			else
			{
				child->setHidden(false);
			}
		}
	
		//parent item must be collapsed and re-expanded so the items will be shown...
		treeview_->collapseItem(item);
		treeview_->expandItem(item);
	}			

	/// Filters hits according to a score @a threshold. Takes the score orientation into account
	void MSMetaDataExplorer::showAllHits_(int tree_item_id)
	{
		// find item in tree belonging to PeptideIdentification object
		QTreeWidgetItem *item = treeview_->findItems(QString::number(tree_item_id),Qt::MatchExactly | Qt::MatchRecursive, 1).first();
	
		//set the items visible or not visible depending to their score and the current threshold
		for(int i=0; i<item->childCount(); ++i)
		{
			item->child(i)->setHidden(false);
		}
	
		//parent item must be collapsed and re-expanded so the items will be shown...
		treeview_->collapseItem(item);
		treeview_->expandItem(item);
	}


}//end of MSMetaDataExplorer






