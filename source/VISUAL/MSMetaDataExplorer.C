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
// $Maintainer: stefan_heess $
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
#include <OpenMS/VISUAL/VISUALIZER/IdentificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/SpectrumSettingsVisualizer.h>

#include <QtGui/QLabel>
#include <QtGui/QTreeWidget>
#include <QtGui/QStackedWidget>
#include <QtGui/QHBoxLayout>
#include <QtGui/QPushButton>

using namespace std;

namespace OpenMS
{
	
	MSMetaDataExplorer::MSMetaDataExplorer(bool editable, QWidget *parent, bool modal)
	: QDialog(parent), 
		editable_(editable)
	{    
	  setModal(modal);
	  
		//basiclayout_ = new QHBoxLayout(this);
		basiclayout_ = new QGridLayout(this);
		basiclayout_->setSpacing(6);
	  basiclayout_->setMargin(11); 
		
		vertlayout_ = new QVBoxLayout();
		vertlayout_->setSpacing(6);
	  vertlayout_->setMargin(11); 
		
		buttonlayout_ = new QHBoxLayout();
		buttonlayout_->setSpacing(6);
		buttonlayout_->setMargin(11); 
	  		
		//Create the tree for exploring data 
		treeview_ = new QTreeWidget(this);
		treeview_->setMinimumWidth(250);
		treeview_->setColumnCount(2);
		treeview_->setHeaderLabel("Browse in Metadata tree");
		treeview_->setRootIsDecorated(true);
	  treeview_->sortByColumn(2, Qt::AscendingOrder);
	  
		basiclayout_->addWidget(treeview_, 0,0,3,1);
		basiclayout_->addLayout(vertlayout_, 0,1,3,1);
		  
		//Create WidgetStack for managing visible metadata
		ws_ = new QStackedWidget(this);
		QLabel* hline = new QLabel(this);
		
		vertlayout_->addWidget(ws_);
		vertlayout_->addWidget(hline);
		vertlayout_->addLayout(buttonlayout_);	
		
		if(isEditable())
		{
		
			saveallbutton_ = new QPushButton("OK", this);
		  cancelbutton_ = new QPushButton("Cancel", this);
			buttonlayout_->addStretch(1);
			buttonlayout_->addWidget(saveallbutton_ );
			buttonlayout_->addWidget(cancelbutton_ );
			connect(saveallbutton_, SIGNAL(clicked()), this, SLOT(saveAll_())  );
			connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject())  );
		}
		else 
		{
		
			closebutton_ = new QPushButton("Close", this);
			buttonlayout_->addStretch(1);
			buttonlayout_->addWidget(closebutton_ );
			connect(closebutton_, SIGNAL(clicked()), this, SLOT(reject())  );
		}
		
	  connect(treeview_, SIGNAL(itemClicked(QTreeWidgetItem*,int)), this, SLOT(showDetails_(QTreeWidgetItem*,int))  );
	  
	}//end of constructor
	
	
	
	bool MSMetaDataExplorer::isEditable()
	{
			return editable_;
	}
	
	
	void MSMetaDataExplorer::showDetails_(QTreeWidgetItem *item,int /*column*/)
	{
	  ws_->setCurrentIndex(item->text(1).toInt());
	}
	
	//Save all changes
	void MSMetaDataExplorer::saveAll_()
	{
		try
		{
			//call internal store function of all active visualizer objects
			const QObjectList& children = ws_->children(); 
			for (int i = 0; i < children.size(); ++i) 
			{
			  dynamic_cast<BaseVisualizer*>(children.at(i))->store();
			}
			
			//close dialog
			this->accept();
		}
		catch(exception& e)
		{
			cout<<"Exception: "<<e.what()<<endl;
		}
	}
	
	//--------------------------------------------------------------------------
	//	Functions modifying the listview
	//--------------------------------------------------------------------------
	
	void MSMetaDataExplorer::updateProteinHits_(ProteinIdentification pid, int tree_item_id)
	{
		float threshold = pid.getProteinSignificanceThreshold();
		
		// find item in tree belonging to ProteinIdentification object
		QTreeWidgetItem* item = treeview_->findItems(QString::number(tree_item_id),Qt::MatchExactly, 1).first();	
		item->sortChildren(2, Qt::AscendingOrder);
		
		//set the items visible or not visible depending to their score and the current threshold
		for(int i=0; i<item->childCount(); ++i)
		{
			QTreeWidgetItem* child = item->child(i);
			if( child->text(2).toFloat() <=threshold)
			{
					child->setHidden(true);
			}
			else
			{
				child->setHidden(false);
			}
		}
	}
	
	
	void MSMetaDataExplorer::updateNonRefPeptideHits_(Identification id, int tree_item_id)
	{
		vector< ProteinHit > prots = id.getProteinHits(); 
		
		String date_time;
		id.getDateTime().get(date_time);
		
		multimap< String, ProteinHit > map;
		for(vector<ProteinHit>::iterator it = prots.begin(); it != prots.end();	it++)
		{
			map.insert(make_pair(date_time, *it));
		}
		
		vector<PeptideHit>* hits = id.getNonReferencingHits(map);
				
		// find item in tree belonging to Identification object
		QTreeWidgetItem *item = treeview_->findItems(QString::number(tree_item_id),Qt::MatchExactly, 1).first();	
			
		//set the items visible or not visible depending of non referencing peptide hits
		for(int i=0; i<item->childCount(); ++i)
		{
			QTreeWidgetItem* child = item->child(i);
			
			for(UnsignedInt i=0; i<hits->size(); ++i)
			{
				String seq = (*hits)[i].getSequence();
				String name = String("Pep ") + seq + " (" + (*hits)[i].getScore() + ")";
				
				if( child->text(0)==name.c_str() && child->text(3)== seq.c_str())
				{
					child->setHidden(false);
					break;
				}
				else
				{
					child->setHidden(true);
				}
			}
		}
	}
	
	void MSMetaDataExplorer::updatePeptideHits_(Identification id, int tree_item_id)
	{
		float threshold = id.getPeptideSignificanceThreshold();
				
		// find item in tree belonging to Identification object
		QTreeWidgetItem *item = treeview_->findItems(QString::number(tree_item_id),Qt::MatchExactly, 1).first();	
		
		//set the items visible or not visible depending to their score and the current threshold
		for(int i=0; i<item->childCount(); ++i)
		{
			QTreeWidgetItem* child = item->child(i);
			
			if( (child->text(2)).toFloat() <= threshold)
			{
				child->setHidden(true);
			}
			else
			{
				child->setHidden(false);
			}
		}
	}
	
	void MSMetaDataExplorer::updateRefPeptideHits_(Identification id, int tree_item_id, String ref_date, String ref_acc)
	{
	  if(ref_date.trim()=="" && ref_acc.trim()=="")
		{
			return;
		}

		id.sort();
				
		// find item in tree belonging to Identification object
		QTreeWidgetItem *item = treeview_->findItems(QString::number(tree_item_id),Qt::MatchExactly, 1).first();	 
		
		//Set all items to not visible
		for(int i=0; i<item->childCount(); ++i)
		{
			 item->child(i)->setHidden(true);
		}
		
		//search all peptide hits
		vector< PeptideHit >& peps= id.getPeptideHits();
		bool hit = false;
		for(UnsignedInt i=0; i<peps.size(); ++i)
		{
			vector<pair<String, String> >& protlist = peps[i].getProteinIndices();
			
			String name = String("Pep ") + peps[i].getSequence() + " (" + peps[i].getScore() + ")";
			//search all protein hits belonging to current peptide hit
			for(vector< pair<String, String> >::iterator it = protlist.begin(); it != protlist.end();	it++)
			{
				//Check if peptide hit refers to a certain protein hit.
				if( (ref_date.trim()=="" && it->second == ref_acc )|| 
				    (it->first==ref_date && ref_acc.trim()=="") ||
						(it->first==ref_date && it->second == ref_acc)   )
				{
					hit = true;
					break;	
				}
			}
			
			if(hit)
			{
				//Search QListviewItem belonging to peptide and set it visible.
				for(int i=0; i<item->childCount(); ++i)
				{
					QTreeWidgetItem* child = item->child(i);
					
					if( child->text(0)==name.c_str() && child->text(3)== peps[i].getSequence().c_str() )
					{
						child->setHidden(false);
					}					
				}	
			}
			hit = false;
		}
	}
	
	
	//-------------------------------------------------------------------------------
	//	overloaded visualize functions to call the corresponding data visualizer
	//-------------------------------------------------------------------------------
	
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
			for(SignedInt i=0; i<meta.countTreatments(); ++i)
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
			for(UnsignedInt i=0; i<v.size(); ++i)    
			{
				visualize_(v[i], item);
			}
		}
		
		//check for metainfo objects
		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		} 
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}  
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
		
		try
		{	
			visualize_(meta.getGradient(), item);
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize Gradient. "<<e.what()<<endl;
		}
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}

		//visualize IonSource object
		try
		{	
			visualize_(meta.getIonSource(), item);
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize IonSource. "<<e.what()<<endl;
		}
		
		//visualize IonDetector object
		try
		{	
			visualize_(meta.getIonDetector(), item);
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize IonDetector. "<<e.what()<<endl;
		}
			
		//Check for MassAnalyzers
		vector<MassAnalyzer>& v= meta.getMassAnalyzers();  
		if(v.size() != 0)
		{
			for(UnsignedInt i=0; i<v.size(); ++i)    
			{
				visualize_(v[i], item);
			}
		}
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}			
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}
	}
	
	//Visualizing MassAnalyzer object
	void MSMetaDataExplorer::visualize_(MassAnalyzer& meta, QTreeWidgetItem* parent)
	{
		MassAnalyzerVisualizer *visualizer = new MassAnalyzerVisualizer(isEditable(), this);  
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}
	}
	
	//Visualizing ProteinIdentification object
	void MSMetaDataExplorer::visualize_(ProteinIdentification& meta, QTreeWidgetItem* parent)
	{
		ProteinIdentificationVisualizer *visualizer = new ProteinIdentificationVisualizer(splitvert_, this); 
		
    QStringList labels;
    int id = ws_->addWidget(visualizer);
    labels << "ProteinIdentification" << QString::number(id);
    
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
		meta.sort();
		vector< ProteinHit > v= meta.getProteinHits();  
		if(v.size() != 0)
		{
			for(UnsignedInt i=0; i<v.size(); ++i)    
			{
				visualize_(v[i], item);
			}
		}
	}
	
	//Visualizing Identification object
	void MSMetaDataExplorer::visualize_(Identification& meta, QTreeWidgetItem* parent)
	{
		IdentificationVisualizer *visualizer = new IdentificationVisualizer(splitvert_, this); 
		
    QStringList labels;
    int id = ws_->addWidget(visualizer);
    labels << "Identification" << QString::number(id);

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
		meta.sort();
		vector< ProteinHit > prots = meta.getProteinHits(); 
		vector< PeptideHit > peps = meta.getPeptideHits();  
			
		//list all peptides hits in the tree
		if(peps.size() != 0)
		{
			for(UnsignedInt i=0; i<peps.size(); ++i)    
			{
				visualize_(peps[i], item, &prots);
			}
		}
	}
	
	//Visualizing ProteinHit object
	void MSMetaDataExplorer::visualize_(ProteinHit& meta, QTreeWidgetItem* parent)
	{
		ProteinHitVisualizer *visualizer = new ProteinHitVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
		String name = String("Prot ") + meta.getAccession() + " (" + meta.getScore() + ')';
    QStringList labels;
    labels << "ProcessingMethod" << QString::number(ws_->addWidget(visualizer)) << QString::number(meta.getScore());
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		} 
	}
	
	//Visualizing PeptideHit object
	void MSMetaDataExplorer::visualize_(PeptideHit& meta, QTreeWidgetItem* parent, vector<ProteinHit> *prots )
	{
		PeptideHitVisualizer *visualizer = new PeptideHitVisualizer(isEditable(), this);  
		visualizer->load(meta);  
		
		String name = String("Pep ") + meta.getSequence() + " (" + meta.getScore() + ')';
    QStringList labels;
    labels << "ProcessingMethod" << QString::number(ws_->addWidget(visualizer)) << QString::number(meta.getScore());
    
    QTreeWidgetItem* item;
		if(parent == 0)
		{
			item = new QTreeWidgetItem(treeview_, labels );
		}
		else
		{
			item = new QTreeWidgetItem(parent, labels );
		} 
			
		//get all protein hits for the peptide hits
		vector<pair<String, String> > protlist = meta.getProteinIndices();
		for(vector< pair<String, String> >::iterator it = protlist.begin(); it != protlist.end();	it++)
		{
			//find all protein hits from the indices and visualize them.						
			for(vector<ProteinHit>::iterator prot_it = prots->begin(); prot_it != prots->end(); prot_it++)
			{
				if(it->second == prot_it->getAccession()  )
				{
					visualize_(*(prot_it), item);
				}
			}
		}
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}
		
		try
		{
			//check for Sample
			visualize_(meta.getSample(), item);
			
			//check for ProteinIdentification
			for(UnsignedInt i=0; i<meta.getProteinIdentifications().size(); ++i)    
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
			for(UnsignedInt i=0; i<meta.getContacts().size(); ++i)    
			{
				visualize_(meta.getContacts()[i], item);
			}
			
			//check for Software
			visualize_(meta.getSoftware(), item);
			
			//check for HPLC
			visualize_(meta.getHPLC(), item);
			
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize ExperimentalSettings. "<<e.what()<<endl;
		}
	}
	
	
	
	//-----------------------------------------------------------------------------
	//			Spectrum Settings
	//-----------------------------------------------------------------------------
	
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
		
		try
		{
			//check for InstrumentSettings
			visualize_(meta.getInstrumentSettings(), item);
			
			//check for Identification
			for(UnsignedInt i=0; i<meta.getIdentifications().size(); ++i)    
			{
				visualize_(meta.getIdentifications()[i], item);
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
		
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize SpectrumSettings. "<<e.what()<<endl;
		}
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}
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

		if(! meta.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}
	}
	
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
		for(UnsignedInt i=0; i< meta.size(); ++i)
		{ 
			visualize_(meta[i], item);
		}
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
			
		//check for metainfo objects
		if(! meta.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(meta), item);
		}  
			
		//Check for source file
		visualize_(meta.getSourceFile(), item);
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
	}
	
}






