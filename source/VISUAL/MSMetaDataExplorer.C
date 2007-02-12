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
#include <OpenMS/config.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
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

//QT
#include <qlistview.h>
#include <qwidget.h>
#include <iostream>
#include <sstream>
#include <qpushbutton.h>


using namespace std;
using namespace OpenMS;

//MSMetaDataExplorer::MSMetaDataExplorer(QWidget *parent, const char *name, bool modal, WFlags fl): QDialog(parent, name, modal, fl)
MSMetaDataExplorer::MSMetaDataExplorer(bool editable, QWidget *parent, const char *name, bool modal, WFlags fl)
: QDialog(parent, name, modal, fl), editable_(editable)
{    
  
	//basiclayout_ = new QHBoxLayout(this);
	basiclayout_ = new QGridLayout(this, 3, 2);
	basiclayout_->setSpacing(6);
  basiclayout_->setMargin(11); 
	
	vertlayout_ = new QVBoxLayout();
	vertlayout_->setSpacing(6);
  vertlayout_->setMargin(11); 
	
	buttonlayout_ = new QHBoxLayout();
	buttonlayout_->setSpacing(6);
	buttonlayout_->setMargin(11); 
	
	//layer_bar_ = new QToolBar(this,"Browse in Metadata tree");
  		
	//Create the tree for exploring data 
	listview_ = new QListView(this);    
	listview_->addColumn(QObject::tr("Browse in Metadata tree"));
	listview_->addColumn(QObject::tr(""));
	listview_->setColumnWidth(1,0);
	listview_->setResizeMode(QListView::NoColumn);
	listview_->setRootIsDecorated(true);
  listview_->setSorting(2, true);

	
	//Vertical splitter for the right part
	//splitvert_ = new QSplitter(Vertical, this); 
  
	basiclayout_->addMultiCellWidget(listview_, 0,3, 0,0);
	//basiclayout_->addWidget(splitvert_);
	basiclayout_->addMultiCellLayout(vertlayout_, 0,3, 1,1);
	  
	//Create WidgetStack for managing visible metadata
  //ws_ = new QWidgetStack(splitvert_);
	ws_ = new QWidgetStack(this);
	//splitvert_ = new QWidgetStack(this);
	QLabel* hline = new QLabel(this);
	hline->setFrameShape(QFrame::HLine); 
	
	vertlayout_->addWidget(ws_);
	//vertlayout_->addWidget(splitvert_);
	vertlayout_->addWidget(hline);
	vertlayout_->addLayout(buttonlayout_);
	
	/*	
	QWidget* wid = new QWidget(splitvert_);
	wid->setFixedHeight(70);
	
	glayout_ = new QGridLayout(wid);
	glayout_->setSpacing(6);
  glayout_->setMargin(11);
	
	vertlayout_ = new QHBoxLayout();
	vertlayout_->setSpacing(6);
  vertlayout_->setMargin(11); 
	
	glayout_->addLayout(vertlayout_, 0,0);
	
	vertlayout_->addStretch(1);
	*/
	
	
	
	if(isEditable())
	{
	
		saveallbutton_ = new QPushButton("OK", this);
	  cancelbutton_ = new QPushButton("Cancel", this);
		buttonlayout_->addStretch(1);
		buttonlayout_->addWidget(saveallbutton_ );
		buttonlayout_->addWidget(cancelbutton_ );
		//closebutton_->hide();
		connect(saveallbutton_, SIGNAL(clicked()), this, SLOT(saveAll_())  );
		connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject())  );
	}
	else 
	{
	
		closebutton_ = new QPushButton("Close", this);
		buttonlayout_->addStretch(1);
		buttonlayout_->addWidget(closebutton_ );
		connect(closebutton_, SIGNAL(clicked()), this, SLOT(reject())  );
		//cancelbutton_->hide();
		//saveallbutton_->hide();
	}
	
  connect(listview_, SIGNAL(clicked(QListViewItem*)), this, SLOT(showDetails_(QListViewItem*))  );
	
  
	//create and initialize id number for objects, so that the widgets may be identified
	//used in function makeID
  obj_id_ = 0;
	
  
}//end of constructor



bool MSMetaDataExplorer::isEditable()
{
		return editable_;
}


void MSMetaDataExplorer::showDetails_(QListViewItem *item)
{
  //Signal clicked from the ListItem arrives here.
  //This SLOT emits further signals, so that the atrributes
  //list on the right side will be updated automatically.
  
	if ( !item )
	{
    return;
  }
		
  String s((const char*) item->text(1));
	int id= s.toInt();
  ws_->raiseWidget(id);
  
}

//Save all changes
void MSMetaDataExplorer::saveAll_()
{
	try
	{
		BaseVisualizer* b_ptr;
			
		//call internal store function of all active visualizer objects
		for(int i=0; i < obj_id_; ++i)
		{
			QWidget* widget = ws_->widget(i);
			b_ptr = dynamic_cast<BaseVisualizer*>(widget);
			b_ptr->store();
		}
		
		//close dialog
		this->accept();
		
	}
	catch(exception& e)
	{
		cout<<"Exception: "<<e.what()<<endl;
	}
	
	
}

int MSMetaDataExplorer::makeID_()
{
  int id= obj_id_;
  obj_id_++;
  
  return id;
}

//--------------------------------------------------------------------------
//	Functions modifying the listview
//--------------------------------------------------------------------------

void MSMetaDataExplorer::updateProteinHits_(ProteinIdentification pid, int tree_item_id)
{
		
		ProteinIdentification protid = pid;
		float threshold = protid.getProteinSignificanceThreshold();
		
		// find item in tree belonging to ProteinIdentification object
		const QString qs = String(tree_item_id);
		QListViewItem *item =listview_->findItem(qs, 1);
		
			
		item->sortChildItems(2, true);
		//item->sort();		
		// get all children of this item
		QListViewItem *first =item->firstChild();
		QListViewItem *current = first;
		
		//set the items visible or not visible depending to their score and the current threshold
		for(int i=0; i<item->childCount(); ++i)
		{
							
				if( (current->text(2)).toFloat() <=threshold)
				{
						current->setVisible(false);
				}
				else
				{
					current->setVisible(true);
				}
				
				QListViewItem *next= current->nextSibling();
				current = next;
		}
		
}


void MSMetaDataExplorer::updateNonRefPeptideHits_(Identification id, int tree_item_id)
{
		
		Identification ident = id;
		vector< ProteinHit > prots= id.getProteinHits(); 
		
		String date_time;
		multimap< String, ProteinHit > map;
		vector<PeptideHit> *hits;
		
		ident.getDateTime().get(date_time);
		
		for(vector<ProteinHit>::iterator it = prots.begin(); it != prots.end();	it++)
		{
			map.insert(make_pair(date_time, *it));
		}
		
		hits = ident.getNonReferencingHits(map);
		
				
		// find item in tree belonging to Identification object
		const QString qs = String(tree_item_id);
		QListViewItem *item =listview_->findItem(qs, 1);
		
		// get all children of this item
		QListViewItem *first =item->firstChild();
		QListViewItem *current = first;
		
			
		//set the items visible or not visible depending of non referencing peptide hits
		for(int i=0; i<item->childCount(); ++i)
		{
				for(UnsignedInt i=0; i<hits->size(); ++i)
				{
					String name("Pep ");
				  String score((*hits)[i].getScore());
					String seq = (*hits)[i].getSequence();
		
					name = name+ seq+" (" +score+")";
					
					
						if( current->text(0)==name && current->text(3)== seq)
						{
								current->setVisible(true);
								break;
						}
						else
						{
							current->setVisible(false);
						}
				}
				QListViewItem *next= current->nextSibling();
				current = next;	
		}

}

void MSMetaDataExplorer::updatePeptideHits_(Identification id, int tree_item_id)
{
		
		Identification& ident = id;
		float threshold = ident.getPeptideSignificanceThreshold();
				
		// find item in tree belonging to Identification object
		const QString qs = String(tree_item_id);
		QListViewItem *item =listview_->findItem(qs, 1);
		
		// get all children of this item
		QListViewItem *first =item->firstChild();
		QListViewItem *current = first;
		
		//set the items visible or not visible depending to their score and the current threshold
		for(int i=0; i<item->childCount(); ++i)
		{
							
				if( (current->text(2)).toFloat() <=threshold)
				{
						current->setVisible(false);
				}
				else
				{
					current->setVisible(true);
				}
				
				QListViewItem *next= current->nextSibling();
				current = next;
		}
		
}

void MSMetaDataExplorer::updateRefPeptideHits_(Identification id, int tree_item_id, String ref_date, String ref_acc)
{
	  if(ref_date.trim()=="" && ref_acc.trim()=="")
		{
			return;
		}
		
		
		Identification& ident = id;
				
		// find item in tree belonging to Identification object
		const QString qs = String(tree_item_id);
		QListViewItem *item =listview_->findItem(qs, 1);
		
		// get all children of this item
		QListViewItem *first =item->firstChild();
		QListViewItem *current = first;
		bool hit = false;
				
		ident.sort();
		
		vector< PeptideHit >& peps= ident.getPeptideHits(); 
		
		//Set all items to not visible
			for(int i=0; i<item->childCount(); ++i)
			{
				 current->setVisible(false);
				 QListViewItem *next= current->nextSibling();
				 current = next;
			}
			current = first;
		
		//search all peptide hits
		for(UnsignedInt i=0; i<peps.size(); ++i)
		{
			vector<pair<String, String> >& protlist = peps[i].getProteinIndices();
			
			String name("Pep ");
			String score(peps[i].getScore());
			String seq = peps[i].getSequence();
			name = name+ seq+ " (" +score+")";
			
			
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
						if( current->text(0)==name && current->text(3)== seq )
						{
							current->setVisible(true);
						}					
						QListViewItem *next= current->nextSibling();
						current = next;
				}
				current = first;	
			}
			
			hit = false;
			
			
		}
			
		
}


//-------------------------------------------------------------------------------
//	overloaded visualize functions to call the corresponding data visualizer
//-------------------------------------------------------------------------------

//Visualizing sample object
void MSMetaDataExplorer::visualize_(Sample &s, QListViewItem* parent )
{
    				
		QListViewItem* item;
    int widgetID = makeID_();
    String w_id(widgetID);
		String name(s.getName());
		String id("Sample ");
		id = id+name;
		
		if(parent == 0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), id , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, id , w_id, "0" );
		}
		
		OpenMS::SampleVisualizer *sv = new OpenMS::SampleVisualizer(isEditable(), this); 
    ws_->addWidget(sv, widgetID);
    sv->load(s);

				
		//check for treatments
		if(s.countTreatments() != 0)
		{	
			for(SignedInt i=0; i<s.countTreatments(); ++i)
			{
				if(s.getTreatment(i).getType()=="Digestion")
				{
					visualize_((const_cast<Digestion&>(dynamic_cast<const Digestion&>(s.getTreatment(i)) ) ), item );
				}
				else if(s.getTreatment(i).getType()== "Modification")
				{
					//Cast SampleTreatment reference to a const modification reference
					visualize_((const_cast<Modification&>(dynamic_cast<const Modification&>(s.getTreatment(i))) ), item );
										
				}
				else if(s.getTreatment(i).getType()=="Tagging")
				{
					visualize_((const_cast<Tagging&>(dynamic_cast<const Tagging&>(s.getTreatment(i)) )), item );
				}
				
			}
		}
		
		//Check for subsamples
		vector<Sample>& v= s.getSubsamples();  
		
		if(v.size() != 0)
		{
			for(UnsignedInt i=0; i<v.size(); ++i)    
			{
				visualize_(v[i], item);
			}
		}
		
		//check for metainfo objects
		if(! s.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(s), item);
		}
		
		
}

//Visualizing MetaInfoInterface object
void MSMetaDataExplorer::visualize_(MetaInfoInterface &m, QListViewItem* parent)
{
    
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			new QListViewItem(listview_, listview_->lastItem(), "MetaInfo" , w_id, "0" );
		}
		else
		{
			new QListViewItem(parent, "MetaInfo" , w_id, "0" );
		}
		OpenMS::MetaInfoVisualizer *meta = new OpenMS::MetaInfoVisualizer(isEditable(),this);  
		ws_->addWidget(meta, widgetID);
		meta->load(m);   
		
}



//Visualizing Digestion object
void MSMetaDataExplorer::visualize_(Digestion &d, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "Digestion" , w_id, "0" );
			 	
		}
		else
		{
			item = new QListViewItem(parent, "Digestion" , w_id, "0" );
		}
		
		//check for metainfo objects
		if(! d.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(d), item);
		}
		OpenMS::DigestionVisualizer *dig = new OpenMS::DigestionVisualizer(isEditable(), this);  
		ws_->addWidget(dig, widgetID);
		dig->load(d);   
}



//Visualizing modification object
void MSMetaDataExplorer::visualize_(Modification &m, QListViewItem* parent)
{
		QListViewItem* item;
    int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "Modification" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "Modification" , w_id, "0" );
		}
		//check for metainfo objects
		if(! m.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(m), item);
		}
		OpenMS::ModificationVisualizer *mod = new OpenMS::ModificationVisualizer(isEditable(), this);  
		ws_->addWidget(mod, widgetID);
		mod->load(m);   
}


//Visualizing tagging object
void MSMetaDataExplorer::visualize_(Tagging &t, QListViewItem* parent)
{
		QListViewItem* item;
    int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "Tagging" , w_id, "0" );
		}
		else
		{
			item =new QListViewItem(parent, "Tagging" , w_id, "0" );  
		}
		//check for metainfo objects
		if(! t.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(t), item);
		}
		OpenMS::TaggingVisualizer *tag = new OpenMS::TaggingVisualizer(isEditable(), this);  
		ws_->addWidget(tag, widgetID);
		tag->load(t);  
}



//Visualizing Software object
void MSMetaDataExplorer::visualize_(Software &s, QListViewItem* parent)
{
    QListViewItem* item;
    int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "Software" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "Software" , w_id, "0" );  
		}
		
		//check for gradient object
		OpenMS::SoftwareVisualizer *software = new OpenMS::SoftwareVisualizer(isEditable(), this);  
		ws_->addWidget(software, widgetID);
		software->load(s);  
}



//Visualizing HPLC object
void MSMetaDataExplorer::visualize_(HPLC &h, QListViewItem* parent)
{
    QListViewItem* item;
    int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "HPLC" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "HPLC" , w_id, "0" );  
		}
		
		//check for gradient object
		OpenMS::HPLCVisualizer *hplc = new OpenMS::HPLCVisualizer(isEditable(), this);  
		ws_->addWidget(hplc, widgetID);
		hplc->load(h);  
			
		try
		{	
		
			visualize_(h.getGradient(), item);
		
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize Gradient. "<<e.what()<<endl;
		}
}

//Visualizing Gradient object
void MSMetaDataExplorer::visualize_(Gradient &gradient, QListViewItem* parent)
{
    int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			new QListViewItem(listview_, listview_->lastItem(), "Gradient" , w_id, "0" );
		}
		else
		{
			new QListViewItem(parent, "Gradient" , w_id, "0" );  
		}
		
		//check for gradient object
		OpenMS::GradientVisualizer *grad = new OpenMS::GradientVisualizer(isEditable(), this); 
		ws_->addWidget(grad, widgetID);
		grad->load(gradient);  
		
		
        
}

//Visualizing SourceFile object
void MSMetaDataExplorer::visualize_(SourceFile &source, QListViewItem* parent)
{
    int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			new QListViewItem(listview_, listview_->lastItem(), "SourceFile" , w_id, "0" );
		}
		else
		{
			new QListViewItem(parent, "SourceFile" , w_id, "0" );  
		}
		
		OpenMS::SourceFileVisualizer *sourcefile = new OpenMS::SourceFileVisualizer(isEditable(), this); 
		ws_->addWidget(sourcefile, widgetID);
		sourcefile->load(source);  
		
		
        
}

//Visualizing ContactPerson object
void MSMetaDataExplorer::visualize_(ContactPerson &person, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "ContactPerson" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "ContactPerson" , w_id, "0" );  
		}
		
		OpenMS::ContactPersonVisualizer *cp = new OpenMS::ContactPersonVisualizer(isEditable(), this); 
		ws_->addWidget(cp, widgetID);
		cp->load(person);    
		
		//check for metainfo objects
		if(! person.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(person), item);
		}
		
}

//Visualizing Instrument object
void MSMetaDataExplorer::visualize_(Instrument &instrument, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "Instrument" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "Instrument" , w_id, "0" );  
		}
		
		OpenMS::InstrumentVisualizer *iv = new OpenMS::InstrumentVisualizer(isEditable(), this); 
		ws_->addWidget(iv, widgetID);
		iv->load(instrument);    
		
		//check for metainfo objects
		if(! instrument.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(instrument), item);
		}
		
		
		//visualize IonSource object
		try
		{	
		
			visualize_(instrument.getIonSource(), item);
		
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize IonSource. "<<e.what()<<endl;
		}
		
		//visualize IonDetector object
		try
		{	
		
			visualize_(instrument.getIonDetector(), item);
		
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize IonDetector. "<<e.what()<<endl;
		}
		
		
		//Check for MassAnalyzers
		vector<MassAnalyzer>& v= instrument.getMassAnalyzers();  
		
		if(v.size() != 0)
		{
			for(UnsignedInt i=0; i<v.size(); ++i)    
			{
				visualize_(v[i], item);
			}
		}
		
		
		
}

//Visualizing IonSource object
void MSMetaDataExplorer::visualize_(IonSource &is, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "IonSource" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "IonSource" , w_id, "0" );  
		}
		
		OpenMS::IonSourceVisualizer *isv = new OpenMS::IonSourceVisualizer(isEditable(), this); 
		ws_->addWidget(isv, widgetID);
		isv->load(is);    
		
		//check for metainfo objects
		if(! is.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(is), item);
		}
		
}

//Visualizing IonDetector object
void MSMetaDataExplorer::visualize_(IonDetector &id, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "IonDetector" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "IonDetector" , w_id, "0" );  
		}
		
		OpenMS::IonDetectorVisualizer *idv = new OpenMS::IonDetectorVisualizer(isEditable(), this); 
		ws_->addWidget(idv, widgetID);
		idv->load(id);    
		
		//check for metainfo objects
		if(! id.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(id), item);
		}
		
}

//Visualizing MassAnalyzer object
void MSMetaDataExplorer::visualize_(MassAnalyzer &ma, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "MassAnalyzer" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "MassAnalyzer" , w_id, "0" );  
		}
		
		OpenMS::MassAnalyzerVisualizer *mav = new OpenMS::MassAnalyzerVisualizer(isEditable(), this); 
		ws_->addWidget(mav, widgetID);
		mav->load(ma);    
		
		//check for metainfo objects
		if(! ma.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(ma), item);
		}
		
}

//Visualizing ProcessingMethod object
void MSMetaDataExplorer::visualize_(ProcessingMethod &pm, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "ProcessingMethod" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "ProcessingMethod" , w_id, "0" );  
		}
		
		OpenMS::ProcessingMethodVisualizer *pmv = new OpenMS::ProcessingMethodVisualizer(isEditable(), this); 
		ws_->addWidget(pmv, widgetID);
		pmv->load(pm);    
		
		//check for metainfo objects
		if(! pm.isMetaEmpty() )
		{
			visualize_(dynamic_cast<MetaInfoInterface&>(pm), item);
		}
		
}

//Visualizing ProteinIdentification object
void MSMetaDataExplorer::visualize_(ProteinIdentification &pid, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "ProteinIdentification" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "ProteinIdentification" , w_id, "0" );  
		}
		
		OpenMS::ProteinIdentificationVisualizer *pidv = new OpenMS::ProteinIdentificationVisualizer(splitvert_, this); 
		ws_->addWidget(pidv, widgetID);
		pidv->load(pid, widgetID);    
		
		//check for proteinhits objects
		pid.sort();
		vector< ProteinHit > v= pid.getProteinHits();  
			
		//item->sortChildItems(2, true);
		
				
		//list all protein hits in the tree
		if(v.size() != 0)
		{
			for(UnsignedInt i=0; i<v.size(); ++i)    
			{
				visualize_(v[i], item);
			}
		}
		//item->sortChildItems(2, true);
		
		
			
}

//Visualizing Identification object
void MSMetaDataExplorer::visualize_(Identification &id, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "Identification" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "Identification" , w_id, "0" );  
		}
		
		OpenMS::IdentificationVisualizer *idv = new OpenMS::IdentificationVisualizer(splitvert_, this); 
		ws_->addWidget(idv, widgetID);
		idv->load(id, widgetID);    
		
		//check for proteins and peptides hits
		id.sort();
		vector< ProteinHit > prots= id.getProteinHits(); 
		vector< ProteinHit > *prots_ptr = &prots; 
		vector< PeptideHit > peps= id.getPeptideHits();  
			
		//list all peptides hits in the tree
		if(peps.size() != 0)
		{
			for(UnsignedInt i=0; i<peps.size(); ++i)    
			{	
					
				visualize_(peps[i], item, prots_ptr);
			}
		}
		
		
			
}



//Visualizing ProteinHit object
void MSMetaDataExplorer::visualize_(ProteinHit &phit, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		String name("Prot ");
		String identifier(phit.getAccession());
		String score(phit.getScore());
		name = name+ identifier+ "  (" +score+")";
		
		if(parent==0)
		{
			//item = new QListViewItem(listview_, listview_->lastItem(), "ProteinHit" , w_id, String(phit.getScore())   );
			item = new QListViewItem(listview_, listview_->lastItem(), name , w_id, String(phit.getScore())   );
		}
		else
		{
			//item = new QListViewItem(parent, "ProteinHit" , w_id, String(phit.getScore()));  
			item = new QListViewItem(parent, name , w_id, String(phit.getScore()));  
		}
		
		OpenMS::ProteinHitVisualizer *phitv = new OpenMS::ProteinHitVisualizer(isEditable(), this); 
		ws_->addWidget(phitv, widgetID);
		phitv->load(phit);    
}

//Visualizing PeptideHit object
void MSMetaDataExplorer::visualize_(PeptideHit &pephit, QListViewItem* parent, vector<ProteinHit> *prots )
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		
		String name("Pep ");
		String score(pephit.getScore());
		String seq(pephit.getSequence());
		name = name+ seq+ " (" +score+")";
		
		if(parent==0)
		{
				item = new QListViewItem(listview_, listview_->lastItem(), name , w_id, String(pephit.getScore()), pephit.getSequence()   );
		}
		else
		{
				//(parent, Col0=name, Col1=w_id ,Col2=score)
				item = new QListViewItem(parent, name , w_id, String(pephit.getScore()), pephit.getSequence());  
		}
		
		OpenMS::PeptideHitVisualizer *pephitv = new OpenMS::PeptideHitVisualizer(isEditable(), this); 
		ws_->addWidget(pephitv, widgetID);
		pephitv->load(pephit); 
		
		//sort all ProteinHits by name in ascending order
		//item->sortChildItems(2, true);
		
		//get all protein hits for the peptide hits
		vector<pair<String, String> > protlist = pephit.getProteinIndices();
		
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
void MSMetaDataExplorer::visualize_(ExperimentalSettings &es, QListViewItem* parent)
{
		
		QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "ExperimentalSettings" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "ExperimentalSettings" , w_id, "0" );
		}
		
		//OpenMS::ExperimentalSettingsVisualizer *esv = new OpenMS::ExperimentalSettingsVisualizer(isEditable(), splitvert_);  
		OpenMS::ExperimentalSettingsVisualizer *esv = new OpenMS::ExperimentalSettingsVisualizer(isEditable(), this);  
		ws_->addWidget(esv, widgetID);
		esv->load(es);
		
		
		//check for metainfo objects
		if(! es.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(es), item);
		}
		
		try
		{
			//check for Sample
			Sample& sample = es.getSample();
			visualize_(sample, item);
			
			
			//check for ProteinIdentification
			vector<ProteinIdentification>& protIDs= es.getProteinIdentifications();  
			
			if(protIDs.size() != 0)
			{
				for(UnsignedInt i=0; i<protIDs.size(); ++i)    
				{
					visualize_(protIDs[i], item);
				}
			}
			
			//check for ProcessingMethod
			ProcessingMethod& pm = es.getProcessingMethod();
			visualize_(pm, item);
			
			//check for Instrument
			Instrument& instrument = es.getInstrument();
			visualize_(instrument, item);
			
			//check for SourceFile
			SourceFile& sf = es.getSourceFile();
			visualize_(sf, item);
			
			
			//check for ContactPersons
			vector<ContactPerson>& cps= es.getContacts();  
			
			if(cps.size() != 0)
			{
				for(UnsignedInt i=0; i<cps.size(); ++i)    
				{
					visualize_(cps[i], item);
				}
			}  
			
			
			//check for Software
			Software& sw = es.getSoftware();
			visualize_(sw, item);
			
			//check for HPLC
			HPLC& hplc = es.getHPLC();
			visualize_(hplc, item);
			
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
void MSMetaDataExplorer::visualize_(SpectrumSettings &ss, QListViewItem* parent)
{
		QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "SpectrumSettings" , w_id, "0" );
		}
		else
		{
			item = new QListViewItem(parent, "SpectrumSettings" , w_id, "0" );
		}
				
		OpenMS::SpectrumSettingsVisualizer *ssv = new OpenMS::SpectrumSettingsVisualizer(isEditable(), this);  
		ws_->addWidget(ssv, widgetID);
		ssv->load(ss);
		
		
		
		try
		{
			
			//check for InstrumentSettings
			InstrumentSettings& insSet = ss.getInstrumentSettings();
			visualize_(insSet, item);
			
			//check for Identification
			vector<Identification>& IDs= ss.getIdentifications();  
			if(IDs.size() != 0)
			{
				for(UnsignedInt i=0; i<IDs.size(); ++i)    
				{
					visualize_(IDs[i], item);
				}
			}
			
			//check for Precursor
			Precursor& pre = ss.getPrecursor();
			visualize_(pre, item);
					
			//check for MetaInfoDescription
			std::map<String, MetaInfoDescription>& mid =ss.getMetaInfoDescriptions();
			//mid=ss.getMetaInfoDescriptions();
			String key("");
			String *key_ptr=&key;
			
			map<String, MetaInfoDescription>::iterator iter;   
  		for( iter = mid.begin(); iter != mid.end(); iter++ ) 
			{
				key=iter->first;
				visualize_(iter->second, item, key_ptr );
      }
			
			 					
			
			//check for AcquisitionInfo
			AcquisitionInfo& ac = ss.getAcquisitionInfo();
			visualize_(ac, item);
		
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to visualize SpectrumSettings. "<<e.what()<<endl;
		}
		
}

//Visualizing InstrumentSettings object
void MSMetaDataExplorer::visualize_(InstrumentSettings &is, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "InstrumentSettings" , w_id, "0" );
			 	
		}
		else
		{
			item = new QListViewItem(parent, "InstrumentSettings" , w_id, "0" );
		}
		
		//check for metainfo objects
		if(! is.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(is), item);
		}
		OpenMS::InstrumentSettingsVisualizer *isv = new OpenMS::InstrumentSettingsVisualizer(isEditable(), this);  
		ws_->addWidget(isv, widgetID);
		isv->load(is);   
}

//Visualizing Acquisition object
void MSMetaDataExplorer::visualize_(Acquisition &a, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		String name("Acquisition Nr. ");
		String number(a.getNumber());
		name = name+ number;
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), name , w_id, "0" );
			 	
		}
		else
		{
			item = new QListViewItem(parent, name , w_id, "0" );
		}
		
		//check for metainfo objects
		if(! a.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(a), item);
		}
		OpenMS::AcquisitionVisualizer *av = new OpenMS::AcquisitionVisualizer(isEditable(), this);  
		ws_->addWidget(av, widgetID);
		av->load(a);   
}

//Visualizing AcquisitionInfo object
void MSMetaDataExplorer::visualize_(AcquisitionInfo &ai, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "Acquisition Info " , w_id, "0" );
			 	
		}
		else
		{
			item = new QListViewItem(parent, "Acquisition Info " , w_id, "0" );
		}
		
		
		OpenMS::AcquisitionInfoVisualizer *aiv = new OpenMS::AcquisitionInfoVisualizer(isEditable(), this);  
		ws_->addWidget(aiv, widgetID);
		aiv->load(ai); 
		
	//Get Aquisition objects
	for(UnsignedInt i=0; i< ai.size(); ++i)
	{ 
		Acquisition &a = ai[i]; 
	  visualize_(a, item);
	}  
	
}

//Visualizing MetaInfoDescription object
void MSMetaDataExplorer::visualize_(MetaInfoDescription &mid,  QListViewItem* parent, String* key)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		String name("MetaInfoDescription ");
		name = name + (*key);
		
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), name , w_id, "0" );
			 	
		}
		else
		{
			item = new QListViewItem(parent, name , w_id, "0" );
		}
		
		//check for metainfo objects
		if(! mid.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(mid), item);
		}
		OpenMS::MetaInfoDescriptionVisualizer *midv = new OpenMS::MetaInfoDescriptionVisualizer(isEditable(), this);  
		ws_->addWidget(midv, widgetID);
		midv->load(mid);   
		
		//Check for source file
		SourceFile& source= mid.getSourceFile();
		visualize_(source, item);
}


//Visualizing Precursor object
void MSMetaDataExplorer::visualize_(Precursor &pre, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID_();
		String w_id(widgetID);
		
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "Precursor " , w_id, "0" );
			 	
		}
		else
		{
			item = new QListViewItem(parent, "Precursor " , w_id, "0" );
		}
		
		
		OpenMS::PrecursorVisualizer *prev = new OpenMS::PrecursorVisualizer(isEditable(), this);  
		ws_->addWidget(prev, widgetID);
		prev->load(pre);   
}








