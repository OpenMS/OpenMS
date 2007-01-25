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
#include <OpenMS/VISUAL/VISUALIZER/AquisitionVisualizer.h>
//#include <OpenMS/VISUAL/VISUALIZER/AquisitionInfoVisualizer.h>

//QT
#include <qlistview.h>
#include <qwidget.h>
#include <iostream>
#include <sstream>
#include <qpushbutton.h>


using namespace std;
using namespace OpenMS;

MSMetaDataExplorer::MSMetaDataExplorer(QWidget *parent, const char *name) : QMainWindow(parent, name)
{
    
  
  
  //Create basic widget  
  split_ = new QSplitter(Horizontal, this); 
  setCentralWidget(split_);
        
  //Create the tree for exploring data 
	listview_ = new QListView(split_);    
	listview_->addColumn(QObject::tr("Object Browser"));
	listview_->setResizeMode(QListView::AllColumns);
  listview_->setRootIsDecorated(true);
  listview_->resize(300,200);
	
	//Vertical splitter for the right part
	splitvert_ = new QSplitter(Vertical, split_); 
  
	  
	//Create WidgetStack for managing visible metadata
  ws_ = new QWidgetStack(splitvert_);
	ws_->resize(300,200);
	
	QWidget* wid = new QWidget(splitvert_);
	wid->setFixedHeight(70);
	
	glayout_ = new QGridLayout(wid);
	glayout_->setSpacing(6);
  glayout_->setMargin(11);
	
	vertlayout_ = new QHBoxLayout();
	vertlayout_->setSpacing(6);
  vertlayout_->setMargin(11); 
	
	glayout_->addLayout(vertlayout_, 0,0);
	
	//saveallbutton_ = new QPushButton("Save All", splitvert_);
	saveallbutton_ = new QPushButton("Save All", wid);
  //saveallbutton_->setFixedSize(120,30);
	vertlayout_->addStretch(1);
	vertlayout_->addWidget(saveallbutton_ );
	
  //show sample attributes
  connect(listview_, SIGNAL(clicked(QListViewItem*)), this, SLOT(showDetails(QListViewItem*))  );
	//save all changes
	connect(saveallbutton_, SIGNAL(clicked()), this, SLOT(saveAll())  );
  
	
  //create and initialize id number for objects, so that the widgets may be identified
	//used in function makeID
  obj_id_ = 0;
  
  //other stuff
  setMinimumSize(800,600);
  
  
}//end of constructor

void MSMetaDataExplorer::showDetails(QListViewItem *item)
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
void MSMetaDataExplorer::saveAll()
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
		
		
	}
	catch(exception& e)
	{
		cout<<"Exception: "<<e.what()<<endl;
	}
	
	
}



//-------------------------------------------------------------------------------
//	overloaded visualize functions to call the corresponding data visualizer
//-------------------------------------------------------------------------------




//Visualizing sample object
void MSMetaDataExplorer::visualize_(Sample &s, QListViewItem* parent )
{
    				
		QListViewItem* item;
    int widgetID = makeID();
    String w_id(widgetID);
		if(parent == 0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "Sample" , w_id);
		}
		else
		{
			item = new QListViewItem(parent, "Sample" , w_id);
		} 
		OpenMS::SampleVisualizer *sv = new OpenMS::SampleVisualizer(splitvert_); 
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
    
		int widgetID = makeID();
		String w_id(widgetID);
		if(parent==0)
		{
			new QListViewItem(listview_, listview_->lastItem(), "MetaInfo" , w_id);
		}
		else
		{
			new QListViewItem(parent, "MetaInfo" , w_id);
		}
		OpenMS::MetaInfoVisualizer *meta = new OpenMS::MetaInfoVisualizer(splitvert_);  
		ws_->addWidget(meta, widgetID);
		meta->load(m);   
		
}



//Visualizing Digestion object
void MSMetaDataExplorer::visualize_(Digestion &d, QListViewItem* parent)
{
    QListViewItem* item;
		int widgetID = makeID();
		String w_id(widgetID);
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "Digestion" , w_id);
			 	
		}
		else
		{
			item = new QListViewItem(parent, "Digestion" , w_id);
		}
		
		//check for metainfo objects
		if(! d.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(d), item);
		}
		OpenMS::DigestionVisualizer *dig = new OpenMS::DigestionVisualizer(splitvert_);  
		ws_->addWidget(dig, widgetID);
		dig->load(d);   
}



//Visualizing modification object
void MSMetaDataExplorer::visualize_(Modification &m, QListViewItem* parent)
{
		QListViewItem* item;
    int widgetID = makeID();
		String w_id(widgetID);
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "Modification" , w_id);
		}
		else
		{
			item = new QListViewItem(parent, "Modification" , w_id);
		}
		//check for metainfo objects
		if(! m.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(m), item);
		}
		OpenMS::ModificationVisualizer *mod = new OpenMS::ModificationVisualizer(splitvert_);  
		ws_->addWidget(mod, widgetID);
		mod->load(m);   
}


//Visualizing tagging object
void MSMetaDataExplorer::visualize_(Tagging &t, QListViewItem* parent)
{
		QListViewItem* item;
    int widgetID = makeID();
		String w_id(widgetID);
		if(parent==0)
		{
			item =new QListViewItem(listview_, listview_->lastItem(), "Tagging" , w_id);
		}
		else
		{
			item =new QListViewItem(parent, "Tagging" , w_id);  
		}
		//check for metainfo objects
		if(! t.isMetaEmpty() )
		{
				visualize_(dynamic_cast<MetaInfoInterface&>(t), item);
		}
		OpenMS::TaggingVisualizer *tag = new OpenMS::TaggingVisualizer(splitvert_);  
		ws_->addWidget(tag, widgetID);
		tag->load(t);  
}



//Visualizing Software object
void MSMetaDataExplorer::visualize_(Software &s, QListViewItem* parent)
{
    QListViewItem* item;
    int widgetID = makeID();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "Software" , w_id);
		}
		else
		{
			item = new QListViewItem(parent, "Software" , w_id);  
		}
		
		//check for gradient object
		OpenMS::SoftwareVisualizer *software = new OpenMS::SoftwareVisualizer(splitvert_);  
		ws_->addWidget(software, widgetID);
		software->load(s);  
}



//Visualizing HPLC object
void MSMetaDataExplorer::visualize_(HPLC &h, QListViewItem* parent)
{
    QListViewItem* item;
    int widgetID = makeID();
		String w_id(widgetID);
		if(parent==0)
		{
			item = new QListViewItem(listview_, listview_->lastItem(), "HPLC" , w_id);
		}
		else
		{
			item = new QListViewItem(parent, "HPLC" , w_id);  
		}
		
		//check for gradient object
		OpenMS::HPLCVisualizer *hplc = new OpenMS::HPLCVisualizer(splitvert_);  
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
    int widgetID = makeID();
		String w_id(widgetID);
		if(parent==0)
		{
			new QListViewItem(listview_, listview_->lastItem(), "Gradient" , w_id);
		}
		else
		{
			new QListViewItem(parent, "Gradient" , w_id);  
		}
		
		//check for gradient object
		OpenMS::GradientVisualizer *grad = new OpenMS::GradientVisualizer(splitvert_); 
		ws_->addWidget(grad, widgetID);
		grad->load(gradient);  
		
		
        
}

//Visualizing SourceFile object
void MSMetaDataExplorer::visualize_(SourceFile &source, QListViewItem* parent)
{
    int widgetID = makeID();
		String w_id(widgetID);
		if(parent==0)
		{
			new QListViewItem(listview_, listview_->lastItem(), "SourceFile" , w_id);
		}
		else
		{
			new QListViewItem(parent, "SourceFile" , w_id);  
		}
		
		//check for gradient object
		OpenMS::SourceFileVisualizer *sourcefile = new OpenMS::SourceFileVisualizer(splitvert_); 
		ws_->addWidget(sourcefile, widgetID);
		grad->load(source);  
		
		
        
}
int MSMetaDataExplorer::makeID()
{
  int id= obj_id_;
  obj_id_++;
  
  return id;
}







