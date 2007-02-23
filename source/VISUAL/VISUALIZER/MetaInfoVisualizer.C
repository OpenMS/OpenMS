// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  this library is distributed in the hope that it will be useful,
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
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoVisualizer.h>

//QT
#include <QtGui/QGridLayout>
#include <QtGui/QPushButton>
#include <QtGui/QLineEdit>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QButtonGroup>

//STL
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

namespace OpenMS
{

//Constructor
MetaInfoVisualizer::MetaInfoVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
  type_="MetaInfo";
	
	buttongroup_ = new QButtonGroup();
	nextrow_=0;
	
	viewlayout_ = new QGridLayout(this);
	viewlayout_->setSpacing(6);
  viewlayout_->setMargin(11);
	addLabel("Modify MetaData information.");
	addSeperator();	
	mainlayout_->addLayout(viewlayout_, row_, row_, 0 ,2);
	//increase row counter for mainlayout_.
	row_++;
		
}


void MetaInfoVisualizer::load(MetaInfoInterface &m)
{	
	
	ptr_ = &m;
	
	//Copy of current object for restoring the original values
	tempmeta_=m;
	  
		
	//keys_ is a vector of indices of type UnsignedInt
	tempmeta_.getKeys(keys_);
		
		
	//Load actual metaInfo Data into viewLayout_
	for(UnsignedInt i=0; i< keys_.size(); ++i)
	{ 
	  loadData_(keys_[i]);
	}
	
	addSeperator();
	addLabel("Add new MetaInfo entry.");
	addLineEdit(newkey_, "Key");
	addLineEdit(newdescription_, "Description");
	addLineEdit(newvalue_, "Value");
	addHorizontalButtons(addbutton_, "Add", clearbutton_, "Clear");
	if(!isEditable())
	{
				addbutton_->setEnabled(false);
				clearbutton_->setEnabled(false);
	}
	addVSpacer();
	
	finishAdding_();
		
	connect(buttongroup_, SIGNAL(clicked(int)), this, SLOT(remove(int)));
	connect(addbutton_, SIGNAL(clicked()), this, SLOT(add()));
	connect(clearbutton_, SIGNAL(clicked()), this, SLOT(clear()));
  
	
}




//----------------------------------------------------------------------------
//			SLOTS
//----------------------------------------------------------------------------

void MetaInfoVisualizer::remove(int index)
{	
  UnsignedInt id=(UnsignedInt)index;
	
	
	//Remove label
	std::vector<std::pair<UnsignedInt,QLabel*> >::iterator iter;
	for(iter = metalabels_.begin(); iter < metalabels_.end(); iter++ ) 
	{	 
     if( (*iter).first == id)
		 {
				viewlayout_->removeWidget((*iter).second);
				delete (*iter).second;
				(*iter).second =0;
				metalabels_.erase(iter);
		 }
   }
	 
	//Remove QLineEdit	
	std::vector<std::pair<UnsignedInt,QLineEdit*> >::iterator iter2; 
	for(iter2 = metainfoptr_.begin(); iter2 < metainfoptr_.end(); iter2++ ) 
	{
     if( (*iter2).first == id)
		 {
				viewlayout_->removeWidget((*iter2).second);
				delete (*iter2).second;
				(*iter2).second =0;
				metainfoptr_.erase(iter2);
		}	
		
		
   }
	 
	//Remove QButton  
	std::vector<std::pair<UnsignedInt,QAbstractButton*> >::iterator iter3; 
	for(iter3 = metabuttons_.begin(); iter3 < metabuttons_.end(); iter3++ ) 
	{
     if( (*iter3).first == id)
		 {	
				viewlayout_->removeWidget((*iter3).second);
				delete (*iter3).second;
				(*iter3).second =0;
				metabuttons_.erase(iter3);
		 }
   }	
		
  //Remove entry from metainfointerface
	tempmeta_.removeMetaValue(id);
	tempmeta_.getKeys(keys_);
		
  this->repaint();
}


void MetaInfoVisualizer::loadData_(UnsignedInt index)
{
  //----------------------------------------------------------------------------	
  //  All metainfo goes into the viewlayout_ 
  //----------------------------------------------------------------------------
	
			QLabel* lab;
			QLineEdit* ptr;
			QPushButton* button;
			
			
			lab = new QLabel(tempmeta_.metaRegistry().getName(index).c_str(), this);
			viewlayout_->addWidget(lab, nextrow_, 0);
			
			ptr = new QLineEdit(this);
			ptr->setText( tempmeta_.getMetaValue(index).toString().c_str() );
			viewlayout_->addWidget(ptr, nextrow_, 1);
							
			button = new QPushButton("Remove", this);
			if(!isEditable())
			{
				button->setEnabled(false);
			}
			viewlayout_->addWidget(button, nextrow_, 2);
			
			//Store information about ID(index) and QWidget
			metalabels_.push_back(make_pair(index,lab));	
			metainfoptr_.push_back(make_pair(index,ptr));	
			metabuttons_.push_back(make_pair(index,button));
				
			//Insert new button with ID into buttongroup
			buttongroup_->addButton(button, index);
			nextrow_++;
			
			lab->show();
			ptr->show();
			button->show();
	
			
}


void MetaInfoVisualizer::add()
{ 
 
	String name(newkey_->text().toStdString()) ;
	String description(newdescription_->text().toStdString());
	String value(newvalue_->text().toStdString());
	
	
	if(name.trim().length() ==0 )
	{ 
		//Must have a name
		return;
	}
	
	//Register new entry and update metainfointerface object
	UnsignedInt newindex = tempmeta_.metaRegistry().registerName(name, description, "");
	
	//Store new data in temporary metainfo object
	tempmeta_.setMetaValue(newindex, value);
	tempmeta_.getKeys(keys_);
	
		
  //Update viewlayout_ 
	try
	{
		//check whether there is already an entry in GUI for added metainfo.
		//If index already exists, return and do nothing. 
		if( buttongroup_->button(newindex) != 0 )
		{
		  return;
		}
		
		loadData_(newindex);	
			
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to create new meta info data. "<<e.what()<<endl;
	} 
}

void MetaInfoVisualizer::clear()
{
		
	newkey_->clear();
	newdescription_->clear();
	newvalue_->clear();
	
}



void MetaInfoVisualizer::store()
{
	try
	{	
		
		//Store QLineEdit information
		std::vector<std::pair<UnsignedInt,QLineEdit*> >::iterator iter2; 
		for(iter2 = metainfoptr_.begin(); iter2 < metainfoptr_.end(); iter2++ ) 
		{
			UnsignedInt index = (*iter2).first;
			String value(((*iter2).second)->text().toStdString());
			tempmeta_.setMetaValue(index, value);
		}
		//copy temporary stored data into metainfo object
		(*ptr_)=tempmeta_;
				
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new metainfo data. "<<e.what()<<endl;
	}
}

void MetaInfoVisualizer::reject()
{
	try
	{
		//Delete all data in GUI		
		//Need a temporary container, because function remove(int index) modifies container keys_.
		std::vector<UnsignedInt> keys_temp= keys_;
		
		for(UnsignedInt i =0; i< keys_temp.size(); ++i)
		{	
					remove(keys_temp[i]);
		}
		
		
		metalabels_.clear();
		metainfoptr_.clear();
		metabuttons_.clear();
		
		//Restore original data and refill GUI
		tempmeta_=(*ptr_);
		nextrow_=0;
		keys_.clear();
		ptr_->getKeys(keys_);
		for(UnsignedInt i =0; i< keys_.size(); ++i)
		{	
			loadData_(keys_[i]);
		}
		
		
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original metainfo data. "<<e.what()<<endl;
	} 
}

}
