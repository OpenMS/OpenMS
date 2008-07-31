// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/GradientVisualizer.h>

//QT
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>
#include <QtGui/QPushButton>
#include <QtGui/QGridLayout>

//STL
#include <iostream>
#include <vector>

using namespace std;

namespace OpenMS
{

//Constructor
GradientVisualizer::GradientVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
{
	
	nextrow_=0;
}



void GradientVisualizer::load(Gradient &g)
{
	ptr_ = &g;
	
	//Copy of current object for restoring the original values
	tempgradient_=g;
	 
	addLabel("Modify Gradient information"); 
	addSeperator();
	
	viewlayout_ = new QGridLayout();
	mainlayout_->addLayout(viewlayout_, row_, 0, 1 ,3);
	row_++;		
	//Get the actuall eluent, timepoint and percentage values.
	loadData_();
			
	addSeperator();
	
	addLineEditButton("Add new Eluent", new_eluent_,  add_eluent_button_, "Add Eluent");
	addLineEditButton( "Add new Timepoint", new_timepoint_, add_timepoint_button_, "Add Timepoint");
	addLabel("Attention: All percentage values at a certain timepoint must add up to 100.");
	addSeperator();
	addLabel("Remove all eluents, timepoints and percentage values.");
	addButton(removebutton_, "Remove");
	
	finishAdding_();
	addSeperator();
	connect(add_timepoint_button_, SIGNAL(clicked()), this, SLOT(addTimepoint()) );	
	connect(add_eluent_button_, SIGNAL(clicked()), this, SLOT(addEluent()) );	
  connect(removebutton_, SIGNAL(clicked()), this, SLOT(deleteData()) );	
	
	//Input validator
	timepoint_vali_= new QIntValidator(new_timepoint_);
	new_timepoint_->setValidator(timepoint_vali_);
	
}




//----------------------------------------------------------------------------
//			SLOTS
//----------------------------------------------------------------------------
void GradientVisualizer::deleteData()
{	
	//Remove entries from Gradient
	tempgradient_.clearEluents();
	tempgradient_.clearTimepoints();
	tempgradient_.clearPercentages();
	
	update_();
}


void GradientVisualizer::addTimepoint()
{
  //Check wether new timepoint is in range
	String m(new_timepoint_->text().toStdString()) ;
	int num_time = timepoints_.size();
	
	if(num_time==0 && m.trim().length() !=0)
	{
		tempgradient_.addTimepoint(m.toInt());
		update_();
	}
	else
	{
		if( m.trim().length() !=0 && timepoints_[num_time-1] < m.toInt())
		{
			tempgradient_.addTimepoint(m.toInt());
			update_();
		}
	}
	
}


void GradientVisualizer::addEluent()
{
	String m(new_eluent_->text().toStdString()) ;
	std::vector<String>::iterator iter; 
	//check if eluent name is empty
	if(m.trim().length() !=0 )
	{	
		//check if eluent name already exists
		for(iter = eluents_.begin(); iter < eluents_.end(); iter++)
		{
			if(*iter == m )
			{
				return;
			}
		}
		tempgradient_.addEluent( m );
		update_();
	}
		
}


void GradientVisualizer::store_()
{
	try
	{	
		//Check, whether sum is 100
		int time_size = timepoints_.size();
		int elu_size = eluents_.size();
		int count = 0;
		int elu_count = 0;
		int sum_check=0;
		for(int i=0; i< time_size; ++i )
		{ 
			elu_count=i;
			for(int j=0; j< elu_size; ++j )
			{
				String value((gradientdata_[elu_count])->text().toStdString());
				elu_count  = elu_count+ time_size;
				sum_check= sum_check + value.toInt();
				if(j== elu_size-1 && sum_check!=100)
				{
					cout<<"The sum does not add up to 100 !"<<endl;
					cout<<"Please check your values."<<endl;
					return;
				}
			}
			sum_check = 0;
		}

       
        
		//Store all values into the gradient object
		for(UInt i=0; i< eluents_.size(); ++i )
		{
			for(UInt j=0; j< timepoints_.size(); ++j )
			{
				String value((gradientdata_[count+j])->text().toStdString());
				tempgradient_.setPercentage(eluents_[i], timepoints_[j], value.toInt()  );      
			}
			count=count+time_size;
			
		}

    //copy temporary stored data into metainfo object
		(*ptr_)=tempgradient_;		
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new gradient data. "<<e.what()<<endl;
	}
}

//----------------------------------------------------------------------------
//			private
//----------------------------------------------------------------------------

void GradientVisualizer::loadData_()
{
	
  nextrow_ =0;
	
	eluents_ =    tempgradient_.getEluents();
  timepoints_ = tempgradient_.getTimepoints();
 	UInt num_timepoints = tempgradient_.getTimepoints().size();	
		
		
	//Add labels to display eluent-timepoint-percentage-triplets.	
	QLabel *label = new QLabel("Eluent names\\Timepoints", this);
  viewlayout_->addWidget(label, 0, 0, 1, num_timepoints);
	label->show();
	nextrow_++;
	gradientlabel_.push_back(label);
	
	//----------------------------------------------------------------------
	//			Actual data
	//----------------------------------------------------------------------
	for(UInt i=0; i<timepoints_.size(); ++i )
	{
		//Add labels to display eluent-timepoint-percentage-triplets.	
		QLabel* label1 = new QLabel(String(timepoints_[i]).c_str(), this);
  	viewlayout_->addWidget(label1, 1, i+1);
		label1->show();
		gradientlabel_.push_back(label1);
	}
	
	nextrow_++;
	
	//Add the percentage values for the eluents and their timepoint.	
	for(UInt i=0; i<eluents_.size(); ++i)
	{
	
	QLabel* eluent = new QLabel(eluents_[i].c_str(), this);
  viewlayout_->addWidget(eluent, nextrow_, 0);
	eluent->show();
	gradientlabel_.push_back(eluent);
	
	  for(UInt j=0; j<timepoints_.size(); ++j)
		{
		  percentage_ = new QLineEdit(this);
			percentage_->setText( String(  tempgradient_.getPercentage( eluents_[i], timepoints_[j]) ).c_str()  );
			viewlayout_->addWidget(percentage_, nextrow_, j+1);
			//Store pointers to the QLineEdits
			gradientdata_.push_back(percentage_);	
			percentage_->show();
		}	
						
		nextrow_++;	
	}
		
}

void GradientVisualizer::removeData_()
{	
  
	//Remove QLineEdits
	std::vector<QLineEdit*>::iterator iter2; 
	
	for(iter2 = gradientdata_.begin(); iter2 < gradientdata_.end(); iter2++ ) 
	{			
				//Delete QLineEdit field from viewlayout_
				viewlayout_->removeWidget((*iter2));
				(*iter2)->hide();
				//Set pointer to 0
				(*iter2) =0;			
				//Free memory of the pointer
				delete (*iter2);
				
  }
	
	//Remove QLabels
	std::vector<QLabel*>::iterator iter_label; 
	
	for(iter_label = gradientlabel_.begin(); iter_label < gradientlabel_.end(); iter_label++ ) 
	{
				viewlayout_->removeWidget((*iter_label));
				(*iter_label)->hide();
				(*iter_label) =0;	
				delete (*iter_label);
						
  }
	
	gradientdata_.clear();
	gradientlabel_.clear();
		
  //this->repaint();
}

void GradientVisualizer::update_()
{
	
	try
	{
	  //Delete all data in GUI		
		removeData_();
	
		//Update GUI
		loadData_();		
		
	}
	catch(exception e)
	{
		cout<<"Error while trying to update GUI. "<<e.what()<<endl;
	}
	
}


void GradientVisualizer::reject_()
{ 
  
	try
	{
	  
		//Delete all data in GUI		
		removeData_();		
		//Restore original data and refill GUI
		tempgradient_=(*ptr_);
		loadData_();		
				
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original gradient data. "<<e.what()<<endl;
	} 
	
}

}
