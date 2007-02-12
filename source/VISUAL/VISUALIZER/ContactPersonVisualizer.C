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
// $Maintainer: stefan_heess   $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/VISUALIZER/ContactPersonVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

//QT
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qpushbutton.h>
#include <iostream>
#include <vector>


//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
ContactPersonVisualizer::ContactPersonVisualizer(bool editable, QWidget *parent, const char *name) 
	: BaseVisualizer(editable, parent, name)
{
  
	addLabel("Modify ContactPerson information");		
	addSeperator();
	addLineEdit(cp_firstname_, "First name" );
	addLineEdit(cp_lastname_, "Last name" );
	addLineEdit(cp_institution_, "Institution" );
	addLineEdit(cp_email_, "Email" );
	addLineEdit(cp_contact_info_, "Contact info" );
	
	finishAdding_();
	
}


void ContactPersonVisualizer::load(ContactPerson &h)
{
  ptr_ = &h;
	
	//Copy of current object for restoring the original values
	tempContactPerson_=h;
  cp_firstname_->setText(h.getFirstName());
  cp_lastname_->setText(h.getLastName());
	cp_institution_->setText(h.getInstitution() );
  cp_email_->setText(h.getEmail() );
	cp_contact_info_->setText(h.getContactInfo() );
		
			
}

void ContactPersonVisualizer::store()
{
	try
	{
				
		(*ptr_).setFirstName(string((const char*)cp_firstname_->text()));

		(*ptr_).setLastName(string((const char*)cp_lastname_->text()));
				
		(*ptr_).setInstitution(string((const char*)cp_institution_->text()) );
		
		(*ptr_).setEmail(string((const char*)cp_email_->text()) );
		
		(*ptr_).setContactInfo(string((const char*)cp_contact_info_->text()) );
		
		tempContactPerson_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new contact person data. "<<e.what()<<endl;
	}
}

void ContactPersonVisualizer::reject()
{
	try
	{
		load(tempContactPerson_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original contact person data. "<<e.what()<<endl;
	} 
}
