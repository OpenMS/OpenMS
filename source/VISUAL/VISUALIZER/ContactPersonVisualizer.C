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
// $Maintainer: Marc Sturm   $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/ContactPersonVisualizer.h>

//QT
#include <QtGui/QLineEdit>

using namespace std;

namespace OpenMS
{

//Constructor
ContactPersonVisualizer::ContactPersonVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
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
  cp_firstname_->setText(h.getFirstName().c_str());
  cp_lastname_->setText(h.getLastName().c_str());
	cp_institution_->setText(h.getInstitution().c_str() );
  cp_email_->setText(h.getEmail().c_str() );
	cp_contact_info_->setText(h.getContactInfo().c_str() );
		
	cp_firstname_->adjustSize();		
	cp_lastname_->adjustSize();		
	cp_lastname_->adjustSize();		
	cp_institution_->adjustSize();		
	cp_email_->adjustSize();		
	cp_email_->repaint();	
	cp_email_->show();
	cp_contact_info_->adjustSize();		
	
}

void ContactPersonVisualizer::store_()
{
	try
	{
				
		(*ptr_).setFirstName(cp_firstname_->text().toStdString());

		(*ptr_).setLastName(cp_lastname_->text().toStdString());
				
		(*ptr_).setInstitution(cp_institution_->text().toStdString());
		
		(*ptr_).setEmail(cp_email_->text().toStdString());
		
		(*ptr_).setContactInfo(cp_contact_info_->text().toStdString());
		
		tempContactPerson_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new contact person data. "<<e.what()<<endl;
	}
}

void ContactPersonVisualizer::reject_()
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

}
