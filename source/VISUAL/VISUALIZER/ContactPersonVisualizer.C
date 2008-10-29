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
	
	ContactPersonVisualizer::ContactPersonVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<ContactPerson>()
	{
	  
		addLabel_("Modify ContactPerson information");		
		addSeparator_();
		addLineEdit_(cp_firstname_, "First name" );
		addLineEdit_(cp_lastname_, "Last name" );
		addLineEdit_(cp_institution_, "Institution" );
		addLineEdit_(cp_email_, "Email" );
		addLineEdit_(cp_contact_info_, "Contact info" );
		finishAdding_();
	}
	
	void ContactPersonVisualizer::update_()
	{
	  cp_firstname_->setText(temp_.getFirstName().c_str());
	  cp_lastname_->setText(temp_.getLastName().c_str());
		cp_institution_->setText(temp_.getInstitution().c_str() );
	  cp_email_->setText(temp_.getEmail().c_str() );
		cp_contact_info_->setText(temp_.getContactInfo().c_str() );
			
		cp_firstname_->adjustSize();		
		cp_lastname_->adjustSize();		
		cp_lastname_->adjustSize();		
		cp_institution_->adjustSize();		
		cp_email_->adjustSize();		
		cp_email_->repaint();	
		cp_email_->show();
		cp_contact_info_->adjustSize();		
	}
	
	void ContactPersonVisualizer::store()
	{
		ptr_->setLastName(cp_lastname_->text().toStdString());
		ptr_->setInstitution(cp_institution_->text().toStdString());
		ptr_->setEmail(cp_email_->text().toStdString());
		ptr_->setContactInfo(cp_contact_info_->text().toStdString());
		
		temp_=(*ptr_);
	}
	
	void ContactPersonVisualizer::undo_()
	{
		update_();
	}

}
