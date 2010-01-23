// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
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
		addLineEdit_(firstname_, "First name" );
		addLineEdit_(lastname_, "Last name" );
		addLineEdit_(institution_, "Institution" );
		addLineEdit_(address_, "Address" );
		addLineEdit_(email_, "Email" );
		addLineEdit_(url_, "URL" );
		addLineEdit_(contact_info_, "Contact info" );
		finishAdding_();
	}
	
	void ContactPersonVisualizer::update_()
	{
	  firstname_->setText(temp_.getFirstName().c_str());
	  lastname_->setText(temp_.getLastName().c_str());
		institution_->setText(temp_.getInstitution().c_str() );
	  email_->setText(temp_.getEmail().c_str() );
		contact_info_->setText(temp_.getContactInfo().c_str() );
	  url_->setText(temp_.getURL().c_str() );
		address_->setText(temp_.getAddress().c_str() );
	}
	
	void ContactPersonVisualizer::store()
	{
		ptr_->setLastName(lastname_->text());
		ptr_->setInstitution(institution_->text());
		ptr_->setEmail(email_->text());
		ptr_->setContactInfo(contact_info_->text());
		ptr_->setURL(url_->text());
		ptr_->setAddress(address_->text());
		
		temp_=(*ptr_);
	}
	
	void ContactPersonVisualizer::undo_()
	{
		update_();
	}

}
