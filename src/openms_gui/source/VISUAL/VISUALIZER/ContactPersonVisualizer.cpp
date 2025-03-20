// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/ContactPersonVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>

using namespace std;

namespace OpenMS
{

  ContactPersonVisualizer::ContactPersonVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<ContactPerson>()
  {

    addLabel_("Modify ContactPerson information");
    addSeparator_();
    addLineEdit_(firstname_, "First name");
    addLineEdit_(lastname_, "Last name");
    addLineEdit_(institution_, "Institution");
    addLineEdit_(address_, "Address");
    addLineEdit_(email_, "Email");
    addLineEdit_(url_, "URL");
    addLineEdit_(contact_info_, "Contact info");
    finishAdding_();
  }

  void ContactPersonVisualizer::update_()
  {
    firstname_->setText(temp_.getFirstName().c_str());
    lastname_->setText(temp_.getLastName().c_str());
    institution_->setText(temp_.getInstitution().c_str());
    email_->setText(temp_.getEmail().c_str());
    contact_info_->setText(temp_.getContactInfo().c_str());
    url_->setText(temp_.getURL().c_str());
    address_->setText(temp_.getAddress().c_str());
  }

  void ContactPersonVisualizer::store()
  {
    ptr_->setLastName(lastname_->text());
    ptr_->setInstitution(institution_->text());
    ptr_->setEmail(email_->text());
    ptr_->setContactInfo(contact_info_->text());
    ptr_->setURL(url_->text());
    ptr_->setAddress(address_->text());

    temp_ = (*ptr_);
  }

  void ContactPersonVisualizer::undo_()
  {
    update_();
  }

}
