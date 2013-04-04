// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/ContactPersonVisualizer.h>

//QT
#include <QtGui/QLineEdit>

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
