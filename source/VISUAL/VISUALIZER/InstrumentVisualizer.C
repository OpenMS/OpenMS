// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/InstrumentVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>
#include <string>

using namespace std;

namespace OpenMS
{

  InstrumentVisualizer::InstrumentVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Instrument>()
  {
    addLabel_("Modify instrument information.");
    addSeparator_();
    addLineEdit_(name_, "Name");
    addLineEdit_(vendor_, "Vendor");
    addLineEdit_(model_, "Model");
    addTextEdit_(customizations_, "Customizations");
    addComboBox_(ion_optics_, "Ion optics");

    finishAdding_();
  }

  void InstrumentVisualizer::update_()
  {
    name_->setText(temp_.getName().c_str());
    vendor_->setText(temp_.getVendor().c_str());
    model_->setText(temp_.getModel().c_str());
    customizations_->setText(temp_.getCustomizations().c_str());

    if (!isEditable())
    {
      fillComboBox_(ion_optics_, &temp_.NamesOfIonOpticsType[temp_.getIonOptics()], 1);
    }
    else
    {
      fillComboBox_(ion_optics_, temp_.NamesOfIonOpticsType, Instrument::SIZE_OF_IONOPTICSTYPE);
    }
  }

  void InstrumentVisualizer::store()
  {
    ptr_->setName(name_->text());
    ptr_->setVendor(vendor_->text());
    ptr_->setModel(model_->text());
    ptr_->setCustomizations(customizations_->toPlainText());
    ptr_->setIonOptics((Instrument::IonOpticsType)ion_optics_->currentIndex());

    temp_ = (*ptr_);
  }

  void InstrumentVisualizer::undo_()
  {
    update_();
  }

}
