// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/VISUAL/VISUALIZER/IonSourceVisualizer.h>

//QT
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  IonSourceVisualizer::IonSourceVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<IonSource>()
  {
    addLabel_("Modify ionsource information.");
    addSeparator_();

    addIntLineEdit_(order_, "Order");
    addComboBox_(inlet_type_, "Inlet type");
    addComboBox_(ionization_method_, "Ionization method");
    addComboBox_(polarity_, "Polarity");

    finishAdding_();
  }

  void IonSourceVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(inlet_type_, &temp_.NamesOfInletType[temp_.getInletType()], 1);
      fillComboBox_(ionization_method_, &temp_.NamesOfIonizationMethod[temp_.getIonizationMethod()], 1);
      fillComboBox_(polarity_, &temp_.NamesOfPolarity[temp_.getPolarity()], 1);
    }
    else
    {
      fillComboBox_(inlet_type_, temp_.NamesOfInletType, IonSource::SIZE_OF_INLETTYPE);
      fillComboBox_(ionization_method_, temp_.NamesOfIonizationMethod, IonSource::SIZE_OF_IONIZATIONMETHOD);
      fillComboBox_(polarity_, temp_.NamesOfPolarity, IonSource::SIZE_OF_POLARITY);

      inlet_type_->setCurrentIndex(temp_.getInletType());
      ionization_method_->setCurrentIndex(temp_.getIonizationMethod());
      polarity_->setCurrentIndex(temp_.getPolarity());
    }

    order_->setText(String(temp_.getOrder()).c_str());
  }

  void IonSourceVisualizer::store()
  {
    ptr_->setOrder(order_->text().toInt());
    ptr_->setInletType((IonSource::InletType)inlet_type_->currentIndex());
    ptr_->setIonizationMethod((IonSource::IonizationMethod)ionization_method_->currentIndex());
    ptr_->setPolarity((IonSource::Polarity)polarity_->currentIndex());

    temp_ = (*ptr_);
  }

  void IonSourceVisualizer::undo_()
  {
    update_();
  }

}
