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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/PrecursorVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>
#include <QtGui/QListWidget>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  PrecursorVisualizer::PrecursorVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Precursor>()
  {
    addLabel_("Modify processing method information.");

    addSeparator_();

    addDoubleLineEdit_(mz_, "m/z");
    addDoubleLineEdit_(int_, "intensity");
    addIntLineEdit_(charge_, "charge");

    addDoubleLineEdit_(window_low_, "Lower offset from target m/z");
    addDoubleLineEdit_(window_up_, "Upper offset from target m/z");

    addListView_(activation_methods_, "Activation methods");
    addDoubleLineEdit_(activation_energy_, "Activation energy");

    finishAdding_();
  }

  void PrecursorVisualizer::update_()
  {
    mz_->setText(String(temp_.getMZ()).c_str());
    int_->setText(String(temp_.getIntensity()).c_str());
    charge_->setText(String(temp_.getCharge()).c_str());

    window_low_->setText(String(temp_.getIsolationWindowLowerOffset()).c_str());
    window_up_->setText(String(temp_.getIsolationWindowUpperOffset()).c_str());

    //actions
    activation_methods_->clear();
    for (Size i = 0; i < Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
    {
      QListWidgetItem * item = new QListWidgetItem(activation_methods_);
      item->setText(QString::fromStdString(Precursor::NamesOfActivationMethod[i]));
      if (temp_.getActivationMethods().count(Precursor::ActivationMethod(i)) == 1)
      {
        item->setCheckState(Qt::Checked);
      }
      else
      {
        item->setCheckState(Qt::Unchecked);
      }
      if (isEditable())
      {
        item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
      }
      else
      {
        item->setFlags(Qt::ItemIsEnabled);
      }
      activation_methods_->addItem(item);
    }

    activation_energy_->setText(String(temp_.getActivationEnergy()).c_str());
  }

  void PrecursorVisualizer::store()
  {
    ptr_->setMZ(mz_->text().toFloat());
    ptr_->setIntensity(int_->text().toFloat());
    ptr_->setCharge(charge_->text().toInt());

    ptr_->setIsolationWindowLowerOffset(window_low_->text().toFloat());
    ptr_->setIsolationWindowUpperOffset(window_up_->text().toFloat());

    ptr_->getActivationMethods().clear();
    for (UInt i = 0; i < Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
    {
      if (activation_methods_->item(i)->checkState() == Qt::Checked)
      {
        ptr_->getActivationMethods().insert(Precursor::ActivationMethod(i));
      }
    }
    ptr_->setActivationEnergy(activation_energy_->text().toFloat());

    temp_ = (*ptr_);
  }

  void PrecursorVisualizer::undo_()
  {
    update_();
  }

}
