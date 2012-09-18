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

#include <OpenMS/VISUAL/VISUALIZER/IonDetectorVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>


//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  IonDetectorVisualizer::IonDetectorVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<IonDetector>()
  {
    addLabel_("Modify iondetector information.");
    addSeparator_();

    addIntLineEdit_(order_, "Order");
    addComboBox_(type_, "Type");
    addComboBox_(ac_mode_, "Acquisition mode");
    addDoubleLineEdit_(res_, "Resolution (in ns)");
    addDoubleLineEdit_(freq_, "ADC sampling frequency (in Hz)");

    finishAdding_();
  }

  void IonDetectorVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(type_, &temp_.NamesOfType[temp_.getType()], 1);
      fillComboBox_(ac_mode_, &temp_.NamesOfAcquisitionMode[temp_.getAcquisitionMode()], 1);
    }
    else
    {
      fillComboBox_(type_, temp_.NamesOfType, IonDetector::SIZE_OF_TYPE);
      fillComboBox_(ac_mode_, temp_.NamesOfAcquisitionMode, IonDetector::SIZE_OF_ACQUISITIONMODE);
      type_->setCurrentIndex(temp_.getType());
      ac_mode_->setCurrentIndex(temp_.getAcquisitionMode());
    }

    order_->setText(String(temp_.getOrder()).c_str());
    res_->setText(String(temp_.getResolution()).c_str());
    freq_->setText(String(temp_.getADCSamplingFrequency()).c_str());
  }

  void IonDetectorVisualizer::store()
  {
    ptr_->setOrder(order_->text().toInt());
    ptr_->setResolution(res_->text().toDouble());
    ptr_->setADCSamplingFrequency(freq_->text().toDouble());
    ptr_->setType((IonDetector::Type)type_->currentIndex());
    ptr_->setAcquisitionMode((IonDetector::AcquisitionMode)ac_mode_->currentIndex());

    temp_ = (*ptr_);
  }

  void IonDetectorVisualizer::undo_()
  {
    update_();
  }

}
