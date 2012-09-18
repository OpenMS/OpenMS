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
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/SpectrumSettingsVisualizer.h>

//QT
#include <QtGui/QComboBox>
#include <QtGui/QTextEdit>
#include <QtGui/QLineEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  SpectrumSettingsVisualizer::SpectrumSettingsVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<SpectrumSettings>()
  {
    addLabel_("Modify the settings of the spectrum.");
    addSeparator_();
    addComboBox_(type_, "Type of spectrum");
    addLineEdit_(native_id_, "Native ID");
    addTextEdit_(comment_, "Comment");

    finishAdding_();
  }

  void SpectrumSettingsVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(type_, &temp_.NamesOfSpectrumType[temp_.getType()], 1);
    }
    else
    {
      fillComboBox_(type_, temp_.NamesOfSpectrumType, SpectrumSettings::SIZE_OF_SPECTRUMTYPE);
      type_->setCurrentIndex(temp_.getType());
    }

    native_id_->setText(temp_.getNativeID().c_str());
    comment_->setText(temp_.getComment().c_str());
  }

  void SpectrumSettingsVisualizer::store()
  {
    ptr_->setType((SpectrumSettings::SpectrumType)type_->currentIndex());
    ptr_->setNativeID(native_id_->text());
    ptr_->setComment(comment_->toPlainText());

    temp_ = (*ptr_);
  }

  void SpectrumSettingsVisualizer::undo_()
  {
    update_();
  }

}
