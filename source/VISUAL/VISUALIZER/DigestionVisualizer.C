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
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/DigestionVisualizer.h>

//QT
#include <QtGui/QValidator>
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  DigestionVisualizer::DigestionVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Digestion>()
  {
    addLabel_("Modify Digestion information");
    addSeparator_();
    addLineEdit_(treatmenttype_, "Treatment type");
    addTextEdit_(treatmentcomment_, "Comment");
    addLineEdit_(digestionenzyme_, "Enzyme");
    addDoubleLineEdit_(digestiontime_, "Digestion time (in minutes)");
    addDoubleLineEdit_(digestiontemperature_, "Temperature (in \xa7 C)");
    addDoubleLineEdit_(digestionPH_, "PH");

    finishAdding_();
  }

  void DigestionVisualizer::update_()
  {
    treatmenttype_->setText(temp_.getType().c_str());
    treatmenttype_->setReadOnly(true);
    treatmentcomment_->setText(temp_.getComment().c_str());
    digestionenzyme_->setText(temp_.getEnzyme().c_str());
    digestiontime_->setText(String(temp_.getDigestionTime()).c_str());
    digestiontemperature_->setText(String(temp_.getTemperature()).c_str());
    digestionPH_->setText(String(temp_.getPh()).c_str());
  }

  void DigestionVisualizer::store()
  {
    ptr_->setComment(treatmentcomment_->toPlainText());
    ptr_->setEnzyme(digestionenzyme_->text());
    ptr_->setDigestionTime(digestiontime_->text().toFloat());
    ptr_->setTemperature(digestiontime_->text().toFloat());
    ptr_->setPh(digestiontime_->text().toFloat());

    temp_ = (*ptr_);
  }

  void DigestionVisualizer::undo_()
  {
    update_();
  }

}
