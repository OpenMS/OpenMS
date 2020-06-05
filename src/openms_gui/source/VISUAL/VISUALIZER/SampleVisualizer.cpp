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
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s


#include <OpenMS/VISUAL/VISUALIZER/SampleVisualizer.h>

//QT
#include <QtWidgets/QTextEdit>
#include <QValidator>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

  SampleVisualizer::SampleVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Sample>()
  {
    addLabel_("Modify Sample information");
    addSeparator_();
    addLineEdit_(samplename_, "Name");
    addLineEdit_(samplenumber_, "Number");
    addLineEdit_(sampleorganism_, "Organism");
    addTextEdit_(samplecomment_, "Comment");
    addComboBox_(samplestate_, "State");
    addDoubleLineEdit_(samplemass_, "Mass (in gram)");
    addDoubleLineEdit_(samplevolume_, "Volume (in ml)");
    addDoubleLineEdit_(sampleconcentration_, "Concentration (in g/l)");

    finishAdding_();
  }

  void SampleVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(samplestate_, &temp_.NamesOfSampleState[temp_.getState()], 1);
    }
    else
    {
      fillComboBox_(samplestate_, temp_.NamesOfSampleState, Sample::SIZE_OF_SAMPLESTATE);
      samplestate_->setCurrentIndex(temp_.getState());
    }

    samplename_->setText(temp_.getName().c_str());
    samplenumber_->setText(temp_.getNumber().c_str());
    sampleorganism_->setText(temp_.getOrganism().c_str());
    samplecomment_->setText(temp_.getComment().c_str());

    samplemass_->setText(String(temp_.getMass()).c_str());
    samplevolume_->setText(String(temp_.getVolume()).c_str());
    sampleconcentration_->setText(String(temp_.getConcentration()).c_str());
  }

  void SampleVisualizer::store()
  {
    ptr_->setName(samplename_->text());
    ptr_->setNumber(samplenumber_->text());
    ptr_->setOrganism(sampleorganism_->text());
    ptr_->setComment(samplecomment_->toPlainText());
    ptr_->setState((Sample::SampleState)samplestate_->currentIndex());
    ptr_->setMass(samplemass_->text().toFloat());
    ptr_->setVolume(samplevolume_->text().toFloat());
    ptr_->setConcentration(sampleconcentration_->text().toFloat());

    temp_ = (*ptr_);
  }

  void SampleVisualizer::undo_()
  {
    update_();
  }

}
