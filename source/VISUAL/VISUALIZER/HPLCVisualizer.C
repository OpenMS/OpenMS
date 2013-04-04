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
#include <OpenMS/VISUAL/VISUALIZER/HPLCVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QValidator>
#include <iostream>

using namespace std;

namespace OpenMS
{

  HPLCVisualizer::HPLCVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<HPLC>()
  {
    addLabel_("Modify HPLC information");
    addSeparator_();
    addLineEdit_(hplcinstrument_, "Instrument");
    addLineEdit_(hplccolumn_, "Column");
    addIntLineEdit_(hplctemperature_, "Temperature (in deg. C)");
    addIntLineEdit_(hplcpressure_, "Pressure (in bar)");
    addIntLineEdit_(hplcflux_, "Flux (in ul/sec)");
    addTextEdit_(hplccomment_, "Comment");

    finishAdding_();
  }

  void HPLCVisualizer::update_()
  {
    hplcinstrument_->setText(temp_.getInstrument().c_str());
    hplccolumn_->setText(temp_.getColumn().c_str());
    hplctemperature_->setText(String(temp_.getTemperature()).c_str());
    hplcpressure_->setText(String(temp_.getPressure()).c_str());
    hplcflux_->setText(String(temp_.getFlux()).c_str());
    hplccomment_->setText(temp_.getComment().c_str());
  }

  void HPLCVisualizer::store()
  {
    ptr_->setInstrument(hplcinstrument_->text());
    ptr_->setColumn(hplccolumn_->text());
    ptr_->setTemperature(hplctemperature_->text().toInt());
    ptr_->setPressure(hplcpressure_->text().toInt());
    ptr_->setFlux(hplcflux_->text().toInt());
    ptr_->setComment(hplccomment_->toPlainText());
    temp_ = (*ptr_);
  }

  void HPLCVisualizer::undo_()
  {
    update_();
  }

}
