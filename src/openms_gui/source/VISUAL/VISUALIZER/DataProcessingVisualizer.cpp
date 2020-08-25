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

#include <OpenMS/VISUAL/VISUALIZER/DataProcessingVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QListWidget>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  DataProcessingVisualizer::DataProcessingVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<DataProcessing>()
  {
    addLabel_("Modify data processing information.");
    addSeparator_();

    addLineEdit_(completion_time_, "Completion time");
    addListView_(actions_, "Processing actions");
    finishAdding_();
  }

  void DataProcessingVisualizer::update_()
  {
    //time
    completion_time_->setText(temp_.getCompletionTime().get().c_str());

    //actions
    actions_->clear();
    for (Size i = 0; i < DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
    {
      QListWidgetItem * item = new QListWidgetItem(actions_);
      item->setText(QString::fromStdString(DataProcessing::NamesOfProcessingAction[i]));
      if (temp_.getProcessingActions().count(DataProcessing::ProcessingAction(i)) == 1)
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
      actions_->addItem(item);
    }
  }

  void DataProcessingVisualizer::store()
  {
    DateTime date;
    try
    {
      date.set(completion_time_->text());
      ptr_->setCompletionTime(date);
    }
    catch (exception & /*e*/)
    {
      if (date.isNull())
      {
        std::string status = "Format of date in DATAPROCESSING is not correct.";
        emit sendStatus(status);
      }
    }

    //actions
    ptr_->getProcessingActions().clear();
    for (UInt i = 0; i < DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
    {
      if (actions_->item(i)->checkState() == Qt::Checked)
      {
        ptr_->getProcessingActions().insert(DataProcessing::ProcessingAction(i));
      }
    }

    temp_ = (*ptr_);
  }

  void DataProcessingVisualizer::undo_()
  {
    update_();
  }

}
