// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
