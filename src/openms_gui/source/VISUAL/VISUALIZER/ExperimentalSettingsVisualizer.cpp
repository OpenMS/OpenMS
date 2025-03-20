// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/ExperimentalSettingsVisualizer.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  ExperimentalSettingsVisualizer::ExperimentalSettingsVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<ExperimentalSettings>()
  {
    addLabel_("Modify the settings of the experiment.");
    addSeparator_();
    addLineEdit_(datetime_, "Date and time of experiment");
    addTextEdit_(comment_, "Comment");
    addLineEdit_(fraction_identifier_, "Fraction identifier");

    finishAdding_();
  }

  void ExperimentalSettingsVisualizer::update_()
  {
    datetime_->setText(temp_.getDateTime().get().c_str());
    comment_->setText(temp_.getComment().c_str());
    fraction_identifier_->setText(temp_.getFractionIdentifier().c_str());
  }

  void ExperimentalSettingsVisualizer::store()
  {
    DateTime date;
    try
    {
      date.set(datetime_->text());
      ptr_->setDateTime(date);
    }
    catch (exception & /*e*/)
    {
      if (date.isNull())
      {
        std::string status = "Format of date in EXPERIMENTALSETTINGS is not correct.";
        emit sendStatus(status);
      }
    }

    ptr_->setComment(comment_->toPlainText());
    ptr_->setFractionIdentifier(fraction_identifier_->text());

    temp_ = (*ptr_);
  }

  void ExperimentalSettingsVisualizer::undo_()
  {
    update_();
  }

}
