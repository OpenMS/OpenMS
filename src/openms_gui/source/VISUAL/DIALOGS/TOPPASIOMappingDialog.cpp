// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <ui_TOPPASIOMappingDialog.h>

#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASSplitterVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>

#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <QtWidgets/QMessageBox>

#include <iostream>
#include <sstream>
#include <OpenMS/VISUAL/TOPPASOutputFolderVertex.h>

namespace OpenMS
{
  TOPPASIOMappingDialog::TOPPASIOMappingDialog(TOPPASEdge* parent)
    : ui_(new Ui::TOPPASIOMappingDialogTemplate)
  {
    ui_->setupUi(this);
    edge_ = parent;
    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(checkValidity_()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));

    fillComboBoxes_();
  }

  TOPPASIOMappingDialog::~TOPPASIOMappingDialog()
  {
    delete ui_;
  }

  int TOPPASIOMappingDialog::firstExec()
  {
    // check if only 1 parameter, if yes: select it
    if (ui_->source_combo->count() == 2)     // <select> + 1 parameter
    {
      ui_->source_combo->setCurrentIndex(1);
    }
    if (ui_->target_combo->count() == 2)
    {
      ui_->target_combo->setCurrentIndex(1);
    }
    // if the target is an output folder and there is no output parameter in the input tool, the egde is invalid
    TOPPASOutputFolderVertex* target_dir = qobject_cast<TOPPASOutputFolderVertex*>(edge_->getTargetVertex());
    if (target_dir && ui_->source_combo->count() == 0)
    {
      return QDialog::Rejected;
    }

    // is there only 1 possible mapping? -> do not show dialog
    if ((ui_->source_combo->count() == 2 || ui_->source_combo->count() == 0) &&
        (ui_->target_combo->count() == 2 || ui_->target_combo->count() == 0))
    {
      checkValidity_();
      return QDialog::Accepted;
    }
    else
    {
      return QDialog::exec();
    }
  }

  void TOPPASIOMappingDialog::fillComboBoxes_()
  {
    TOPPASVertex* source = edge_->getSourceVertex();
    TOPPASVertex* target = edge_->getTargetVertex();

    TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
    TOPPASToolVertex* target_tool = qobject_cast<TOPPASToolVertex*>(target);
    TOPPASMergerVertex* source_merger = qobject_cast<TOPPASMergerVertex*>(source);
    TOPPASMergerVertex* target_merger = qobject_cast<TOPPASMergerVertex*>(target);
    TOPPASSplitterVertex* source_splitter = qobject_cast<TOPPASSplitterVertex*>(source);
    TOPPASSplitterVertex* target_splitter = qobject_cast<TOPPASSplitterVertex*>(target);
    TOPPASInputFileListVertex* source_list = qobject_cast<TOPPASInputFileListVertex*>(source);
    TOPPASOutputFileListVertex* target_list = qobject_cast<TOPPASOutputFileListVertex*>(target);
    TOPPASOutputFolderVertex* target_dir = qobject_cast<TOPPASOutputFolderVertex*>(target);

    // an output folder can only be connected to a tool
    if (target_dir)
    {
      if (!source_tool)
      { // bad news: no source tool, hence no connection possible
        return;
      }
      const auto source_output_dirs = source_tool->getOutputParameters();
      ui_->source_combo->addItem("<select>");
      for (const auto& info : source_output_dirs)
      {
        if (info.type != TOPPASToolVertex::IOInfo::IOT_DIR) continue;

        String item_name = "Directory: " + info.param_name + " ";
        ui_->source_combo->addItem(item_name.toQString(), source_output_dirs.indexOf(info));
      }
      if (ui_->source_combo->count() == 1) // no directories found; return empty to signal invalid edge
      {
        ui_->source_combo->clear();
        return;
      }
    }
    else if (source_tool)
    {
      QVector<TOPPASToolVertex::IOInfo> source_output_files = source_tool->getOutputParameters();
      ui_->source_label->setText(source_tool->getName().toQString());
      if (!source_tool->getType().empty())
      {
        ui_->source_type_label->setText("(" + source_tool->getType().toQString() + ")");
      }
      else
      {
        ui_->source_type_label->setVisible(false);
      }
      ui_->source_combo->addItem("<select>");
      for (TOPPASToolVertex::IOInfo info : source_output_files)
      {
        if (info.type == TOPPASToolVertex::IOInfo::IOT_DIR) continue;
        String item_name;
        if (info.type == TOPPASToolVertex::IOInfo::IOT_FILE)
        {
          if (target_splitter) continue; // inputs for splitters must be lists
          item_name = "File: ";
        }
        else
        {
          item_name = "List: ";
        }
        item_name += info.param_name + " ";
        std::ostringstream ss;
        ss << info.valid_types;
        item_name += ss.str();

        ui_->source_combo->addItem(item_name.toQString(), source_output_files.indexOf(info));
      }
    }
    else if (source_list || source_merger || source_splitter)
    {
      if (source_list)
      {
        ui_->source_label->setText("List");
      }
      else if (source_merger)
      {
        ui_->source_label->setText(source_merger->roundBasedMode() ? "Merger" : "Collector");
      }
      else if (source_splitter)
      {
        ui_->source_label->setText("Splitter");
      }
      ui_->source_type_label->setVisible(false);
      ui_->source_combo->setVisible(false);
      ui_->source_parameter_label->setVisible(false);
    }

    if (target_tool)
    {
      QVector<TOPPASToolVertex::IOInfo> target_input_files = target_tool->getInputParameters();
      ui_->target_label->setText(target_tool->getName().toQString());
      if (!target_tool->getType().empty())
      {
        ui_->target_type_label->setText("(" + target_tool->getType().toQString() + ")");
      }
      else
      {
        ui_->target_type_label->setVisible(false);
      }
      ui_->target_combo->addItem("<select>");
      for (TOPPASToolVertex::IOInfo info : target_input_files)
      {
        // check if parameter occupied by another edge already
        bool occupied = false;
        for (TOPPASVertex::ConstEdgeIterator it = target->inEdgesBegin(); it != target->inEdgesEnd(); ++it)
        {
          int param_index = (*it)->getTargetInParam();
          if (*it != edge_ && param_index >= 0 && param_index < target_input_files.size())
          {
            if (info.param_name == target_input_files[param_index].param_name)
            {
              occupied = true;
              break;
            }
          }
        }
        if (occupied)
        {
          continue;
        }

        String item_name;
        if (info.type == TOPPASToolVertex::IOInfo::IOT_FILE)
        {
          if (source_merger && !source_merger->roundBasedMode()) continue; // collectors produce lists
          item_name = "File: ";
        }
        else
        {
          item_name = "List: ";
        }
        item_name += info.param_name + " ";
        std::ostringstream ss;
        ss << info.valid_types;
        item_name += ss.str();

        ui_->target_combo->addItem(item_name.toQString(), target_input_files.indexOf(info));
      }
    }
    else if (target_list || target_dir || target_merger || target_splitter)
    {
      if (target_list)
      {
        ui_->target_label->setText("List");
      }
      else if (target_dir)
      {
        ui_->target_label->setText("Directory");
      }
      else if (target_merger)
      {
        ui_->target_label->setText(target_merger->roundBasedMode() ? "Merger" : "Collector");
      }
      else if (target_splitter)
      {
        ui_->target_label->setText("Splitter");
      }
      ui_->target_type_label->setVisible(false);
      ui_->target_combo->setVisible(false);
      ui_->target_parameter_label->setVisible(false);
    }

    // pre-select the current mapping for existing edges
    // note: for new edges, @p edge_index is -1
    auto find_index = [](QComboBox* combo, int edge_index) -> int {
      for (int i = 1; i < combo->count(); ++i)
      {
        if (combo->itemData(i).toInt() == edge_index) return i;
      }
      // no mapping found
      if (combo->count() == 2) // only 1 parameter (+ <select>)
      {
        return 1; // pre-select the only parameter
      }
      return 0; // use '<select>'
    };
    ui_->source_combo->setCurrentIndex(find_index(ui_->source_combo, edge_->getSourceOutParam()));
    ui_->target_combo->setCurrentIndex(find_index(ui_->target_combo, edge_->getTargetInParam()));
    
    resize(width(), 0);
  }

  void TOPPASIOMappingDialog::checkValidity_()
  {
    const QString& source_text = ui_->source_combo->currentText();
    const QString& target_text = ui_->target_combo->currentText();

    TOPPASVertex* source = edge_->getSourceVertex();
    TOPPASVertex* target = edge_->getTargetVertex();
    TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
    TOPPASToolVertex* target_tool = qobject_cast<TOPPASToolVertex*>(target);

    if (source_text == "<select>")
    {
      QMessageBox::warning(nullptr, "Invalid selection", "You must specify the source output parameter!");
      return;
    }
    if (target_text == "<select>")
    {
      QMessageBox::warning(nullptr, "Invalid selection", "You must specify the target input parameter!");
      return;
    }

    if (source_tool)
    {
      edge_->setSourceOutParam(ui_->source_combo->currentData().toInt());
    }
    if (target_tool)
    {
      edge_->setTargetInParam(ui_->target_combo->currentData().toInt());
    }
    edge_->updateColor();

    TOPPASEdge::EdgeStatus es = edge_->getEdgeStatus();
    if (es == TOPPASEdge::ES_VALID || es == TOPPASEdge::ES_NOT_READY_YET)
    {
      accept();
    }
    else
    {
      if (es == TOPPASEdge::ES_NO_TARGET_PARAM)
      {
        QMessageBox::warning(nullptr, "Invalid selection", "You must specify the target input parameter!");
      }
      else if (es == TOPPASEdge::ES_NO_SOURCE_PARAM)
      {
        QMessageBox::warning(nullptr, "Invalid selection", "You must specify the source output parameter!");
      }
      else if (es == TOPPASEdge::ES_FILE_EXT_MISMATCH)
      {
        QMessageBox::warning(nullptr, "Invalid selection", "The file types of source output and target input parameter do not match!");
      }
      else if (es == TOPPASEdge::ES_MERGER_EXT_MISMATCH)
      {
        QMessageBox::warning(nullptr, "Invalid selection", "The file types of source output and the target input parameter do not match!");
      }
      else if (es == TOPPASEdge::ES_MERGER_WITHOUT_TOOL)
      {
        // this should be prevented already by "TOPPASScene::isEdgeAllowed_":
        QMessageBox::warning(nullptr, "Invalid selection", "Mergers or splitters connecting input and output files directly are not allowed!");
      }
      else
      {
        QMessageBox::warning(nullptr, "Ooops", "This should not have happened. Please contact the OpenMS mailing list and report this bug.");
      }
    }
  }

} // namespace
