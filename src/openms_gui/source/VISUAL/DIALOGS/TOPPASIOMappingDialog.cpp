// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASSplitterVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <QtGui/QMessageBox>

#include <iostream>
#include <sstream>

namespace OpenMS
{
  TOPPASIOMappingDialog::TOPPASIOMappingDialog(TOPPASEdge* parent)
  {
    edge_ = parent;
    setupUi(this);
    connect(ok_button, SIGNAL(clicked()), this, SLOT(checkValidity_()));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT(reject()));

    fillComboBoxes_();
  }

  int TOPPASIOMappingDialog::firstExec()
  {
    // check if only 1 parameter, if yes: select it
    if (source_combo->count() == 2)     // <select> + 1 parameter
    {
      source_combo->setCurrentIndex(1);
    }
    if (target_combo->count() == 2)
    {
      target_combo->setCurrentIndex(1);
    }

    // is there only 1 possible mapping? -> do not show dialog
    if ((source_combo->count() == 2 || source_combo->count() == 0) &&
        (target_combo->count() == 2 || target_combo->count() == 0))
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
    target_input_param_indices_.clear();

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


    if (source_tool)
    {
      QVector<TOPPASToolVertex::IOInfo> source_output_files;
      source_tool->getOutputParameters(source_output_files);
      source_label->setText(source_tool->getName().toQString());
      if (source_tool->getType() != "")
      {
        source_type_label->setText("(" + source_tool->getType().toQString() + ")");
      }
      else
      {
        source_type_label->setVisible(false);
      }
      source_combo->addItem("<select>");
      foreach(TOPPASToolVertex::IOInfo info, source_output_files)
      {
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

        source_combo->addItem(item_name.toQString());
      }
      if (source_combo->count() == 2) // only 1 parameter
      {
        source_combo->setCurrentIndex(1);
      }
    }
    else if (source_list || source_merger || source_splitter)
    {
      if (source_list)
      {
        source_label->setText("List");
      }
      else if (source_merger)
      {
        source_label->setText(source_merger->roundBasedMode() ? "Merger" : "Collector");
      }
      else if (source_splitter)
      {
        source_label->setText("Splitter");
      }
      source_type_label->setVisible(false);
      source_combo->setVisible(false);
      source_parameter_label->setVisible(false);
    }

    if (target_tool)
    {
      QVector<TOPPASToolVertex::IOInfo> target_input_files;
      target_tool->getInputParameters(target_input_files);
      target_label->setText(target_tool->getName().toQString());
      if (target_tool->getType() != "")
      {
        target_type_label->setText("(" + target_tool->getType().toQString() + ")");
      }
      else
      {
        target_type_label->setVisible(false);
      }
      target_combo->addItem("<select>");
      int param_counter = -1;
      foreach(TOPPASToolVertex::IOInfo info, target_input_files)
      {
        ++param_counter;
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

        target_combo->addItem(item_name.toQString());
        target_input_param_indices_.push_back(param_counter);
      }
      if (target_combo->count() == 2) // only 1 parameter
      {
        target_combo->setCurrentIndex(1);
      }
    }
    else if (target_list || target_merger || target_splitter)
    {
      if (target_list)
      {
        target_label->setText("List");
      }
      else if (target_merger)
      {
        target_label->setText(target_merger->roundBasedMode() ? "Merger" : "Collector");
      }
      else if (target_splitter)
      {
        target_label->setText("Splitter");
      }
      target_type_label->setVisible(false);
      target_combo->setVisible(false);
      target_parameter_label->setVisible(false);
    }

    int source_out = edge_->getSourceOutParam();
    int target_in = edge_->getTargetInParam();
    int combo_index = target_input_param_indices_.indexOf(target_in) + 1;
    if (source_out != -1)
    {
      source_combo->setCurrentIndex(source_out + 1);
    }
    if (combo_index != 0)
    {
      target_combo->setCurrentIndex(combo_index);
    }

    resize(width(), 0);
  }

  void TOPPASIOMappingDialog::checkValidity_()
  {
    const QString& source_text = source_combo->currentText();
    const QString& target_text = target_combo->currentText();

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
      edge_->setSourceOutParam(source_combo->currentIndex() - 1);
    }
    if (target_tool)
    {
      int target_index;
      int tci = target_combo->currentIndex() - 1;
      if (0 <= tci && tci < target_input_param_indices_.size())
      {
        target_index = target_input_param_indices_[tci];
      }
      else
      {
        std::cerr << "Parameter index out of bounds!" << std::endl;
        return;
      }
      edge_->setTargetInParam(target_index);
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
