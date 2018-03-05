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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

#include <QtGui/QLayout>
#include <QtGui/QWidget>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>
#include <QtGui/QLabel>
#include <QtGui/QTextEdit>
#include <QtGui/QGridLayout>
#include <QtGui/QPushButton>
#include <QtGui/QHBoxLayout>
#include <QtGui/QListWidget>

namespace OpenMS
{
  BaseVisualizerGUI::BaseVisualizerGUI(bool editable, QWidget * parent) :
    QWidget(parent),
    undo_button_(nullptr),
    mainlayout_(nullptr),
    row_(0),
    editable_(editable)
  {
    mainlayout_ = new QGridLayout(this);
    mainlayout_->setMargin(0);
  }

  bool BaseVisualizerGUI::isEditable() const
  {
    return editable_;
  }

  void BaseVisualizerGUI::finishAdding_()
  {
    if (isEditable())
    {
      addSeparator_();
      addButton_(undo_button_, "Undo");
      connect(undo_button_, SIGNAL(clicked()), this, SLOT(undo_()));
    }
    addVSpacer_();
  }

  void BaseVisualizerGUI::addLabel_(QString label, UInt row)
  {
    QLabel * label_item = new QLabel(label, this);
    mainlayout_->addWidget(label_item, row, 0);
  }

  void BaseVisualizerGUI::addLabel_(QString label)
  {
    QLabel * label_item = new QLabel(label, this);
    mainlayout_->addWidget(label_item, row_, 0, 1, 3);
    row_++;
  }

  void BaseVisualizerGUI::addLineEdit_(QLineEdit * & ptr, QString label)
  {
    ptr = new QLineEdit(this);
    ptr->setMinimumWidth(180);
    addLabel_(label, row_);
    mainlayout_->addWidget(ptr, row_, 1, 1, 2);
    ptr->setReadOnly(!isEditable());
    row_++;
  }

  void BaseVisualizerGUI::addIntLineEdit_(QLineEdit * & ptr, QString label)
  {
    ptr = new QLineEdit(this);
    ptr->setMinimumWidth(180);
    QIntValidator * vali = new QIntValidator(ptr);
    ptr->setValidator(vali);
    addLabel_(label, row_);
    mainlayout_->addWidget(ptr, row_, 1, 1, 2);
    ptr->setReadOnly(!isEditable());
    row_++;

  }

  void BaseVisualizerGUI::addDoubleLineEdit_(QLineEdit * & ptr, QString label)
  {
    ptr = new QLineEdit(this);
    ptr->setMinimumWidth(180);
    QDoubleValidator * vali = new QDoubleValidator(ptr);
    ptr->setValidator(vali);
    addLabel_(label, row_);
    mainlayout_->addWidget(ptr, row_, 1, 1, 2);
    ptr->setReadOnly(!isEditable());
    row_++;

  }

  void BaseVisualizerGUI::addLineEditButton_(QString label, QLineEdit * & ptr1, QPushButton * & ptr2, QString buttonlabel)
  {
    QLabel * label_item = new QLabel(label, this);
    ptr1 = new QLineEdit(this);
    ptr1->setMinimumWidth(180);
    ptr2 = new QPushButton(buttonlabel, this);
    mainlayout_->addWidget(label_item, row_, 0);
    mainlayout_->addWidget(ptr1, row_, 1);
    mainlayout_->addWidget(ptr2, row_, 2);

    ptr1->setReadOnly(!isEditable());
    ptr2->setEnabled(isEditable());
    row_++;
  }

  void BaseVisualizerGUI::addTextEdit_(QTextEdit * & ptr, QString label)
  {
    ptr = new QTextEdit(this);
    addLabel_(label, row_);
    mainlayout_->addWidget(ptr, row_, 1, 1, 2);
    ptr->setReadOnly(!isEditable());
    row_++;
  }

  void BaseVisualizerGUI::addComboBox_(QComboBox * & ptr, QString label)
  {
    ptr = new QComboBox(this);
    addLabel_(label, row_);
    mainlayout_->addWidget(ptr, row_, 1, 1, 2);
    ptr->blockSignals(true);
    row_++;
  }

  void BaseVisualizerGUI::addBooleanComboBox_(QComboBox * & ptr, QString label)
  {
    ptr = new QComboBox(this);
    ptr->insertItem(0, "false");
    ptr->insertItem(1, "true");
    addLabel_(label, row_);
    mainlayout_->addWidget(ptr, row_, 1, 1, 2);
    ptr->blockSignals(true);
    row_++;
  }

  void BaseVisualizerGUI::fillComboBox_(QComboBox * & ptr, const std::string * items, int agr)
  {
    for (int i = 0; i < agr; ++i)
    {
      ptr->insertItem(i, QString(items[i].c_str()));
    }
  }

  void BaseVisualizerGUI::addButton_(QPushButton * & ptr, QString label)
  {
    ptr = new QPushButton(label, this);
    QHBoxLayout * box = new QHBoxLayout();
    box->addStretch(1);
    box->addWidget(ptr);
    mainlayout_->addLayout(box, row_, 0, 1, 3);

    ptr->setEnabled(isEditable());
    row_++;
  }

  void BaseVisualizerGUI::addVSpacer_()
  {
    mainlayout_->setRowStretch(row_, 1);
    row_++;
  }

  void BaseVisualizerGUI::add2Buttons_(QPushButton * & ptr1, QString label1, QPushButton * & ptr2, QString label2)
  {
    ptr1 = new QPushButton(label1, this);
    ptr2 = new QPushButton(label2, this);
    QHBoxLayout * box = new QHBoxLayout();
    box->addStretch(1);
    box->addWidget(ptr1);
    box->addWidget(ptr2);
    mainlayout_->addLayout(box, row_, 0, 1, 3);
    row_++;
  }

  void BaseVisualizerGUI::addSeparator_()
  {
    QLabel * pLabel = new QLabel(this);
    pLabel->setFrameShape(QFrame::HLine);
    mainlayout_->addWidget(pLabel, row_, 0, 1, 3);
    row_++;
  }

  void BaseVisualizerGUI::addListView_(QListWidget * & ptr, QString label)
  {
    ptr = new QListWidget(this);
    addLabel_(label, row_);
    mainlayout_->addWidget(ptr, row_, 1, 1, 2);
    row_++;
  }

}
