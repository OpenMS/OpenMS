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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/ListFilterDialog.h>

#include <OpenMS/VISUAL/MISC/FilterableList.h>

#include <ui_ListFilterDialog.h>

#include <QCloseEvent>

using namespace std;

namespace OpenMS
{

  /// C'tor with items to show and select from
  ListFilterDialog::ListFilterDialog(QWidget* parent, const QStringList& items, const QStringList& items_prechosen)
      : QDialog(parent),
      ui_(new Ui::ListFilterDialog)
  {
    ui_->setupUi(this);

    connect(ui_->ok_button_, &QPushButton::clicked, this, &QDialog::accept);
    connect(ui_->cancel_button_, &QPushButton::clicked, this, &QDialog::reject);

    connect(ui_->btn_left_right, &QPushButton::clicked, this, &ListFilterDialog::BtnLRClicked_);
    connect(ui_->btn_left_right_all, &QPushButton::clicked, this, &ListFilterDialog::BtnLRAllClicked_);
    connect(ui_->btn_right_left, &QPushButton::clicked, this, &ListFilterDialog::BtnRLClicked_);
    connect(ui_->btn_right_left_all, &QPushButton::clicked, this, &ListFilterDialog::BtnRLAllClicked_);
    // if something was double clicked, it's also selected. So just forward to '>>' button
    connect(ui_->list_in, &FilterableList::itemDoubleClicked, this, &ListFilterDialog::BtnLRClicked_);
    connect(ui_->list_out, &QListWidget::itemDoubleClicked, this, &ListFilterDialog::BtnRLClicked_);

    setItems(items);
    setPrechosenItems(items_prechosen);
  }     
  

  ListFilterDialog::~ListFilterDialog()
  {
    delete ui_;
  }

  void ListFilterDialog::closeEvent(QCloseEvent* event)
  {
    event->accept();
    reject();
  }

  void ListFilterDialog::setItems(const QStringList& items)
  {
    ui_->list_in->setItems(items);
  }

  void ListFilterDialog::setPrechosenItems(const QStringList& items_pre)
  {
    ui_->list_in->setBlacklistItems(items_pre);
    ui_->list_out->clear();
    ui_->list_out->addItems(items_pre);
  }

  QStringList ListFilterDialog::getChosenItems() const
  {
    QStringList items;
    for (int row = 0; row < ui_->list_out->count(); ++row) items << ui_->list_out->item(row)->text();
    return items;
  }
  
  void ListFilterDialog::BtnLRClicked_()
  {
    auto selected = ui_->list_in->getSelectedItems();
    ui_->list_out->insertItems(ui_->list_out->count(), selected);
    ui_->list_in->addBlackListItems(selected);
  }

  void ListFilterDialog::BtnLRAllClicked_()
  {
    auto selected = ui_->list_in->getAllVisibleItems();
    ui_->list_out->insertItems(ui_->list_out->count(), selected);
    ui_->list_in->addBlackListItems(selected);
  }

  void ListFilterDialog::BtnRLClicked_()
  {
    QStringList selected;
    auto widget_selected = ui_->list_out->selectedItems();
    for (auto& item : widget_selected) selected << item->text();
    qDeleteAll(widget_selected);
    ui_->list_in->removeBlackListItems(selected);
  }

  void ListFilterDialog::BtnRLAllClicked_()
  {
    QStringList selected = getChosenItems();
    ui_->list_out->clear();
    ui_->list_in->removeBlackListItems(selected);
  }



}
