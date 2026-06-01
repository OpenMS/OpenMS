// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
