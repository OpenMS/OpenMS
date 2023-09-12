// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/MISC/FilterableList.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <ui_FilterableList.h>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    FilterableList::FilterableList(QWidget *parent) :
      QWidget(parent),
      ui_(new Ui::FilterableList)
    {
      ui_->setupUi(this);
      connect(ui_->filter_text, &QLineEdit::textChanged, this, &FilterableList::filterEdited_);
      // forward double-clicked signal to outside
      connect(ui_->list_items, &QListWidget::itemDoubleClicked, [&](QListWidgetItem* item) {
        emit itemDoubleClicked(item);
      });
    }

    FilterableList::~FilterableList()
    {
      delete ui_;
    }

    void FilterableList::setItems(const QStringList& items)
    {
      items_ = items;
      updateInternalList_();
    }

    void FilterableList::setBlacklistItems(const QStringList& bl_items)
    {
      /*
       * Suppressing warning toSet() deprecated till Qt 5.14
       */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      blacklist_ = bl_items.toSet();
#pragma GCC diagnostic pop
      updateInternalList_();
    }

    void FilterableList::addBlackListItems(const QStringList& items)
    {
      /*
       * Suppressing warning toSet() deprecated till Qt 5.14
       */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      blacklist_.unite(items.toSet());
#pragma GCC diagnostic pop
      updateInternalList_();
    }

    void FilterableList::removeBlackListItems(const QStringList& outdated_blacklist_items)
    {
      // quadratic runtime, but maintains order of items (as opposed to converting to set)
      /*
       * Suppressing warning toSet() deprecated till Qt 5.14
       */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      for (const auto& bl : outdated_blacklist_items.toSet())
#pragma GCC diagnostic pop
      {
        if (blacklist_.remove(bl) == 0)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value cannot be taken from blacklist. Does not belong to set!", bl.toStdString());
        }
      }
      updateInternalList_();
    }

    QStringList FilterableList::getSelectedItems() const
    {
      QStringList items;
      for (const auto& item : ui_->list_items->selectedItems()) items << item->text();
      return items;
    }

    QStringList FilterableList::getAllVisibleItems() const
    {
      QStringList items;
      for (int row = 0; row < ui_->list_items->count(); ++row) items << ui_->list_items->item(row)->text();
      return items;
    }

    void FilterableList::filterEdited_(const QString& filter_text)
    {
      // update list of visible items
      updateVisibleList_();
      // let outside world know about it
      emit filterChanged(filter_text);
    }

    void FilterableList::updateInternalList_()
    {
      items_wo_bl_ = items_;
      // quadratic runtime, but maintains order of items (as opposed to converting to set)
      for (const auto& bl : blacklist_)
      {
        if (items_wo_bl_.removeAll(bl) == 0)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value does not belong to set!", bl.toStdString());
        }
      }
      updateVisibleList_();
    }

    void FilterableList::updateVisibleList_()
    {
      QRegExp regex(ui_->filter_text->text(), Qt::CaseInsensitive, QRegExp::WildcardUnix);
      ui_->list_items->clear();
      ui_->list_items->addItems(items_wo_bl_.filter(regex));
    }

  } //namespace Internal
} //namspace OpenMS

