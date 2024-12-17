// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/MISC/FilterableList.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Qt5Port.h>
#include <OpenMS/DATASTRUCTURES/String.h>

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
      blacklist_ = toQSet(bl_items);
      updateInternalList_();
    }

    void FilterableList::addBlackListItems(const QStringList& items)
    {
      blacklist_.unite(toQSet(items));
      updateInternalList_();
    }

    void FilterableList::removeBlackListItems(const QStringList& outdated_blacklist_items)
    {
      for (const auto& bl : outdated_blacklist_items)
      {
        if (!blacklist_.contains(bl))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value '" + String(bl) + "' cannot be taken from blacklist. Does not belong to set!", bl.toStdString());
        }
      }
      // remove all items from blacklist
      blacklist_.subtract(toQSet(outdated_blacklist_items));
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
      QRegularExpression regex(QRegularExpression::wildcardToRegularExpression(ui_->filter_text->text()),
                               QRegularExpression::CaseInsensitiveOption);
      ui_->list_items->clear();
      ui_->list_items->addItems(items_wo_bl_.filter(regex));
    }

  } //namespace Internal
} //namspace OpenMS

