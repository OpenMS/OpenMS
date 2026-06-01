// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QSet>
#include <QWidget>

namespace Ui
{
  class FilterableList;
}

class QListWidgetItem;

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief A widget which shows a list of text items, which can be filtered.

      A text field serves as filter expression (unix wildcard regEx support - 
      see https://doc.qt.io/archives/qt-5.10/qregexp.html#wildcard-matching), which hides any items in the list which do not
      match the currently typed text.

      You can also specify a blacklist of items which should never be shown, even though they are in the list.
      This is useful for only showing items which are not yet chosen elsewhere.
      Make sure all blacklist items are actually valid text items, which are contained in the list and can actually be blacklisted.
      Otherwise an exception is thrown.

      To query the currently selected items, use getSelectedItems() or connect to the 'itemDoubleClicked()' signal.

    */
    class FilterableList : public QWidget
    {
        Q_OBJECT

    public:
      /// C'tor
      explicit FilterableList(QWidget* parent);
      ~FilterableList();

      /// get the currently selected items of all visible items, i.e. must pass the filter and be selected
      QStringList getSelectedItems() const;

      /// get all items which are visible (i.e. excludes the ones which are hidden by the filter)
      QStringList getAllVisibleItems() const;

    public slots:
      /// Provide new items to the widget. 
      /// Be careful if blacklisted items have already been set. If in doubt, clear blacklisted items first to avoid an exception.
      /// @throws Exception::InvalidValue if any of the internal @p blacklist_items_ is not contained in the given @p items
      void setItems(const QStringList& items);

      /// sets items of a blacklist, i.e. they are not shown, even if they are in the item list
      /// @throws Exception::InvalidValue if any of @p blacklist_items is not contained in current items
      void setBlacklistItems(const QStringList& blacklist_items);

      /// adds items to a blacklist, i.e. they are not shown, even if they are in the item list
      /// @throws Exception::InvalidValue if any of @p additional_blacklist_items is not contained in current items
      void addBlackListItems(const QStringList& additional_blacklist_items);
     
      /// removes items from blacklist, which should not be shown, even if they are in the item list
      /// @throws Exception::InvalidValue if any of @p outdated_blacklist_items is not contained in current blacklist
      void removeBlackListItems(const QStringList& outdated_blacklist_items);
      

    signals:
      /// emitted when the user has edited the filter
      void filterChanged(const QString& filter_text);

      /// emitted when this item was double clicked
      void itemDoubleClicked(QListWidgetItem* item);
    
    private slots:
      /// internally invoked when the filter was changed by the user
      /// emits 'filterChanged' signal and updates the current list of visible items
      void filterEdited_(const QString& filter_text);

      /// recompute @p items_wo_bl_, whenever items_ or blacklist_ changed.
      /// and call updateVisibleList_()
      void updateInternalList_();
      
      /// update shown items, based on current @p items_wo_bl_ and current filter
      void updateVisibleList_();

    private:
      Ui::FilterableList* ui_;
      QStringList items_; ///< full list of items to show; when filtered only a subset is shown
      // this needs to be a Set, otherwise items might be blacklisted twice, which will lead to an Exception (because it cannot find the item on the second round)
      QSet<QString> blacklist_; ///< blacklisted items, which are never shown, even if in @p items_;
      QStringList items_wo_bl_; ///< items from @p item_ with blacklisted items removed
    };
  } // ns Internal
} // ns OpenMS

// this is required to allow parent widgets (auto UIC'd from .ui) to have a FilterableList member
using FilterableList = OpenMS::Internal::FilterableList;
