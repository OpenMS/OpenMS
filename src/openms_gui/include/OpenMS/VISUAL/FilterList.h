// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/PROCESSING/MISC/DataFilters.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QWidget>

namespace Ui
{
  class FilterList;
}

class QListWidgetItem;

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief A widget which shows a list of DataFilter items.

      Filters can be added, edited and removed.
      A checkbox allows to switch them all on/off.

    */
    class FilterList : public QWidget
    {
        Q_OBJECT

    public:
      /// C'tor
      explicit FilterList(QWidget* parent);
      ~FilterList() override;

    public slots:
      /// provide new filters to the widget
      /// does invoke the 'filterChanged' signal
      void set(const DataFilters& filters);

    signals:
      /// emitted when the user has edited/added/removed a filter
      void filterChanged(const DataFilters& filters);
    
    private slots:
      /// the user wants to edit a filter (by double-clicking it)
      /// emits 'filterChanged' signal if filter was modified
      void filterEdit_(QListWidgetItem* item);

      /// right-clicking on the QListWidget 'filter' will call this slot
      void customContextMenuRequested_(const QPoint &pos);

    private:
      Ui::FilterList *ui_;
      DataFilters filters_; ///< internal representation of filters
    };
  } // ns Internal
} // ns OpenMS

// this is required to allow parent widgets (auto UIC'd from .ui) to have a FilterList member
using FilterList = OpenMS::Internal::FilterList;
