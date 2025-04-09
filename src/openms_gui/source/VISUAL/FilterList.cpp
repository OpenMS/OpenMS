// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/FilterList.h>
#include <ui_FilterList.h>

#include <OpenMS/VISUAL/DIALOGS/DataFilterDialog.h>

#include <QMenu>

using namespace std;

namespace OpenMS::Internal
{

    FilterList::FilterList(QWidget *parent) :
      QWidget(parent),
      ui_(new Ui::FilterList)
    {
      ui_->setupUi(this);
      connect(ui_->filter, &QListWidget::itemDoubleClicked, this, &FilterList::filterEdit_);
      connect(ui_->filter, &QListWidget::customContextMenuRequested, this, &FilterList::customContextMenuRequested_);
      connect(ui_->check, &QCheckBox::clicked, [&]() // only on user interaction; not when calling setChecked()!
      {
        filters_.setActive(!filters_.isActive()); // invert internal representation
        emit filterChanged(filters_);             // make it public
      });
    }

    FilterList::~FilterList()
    {
      delete ui_;
    }

    void FilterList::filterEdit_(QListWidgetItem* item)
    {
      auto row = ui_->filter->row(item);
      DataFilters::DataFilter filter = filters_[row];
      DataFilterDialog dlg(filter, this);
      if (dlg.exec())
      {
        filters_.replace(row, filter);
        set(filters_);
      }
    }

    void FilterList::set(const DataFilters& filters)
    {
      filters_ = filters;

      ui_->filter->clear();
      for (Size i = 0; i < filters.size(); ++i)
      {
        QListWidgetItem* item = new QListWidgetItem(ui_->filter);
        item->setText(filters[i].toString().toQString());
      }
      // update check box
      ui_->check->setChecked(filters.isActive());

      emit filterChanged(filters_);
    }

    void FilterList::customContextMenuRequested_(const QPoint& pos)
    {
      QMenu context_menu;

      // add actions
      QListWidgetItem* item = ui_->filter->itemAt(pos);
      if (item)
      {
        context_menu.addAction("Edit", [&]() 
        {
          filterEdit_(item);
        });
        context_menu.addAction("Delete", [&]() 
        {
          filters_.remove(ui_->filter->row(item));
          set(filters_);
        });
      }
      context_menu.addAction("Add filter", [&]()
      {
        DataFilters::DataFilter filter;
        DataFilterDialog dlg(filter, this);
        if (dlg.exec())
        {
          filters_.add(filter);
          set(filters_);
        }
      });

      context_menu.exec(ui_->filter->mapToGlobal(pos));
    }
    
} //namspace OpenMS //namespace Internal

