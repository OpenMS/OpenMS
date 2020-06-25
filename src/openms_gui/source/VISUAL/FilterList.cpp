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

#include <OpenMS/VISUAL/FilterList.h>
#include <ui_FilterList.h>

#include <OpenMS/VISUAL/DIALOGS/DataFilterDialog.h>

#include <QMenu>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    FilterList::FilterList(QWidget *parent) :
      QWidget(parent),
      ui_(new Ui::FilterList)
    {
      ui_->setupUi(this);
      connect(ui_->filter, &QListWidget::itemDoubleClicked, this, &FilterList::filterEdit_);
      connect(ui_->filter, &QListWidget::customContextMenuRequested, this, &FilterList::customContextMenuRequested_);
      connect(ui_->check, &QCheckBox::toggled, [&]()
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

  } //namespace Internal
} //namspace OpenMS

