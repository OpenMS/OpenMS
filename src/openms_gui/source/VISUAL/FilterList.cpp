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

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/DIALOGS/DataFilterDialog.h>


using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    FilterList::FilterList(QWidget *parent, TOPPViewBase* const base) :
      QWidget(parent),
      ui_(new Ui::FilterList),
      base_(base)
    {
      ui_->setupUi(this);
      connect(ui_->filter, &QListWidget::itemDoubleClicked, this, &FilterList::filterEdit);
      connect(ui_->filter, &QListWidget::customContextMenuRequested, this, &FilterList::customContextMenuRequested);
      connect(ui_->check, &QCheckBox::toggled, base_, &TOPPViewBase::layerFilterVisibilityChange);
    }

    FilterList::~FilterList()
    {
      delete ui_;
    }

    void FilterList::filterEdit(QListWidgetItem* item)
    {
      auto row = ui_->filter->row(item);
      DataFilters filters = base_->getActiveCanvas()->getCurrentLayer().filters;
      DataFilters::DataFilter filter = filters[row];
      DataFilterDialog dlg(filter, this);
      if (dlg.exec())
      {
        filters.replace(row, filter);
        base_->getActiveCanvas()->setFilters(filters);
        update(filters);
      }
    }

    void FilterList::update(const DataFilters& filters)
    {
      ui_->filter->clear();
      for (Size i = 0; i < filters.size(); ++i)
      {
        QListWidgetItem* item = new QListWidgetItem(ui_->filter);
        item->setText(filters[i].toString().toQString());
      }
      // update check box
      ui_->check->setChecked(filters.isActive());
    }

    void FilterList::customContextMenuRequested(const QPoint& pos)
    {
      SpectrumCanvas* canvas = base_->getActiveCanvas();
      // do nothing if no window is open
      if (canvas == nullptr)
        return;

      // do nothing if no layer is loaded into the canvas
      if (canvas->getLayerCount() == 0)
        return;

      QMenu context_menu;

      // warn if the current layer is not visible
      String layer_name = String("Layer: ") + base_->getActiveCanvas()->getCurrentLayer().name;
      if (!canvas->getCurrentLayer().visible)
      {
        layer_name += " (invisible)";
      }
      context_menu.addAction(layer_name.toQString())->setEnabled(false);
      context_menu.addSeparator();

      // add actions
      QListWidgetItem* item = ui_->filter->itemAt(pos);
      if (item)
      {
        context_menu.addAction("Edit");
        context_menu.addAction("Delete");
      }
      else
      {
        context_menu.addAction("Add filter");
      }
      
      // results
      QAction* selected = context_menu.exec(ui_->filter->mapToGlobal(pos));
      
      if (selected == nullptr) return;

      if (selected->text() == "Delete")
      {
        DataFilters filters = canvas->getCurrentLayer().filters;
        filters.remove(ui_->filter->row(item));
        canvas->setFilters(filters);
        update(filters);
      }
      else if (selected->text() == "Edit")
      {
        filterEdit(item);
      }
      else if (selected->text() == "Add filter")
      {
        DataFilters::DataFilter filter;
        DataFilterDialog dlg(filter, this);
        if (dlg.exec())
        {
          DataFilters filters = canvas->getCurrentLayer().filters;
          filters.add(filter);
          canvas->setFilters(filters);
          update(filters);
        }
      }
    }

  } //namespace Internal
} //namspace OpenMS

