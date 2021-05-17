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
      blacklist_ = bl_items.toSet();
      updateInternalList_();
    }

    void FilterableList::addBlackListItems(const QStringList& items)
    {
      blacklist_.unite(items.toSet());
      updateInternalList_();
    }

    void FilterableList::removeBlackListItems(const QStringList& outdated_blacklist_items)
    {
      // quadratic runtime, but maintains order of items (as opposed to converting to set)
      for (const auto& bl : outdated_blacklist_items.toSet())
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

