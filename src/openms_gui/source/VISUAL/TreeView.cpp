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

#include <OpenMS/VISUAL/TreeView.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QHeaderView>
#include <QMenu>

using namespace std;

///@improvement write the visibility-status of the columns in toppview.ini and read at start


namespace OpenMS
{
  TreeView::TreeView(QWidget* parent) :
    QTreeWidget(parent)
  {
    this->setObjectName("tree_widget");

    this->header()->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this->header(), &QHeaderView::customContextMenuRequested, this, &TreeView::headerContextMenu_);
  }


  void TreeView::headerContextMenu_(const QPoint& pos)
  {
    // allows to hide/show columns
    QMenu context_menu(this->header());
    const auto& header = this->headerItem();

    for (int i = 0; i < header->columnCount(); ++i)
    {
      auto action = context_menu.addAction(header->text(i), [i, this]() {
        this->setColumnHidden(i, !this->isColumnHidden(i));
        });
      action->setCheckable(true);
      action->setChecked(!this->isColumnHidden(i));
    }

    // show and execute menu
    context_menu.exec(this->mapToGlobal(pos));
  }

  void TreeView::setHeaders(const QStringList& headers)
  {
    setColumnCount(headers.size());
    setHeaderLabels(headers);
  }

  void TreeView::hideColumns(const QStringList& header_names)
  {
    auto hset = header_names.toSet();
    // add actions which show/hide columns
    const auto& header = this->headerItem();

    for (int i = 0; i < header->columnCount(); ++i)
    {
      if (hset.contains(header->text(i)))
      {
        setColumnHidden(i, true);
        hset.remove(header->text(i));
      }
    }
    if (!hset.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "header_names contains a column name which is unknown: " + String(hset.toList().join(", ")));
    }
  }

  QStringList TreeView::getHeaderNames(const WidgetHeader which) const
  {
    QStringList header_labels;
    for (int i = 0; i != columnCount(); ++i)
    {
      // do not export hidden columns
      if (which == WidgetHeader::VISIBLE_ONLY && isColumnHidden(i))
      {
        continue;
      }
      header_labels << getHeaderName(i);
    }
    return header_labels;
  }

  /// get the displayed name of the header in column with index @p header_column
  /// @throws Exception::ElementNotFound if header at index @p header_column is not valid

  QString TreeView::getHeaderName(const int header_column) const
  {
    const auto& header = this->headerItem();
    if (header->columnCount() <= header_column) throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header index " + String(header_column) + " is too large. There are only " + String(header->columnCount()) + " columns!");

    return header->text(header_column);
  }

}
