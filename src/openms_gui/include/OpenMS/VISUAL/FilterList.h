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

#pragma once

#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>
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
      ~FilterList();

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
