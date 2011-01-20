// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRAIDENTIFICATIONVIEWWIDGET_H
#define OPENMS_VISUAL_SPECTRAIDENTIFICATIONVIEWWIDGET_H

#include <QWidget>
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>
#include <QtGui/QTableWidget>
#include <QtGui/QCheckBox>

#include <OpenMS/VISUAL/LayerData.h>

namespace OpenMS
{
  class SpectraIdentificationViewWidget
    : public QWidget,
      public DefaultParamHandler
  {
    Q_OBJECT
    public:
      /// Constructor
      SpectraIdentificationViewWidget(const Param& preferences, QWidget* parent = 0);
      /// Destructor
      virtual ~SpectraIdentificationViewWidget();
      /// Attach model
      void attachLayer(LayerData* model);
      /// Helper function to block outgoing signals
      bool ignore_update;
    public slots:
      /// Rebuild table entries
      void updateEntries();
    signals:
      void spectrumSelected(int);
      void spectrumDeselected(int);
      void spectrumDoubleClicked(int);
      void showSpectrumAs1D(int);
      void showSpectrumMetaData(int);
      void requestVisibleArea1D(DoubleReal, DoubleReal);
    private:
      LayerData* layer_;
      QCheckBox* hide_no_identification_;
      QCheckBox* create_rows_for_commmon_metavalue_;
      QTableWidget* table_widget_;
      bool is_ms1_shown_;
    private slots:
      /// Emits spectrumSelected with the current spectrum index
      void spectrumSelectionChange_(QTableWidgetItem*, QTableWidgetItem*);
      /// Export table entries as csv
      void exportEntries_();
      /// Saves the (potentially filtered) idXML
      void saveIdXML_();
      /// Display header context menu
      void headerContextMenu_(const QPoint&);
      /// Cell clicked in table_widget
      void cellClicked_(int row, int column);
   };
}

#endif // OPENMS_VISUAL_SPECTRAIDENTIFICATIONVIEWWIDGET_H
