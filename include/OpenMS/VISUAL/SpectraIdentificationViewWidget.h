// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
      void spectrumDoubleClicked(int);
      void showSpectrumAs1D(int);
      void showSpectrumMetaData(int);
    private:
      LayerData* layer_;
      QCheckBox* hide_no_identification_;
      QCheckBox* hide_ms1_;
      QTableWidget* table_widget_;
    private slots:
      void spectrumSelectionChange_(QTableWidgetItem*, QTableWidgetItem*);
   };
}

#endif // OPENMS_VISUAL_SPECTRAIDENTIFICATIONVIEWWIDGET_H
