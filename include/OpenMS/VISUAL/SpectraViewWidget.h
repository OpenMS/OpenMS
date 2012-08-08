// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_VISUAL_SPECTRAVIEWWIDGET_H
#define OPENMS_VISUAL_SPECTRAVIEWWIDGET_H

#include <QWidget>
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>
#include <QtGui/QTreeWidget>

#include <OpenMS/VISUAL/LayerData.h>

namespace OpenMS
{
  /**
    @brief Hierarchical visualization and selection of spectra.

    @ingroup SpectrumWidgets
  */
  class SpectraViewWidget :
    public QWidget
  {
    Q_OBJECT
public:
    /// Constructor
    SpectraViewWidget(QWidget * parent = 0);
    /// Destructor
    virtual ~SpectraViewWidget();
    QTreeWidget * getTreeWidget();
    QComboBox * getComboBox();
    void updateEntries(const LayerData & cl);
signals:
    void spectrumSelected(int);
    void spectrumSelected(std::vector<int, std::allocator<int> > indices);
    void spectrumDoubleClicked(int);
    void showSpectrumAs1D(int);
    void showSpectrumAs1D(std::vector<int, std::allocator<int> > indices);
    void showSpectrumMetaData(int);
private:
    QLineEdit * spectra_search_box_;
    QComboBox * spectra_combo_box_;
    QTreeWidget * spectra_treewidget_;
    // cache to store mapping of chromatogram precursors to chromatogram indices
    std::map<int, std::map<Precursor, std::vector<Size>, Precursor::MZLess> > map_precursor_to_chrom_idx_cache;
private slots:
    void spectrumSelected_(const QString & text);
    void spectrumBrowserHeaderContextMenu_(const QPoint &);
    void spectrumSelectionChange_(QTreeWidgetItem *, QTreeWidgetItem *);
    void spectrumDoubleClicked_(QTreeWidgetItem *, int);
    void spectrumContextMenu_(const QPoint &);
  };
}

#endif // OPENMS_VISUAL_SPECTRAVIEWWIDGET_H
