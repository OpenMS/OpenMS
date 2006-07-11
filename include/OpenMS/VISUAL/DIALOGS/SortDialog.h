// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  $
// --------------------------------------------------------------------------
#ifndef OPENMS_VISUAL_DIALOGS_SORTDIALOG_H
#define OPENMS_VISUAL_DIALOGS_SORTDIALOG_H

#include <qdialog.h>
#include <qcombobox.h>
#include <qlabel.h>
#include <qpushbutton.h>
#include <qlayout.h>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

#include <vector>

namespace OpenMS
{

  /**
  lets the user choose an AnalysisFunctor and a measure to sort <br>
  the ClusterRun`s in a ClusterExperiment <br>
  */
  class SortDialog : public QDialog
  {
    Q_OBJECT
  public:
    SortDialog(std::vector<const ClusterExperiment::Analysis*>, QWidget* = 0, const char* = 0 );

    /** @brief return chosen Analysis <br> */
    QString analysis();
    /** @brief return chosen measure <br>*/
    QString measure();
  public slots:
    void ok();
    /** @brief look for results in the chosen Analysis <br> */
    void fillmeasure();
  private:
    /** @brief look for AnalysisFunctors to sort by <br> */
    void fillanalysis_();

    QGridLayout* gridlayout_;
    QLabel* analysislabel_;
    QLabel* measurelabel_;
    QComboBox* analysisbox_;
    QComboBox* measurebox_;
    QPushButton* ok_;

    std::vector<const ClusterExperiment::Analysis*> anafuncs_;
  };

}
#endif // OPENMS_VISUAL_SORTDIALOG_H
