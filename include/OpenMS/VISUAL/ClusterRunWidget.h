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
#ifndef OPENMS_VISUAL_CLUSTERRUNWIDGET_H
#define OPENMS_VISUAL_CLUSTERRUNWIDGET_H

#include <qdialog.h>

#include <vector>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

class QPushButton;
class QComboBox;
class QLabel;
class QLayout;
class QGridLayout;

namespace OpenMS
{
  class CompareFunctor;
  class PreprocessingFunctor;

  /**
  allows the creation of (half) ClusterRun`s<br>
  half meaning only Preprocessing and Comparison <br>
  */
  class ClusterRunWidget : public QDialog
  {
    Q_OBJECT
  public:
    ClusterRunWidget(QWidget* = 0, const char* = 0 );

    ClusterExperiment::ClusterRun* getClusterRun();
  public slots:
    void ok();
    void addpp();
    void clearpp();
    void usecf();
  private:
    void fillbox_(String type, QComboBox* box);
    void configure_(FactoryProduct* cp);

    QGridLayout* gridlayout_;
    QGridLayout* ppsublayout_;
    QLabel* cflabel_;
    QLabel* pplabel_;
    QLabel* ppnames_;

    QPushButton* addpp_;
    QPushButton* clearpp_;
    QPushButton* usecf_;

    QComboBox* cfbox_;
    QComboBox* ppbox_;

    PeakSpectrumCompareFunctor* cfp_;
    std::vector<PreprocessingFunctor*> mowers_;
    double binsize_;
		uint binspread_;
    QPushButton* ok_;
  };
}
#endif // OPENMS_VISUAL_CLUSTERRUNWIDGET_H
