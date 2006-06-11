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
// $Id: ClusterRunWidget.h,v 1.5 2006/03/28 08:03:26 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#ifndef OPENMS_VISUAL_CLUSTERRUNWIDGET_H
#define OPENMS_VISUAL_CLUSTERRUNWIDGET_H

#include <qdialog.h>

#include <qpushbutton.h>
#include <qcombobox.h>
#include <qlabel.h>
#include <qlayout.h>

#include <vector>

#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/MowerFunctor.h>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

namespace OpenMS
{

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

    CompareFunctor* cfp_;
    std::vector<MowerFunctor*> mowers_;
    double binsize_;
		uint binspread_;
    QPushButton* ok_;
  };
}
#endif // OPENMS_VISUAL_CLUSTERRUNWIDGET_H
