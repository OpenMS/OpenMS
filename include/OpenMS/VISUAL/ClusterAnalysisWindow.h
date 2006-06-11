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
// $Id: ClusterAnalysisWindow.h,v 1.10 2006/06/06 12:53:32 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#ifndef OPENMS_VISUAL_CLUSTERANALYSISWINDOW_H
#define OPENMS_VISUAL_CLUSTERANALYSISWINDOW_H

#include <qwidget.h>
#include <qsplitter.h>
#include <qlayout.h>
#include <qvbox.h>
#include <qmainwindow.h>
#include <qlistview.h>
#include <qmenubar.h>
#include <qstatusbar.h>
#include <qtabwidget.h>

#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>
#include <OpenMS/VISUAL/BinnedRepWidget.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/VISUAL/DIALOGS/DBDialog.h>
#include <OpenMS/VISUAL/SimMatrixWidget.h>
#include <OpenMS/VISUAL/ScaleWidget.h>
#include <OpenMS/VISUAL/FactoryProductView.h>
#include <OpenMS/VISUAL/ResultView.h>
#include <OpenMS/VISUAL/ClusterExperimentView.h>

#include <vector>
#include <map>

namespace OpenMS
{
	class DBAdapter;
	
  /**
  	@brief main window for analyzing tandem ms spectrum clusters
  
  	 
  */
  class ClusterAnalysisWindow : public QVBox
  {
    Q_OBJECT
  public:
    ClusterAnalysisWindow(QWidget* = 0, const char* = 0);
    ~ClusterAnalysisWindow();
  public slots:
    /** @brief display spectra <br>*/
    void displaySpectra(int,int);
    /** @brief show Cluster <br> */
    void showCluster(const std::vector<int>& ,QListViewItem* = 0);
    /** @brief add single ClusterSpectrum to cache */
    void addClusterSpectrum(const ClusterSpectrum& ,bool = 1);
    /** @brief display cluster in SimMatrixWidget or highlight peptide <br> */
    void clusterListViewDoubleClick(QListViewItem*);
    /** @brief display two spectra <br> */
    void inspect(int = 0,int = 0,bool = 0, bool = 0);
    /** @brief show spectrum information in the status bar <br> */
    void updateStatusBar(int,int);
    /** @brief let user input DB connection <br> */
    void connect2DB();
    /** @brief let user specify location of ClusterExperiment xml <br> */
    void loadClusterExperiment();
    /** @brief show clustering in clusterlistview_ <br> */
    void showClustering(const std::map<int,ClusterNode*>& clustering);
    /** @brief raise dialog to specify Preprocessing and Comparison functors <br> */
    void choosefunctors();
    /** @brief set ClusterRun containing Preprocessing and Comparison functors <br> */
    void getClusterRun(const ClusterExperiment::ClusterRun&);
  private:
    void createLayout_();
    void createWidgets_();
    void init_();
    void createMatrix_();
    void doLayout_();
    void connect_();
    void fillspecListView_(int,int);

    QTabWidget* infotab_;
    QLabel* specinfoview_;
    QListView* clusterlistview_;
    QVBox* mainbox_;
    QWidget* centralwidget_;
    QWidget* rightwidget_;
    QMenuBar* menubar_;
    QSplitter* mainsplit_;
    QSplitter* rightsplit_;
    QVBox* leftbox_;
    QVBox* specbox_;
    QHBox* infobox_;
    QStatusBar* statusbar_;
    QLabel* position_;
    QLabel* info_;
    QLabel* score_;
    QBoxLayout* mainlayout_;
    QBoxLayout* rightlayout_;
    DBDialog* dbdialog_;

    //data
    std::vector<std::vector<double > >* current_matrix_;
    std::vector<OpenMS::ClusterSpectrum> current_cspectra_;
    std::vector<int> current_ids_;

    const std::map<int,ClusterNode*>* clusterp_;
    std::vector<int> ids_;

    ClusterExperiment::ClusterRun clusterrun_;

    SimMatrixWidget* clmatrix_;
    ScaleWidget* scale_;
    ClusterExperimentView* clev_;
    FactoryProductView* cfigview_;
    ResultView* resview_;
    OpenMS::BinnedRepWidget* topbin_;
    OpenMS::BinnedRepWidget* bottombin_;
    OpenMS::Spectrum1DWidget* topspec_;
    OpenMS::Spectrum1DWidget* bottomspec_;
    OpenMS::DBAdapter* adapter_;
  };
}
#endif // OPENMS_VISUAL_CLUSTERANALYSISWINDOW_H
