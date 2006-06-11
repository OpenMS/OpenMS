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
// $Id: ClusterExperimentView.h,v 1.8 2006/03/29 12:30:29 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer:  $
// --------------------------------------------------------------------------
#ifndef OPENMS_VISUAL_CLUSTEREXPERIMENTVIEW_H
#define OPENMS_VISUAL_CLUSTEREXPERIMENTVIEW_H

#include <qtable.h>

#include <vector>
#include <map>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{
	class DBAdapter;
	
  /**
   displays the contents of a ClusterExperiment<br>
   behaves similar to a TreeViewer <br>
   */
  class ClusterExperimentView
    :public QTable
  {
    Q_OBJECT
  public:
    ClusterExperimentView( QWidget* parent = 0 ,const char* name = 0);
    ~ClusterExperimentView();

    //uses appropriate colorgroup, if found
    void paintCell ( QPainter * p, int row, int col, const QRect & cr, bool selected , const QColorGroup & cg );
    void paintCell ( QPainter * p, int row, int col, const QRect & cr, bool selected );

    void setDBAdapter(DBAdapter* adapterp){adapterp_ = adapterp;}
  public slots:
    void doubleClick(int row, int col, int button);

    //create clex_ from file name
    void load(String name);

    //sorts clex_, not just the display
    void sort();

    //access to ClusterRuns
    const ClusterExperiment::ClusterRun& getCR(uint pos);
  signals:
    void clustering(const std::map<int,ClusterNode*>&);
    //used to show params of FactoryProduct in ClMain
    void showFactoryProduct(const FactoryProduct* const);
    //used to show result in ClMain
    void showResult(const std::map<String,double>&);
    //current ClusterRun
    void clusterRun(const ClusterExperiment::ClusterRun&);
  protected:
    //enables the user to sort
    void contextMenuEvent(QContextMenuEvent*);
    //navigation by key
    void keyPressEvent(QKeyEvent* event);
  private:
    enum RowType {DEFAULT,DATA,INFO,CLUSTERRUN,SIMFUNC,PREPROCESS,PREPROCESSFUNC,ANALYSIS,ANALYSISOBJECT,CLUSTERING,CLUSTERFUNC};

    void init_();
    void refresh_();
    void createColorGroups_();

    std::vector<const ClusterExperiment::Analysis*> getAnalysisvec_();

    uint findParentClusterRun_(uint row, RowType = DEFAULT);
    ClusterExperiment clex_;
    std::map<String,std::vector<int> > referencemap_;

    DBAdapter* adapterp_;

    //important members
    QTableItem* data_;
    QTableItem* info_;
    std::vector<QTableItem*> crun_;

    //to simulate a tree structure
    //the view has to know which items are expanded so it can close them
    bool dataexpanded_;
    bool infoexpanded_;
    std::vector<bool> crunexpanded_;
    std::vector<bool> preprocessexpanded_;
    std::vector<bool> analysisexpanded_;

    std::vector<int> oldcolumnwidth_;
    std::map<RowType,QColorGroup> colorgroups_;
    std::map<QTableItem*,RowType> rowtypes_;
  };
}
#endif //OPENMS_VISUAL_CLUSTEREXPERIMENTVIEW_H
