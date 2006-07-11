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
#include <OpenMS/VISUAL/ClusterExperimentView.h>

#include <OpenMS/VISUAL/DIALOGS/SortDialog.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/METADATA/Identification.h>
#include <qpopupmenu.h>
#include <qaction.h>
#include <qdialog.h>
#include <qlineedit.h>
#include <qvbox.h>
#include <qmessagebox.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace std;
namespace OpenMS
{
  ClusterExperimentView::ClusterExperimentView( QWidget* parent ,const char* name)
    :QTable(parent,name),clex_(),referencemap_(),data_(0),info_(0),crun_(),dataexpanded_(0),infoexpanded_(0),crunexpanded_(),preprocessexpanded_(),analysisexpanded_(),oldcolumnwidth_(2),colorgroups_(),rowtypes_()
  {
    setTopMargin(0);
    setLeftMargin(0);
    setShowGrid(0);
    createColorGroups_();
    init_();
    setSelectionMode(QTable::NoSelection);
    setPaletteBackgroundColor (Qt::darkGray );
    connect(this,SIGNAL(doubleClicked(int,int,int,const QPoint&)),this,SLOT(doubleClick(int,int,int)));
  }

  void ClusterExperimentView::createColorGroups_()
  {
    if (colorgroups_.size()) return;
    else
    {
      QColor foreground = Qt::black; //for example text in QLabel
      QColor button = Qt::gray; //background on buttons
      QColor light = Qt::lightGray; //3D effects like on the side of buttons
      QColor dark = Qt::darkGray; //see above (light)
      QColor mid = Qt::gray; // dont know, dont care, most of these colors are not used in ClusterExperimentView, but QColorGroup wants it anyway
      QColor text = Qt::white; //text color
      QColor brighttext = Qt::black; // inverted text color
      QColor base = Qt::darkGray; //background for text widgets
      QColor background = Qt::darkGray; // actually just base and text are used in this Class. So the other Values might not make much sense
      QColor lightblue = QColor(200,200,255);
      QColor mediumblue = QColor(180,180,255);
      colorgroups_.insert(make_pair(DEFAULT,        QColorGroup(foreground,button,light,dark,mid, text,             brighttext,base,          background)));
      colorgroups_.insert(make_pair(DATA,           QColorGroup(foreground,button,light,dark,mid, text,             brighttext,base,          background)));
      colorgroups_.insert(make_pair(INFO,           QColorGroup(foreground,button,light,dark,mid, text,             brighttext,base,          background)));
      colorgroups_.insert(make_pair(CLUSTERRUN,     QColorGroup(foreground,button,light,dark,mid, /*text*/lightblue,brighttext,base,          background)));
      //DEFAULT,DATA,INFO,CLUSTERRUN,SIMFUNC,PREPROCESS,ANALYSIS,CLUSTERING,CLUSTERFUNC
      colorgroups_.insert(make_pair(SIMFUNC,        QColorGroup(foreground,button,light,dark,mid, Qt::blue,        brighttext,Qt::gray,      background)));
      colorgroups_.insert(make_pair(CLUSTERFUNC,    QColorGroup(foreground,button,light,dark,mid, Qt::blue,        brighttext,Qt::gray,      background)));
      colorgroups_.insert(make_pair(PREPROCESS,     QColorGroup(foreground,button,light,dark,mid, Qt::blue,          brighttext,Qt::gray,      background)));
      colorgroups_.insert(make_pair(PREPROCESSFUNC, QColorGroup(foreground,button,light,dark,mid, Qt::darkBlue,      brighttext,Qt::lightGray, background)));
      colorgroups_.insert(make_pair(ANALYSIS,       QColorGroup(foreground,button,light,dark,mid, Qt::blue,      brighttext,Qt::gray,      background)));
      colorgroups_.insert(make_pair(ANALYSISOBJECT, QColorGroup(foreground,button,light,dark,mid, Qt::darkBlue,      brighttext,Qt::lightGray, background)));
      colorgroups_.insert(make_pair(CLUSTERING,     QColorGroup(foreground,button,light,dark,mid, Qt::blue,         brighttext,Qt::gray,      background)));
    }
  }

  ClusterExperimentView::~ClusterExperimentView()
  {
  }

  void ClusterExperimentView::load(String name)
  {
    clex_.load(name);
    refresh_();
  }

  void ClusterExperimentView::refresh_()
  {
    init_();
    setNumRows(2+clex_.size());
    QTableItem* temp;
    rowtypes_.clear();
    crun_.clear();
    crun_.resize(clex_.size());
    crunexpanded_.clear();
    crunexpanded_.resize(clex_.size());
    preprocessexpanded_.clear();
    preprocessexpanded_.resize(clex_.size());
    analysisexpanded_.clear();
    analysisexpanded_.resize(clex_.size());
    rowtypes_[item(0,0)] = DATA;
    rowtypes_[item(1,0)] = INFO;
    for (uint i = 0; i < clex_.size(); ++i)
    {
      temp = new QTableItem(this,QTableItem::Never,"ClusterRun");
      crun_[i] = temp;
      setItem(i+2,0,temp);
      //the ClusterRun Numbers are not just pretty, they do something usefull, too
      //so if you change the text of those Cells, doubleClick needs to be changed as well
      setItem(i+2,1,new QTableItem(this,QTableItem::Never,QString("%1").arg(i)));
      rowtypes_[item(i+2,0)] = CLUSTERRUN;
    }
  }

  void ClusterExperimentView::init_()
  {
    setNumRows(0);
    setNumCols(0);
    setNumRows(2);
    setNumCols(2);
    data_ = new QTableItem(this,QTableItem::Never,"Data");
    info_ = new QTableItem(this,QTableItem::Never,"Info");
    setItem(0,0,data_);
    setItem(1,0,info_);
    setItem(0,1,new QTableItem(this,QTableItem::Never,clex_.datasetname().c_str()));
    setItem(1,1,new QTableItem(this,QTableItem::Never,clex_.infodate().c_str()));
    rowtypes_[item(0,0)] = DATA;
    rowtypes_[item(1,0)] = INFO;
  }

  vector<const ClusterExperiment::Analysis*> ClusterExperimentView::getAnalysisvec_()
  {
    map<String,const ClusterExperiment::Analysis*> differentanas;
    for ( uint i = 0; i < clex_.size(); ++i )
    {
      for ( uint j = 0; j < clex_[i].size(); ++j )
      {
        const ClusterExperiment::Analysis* ana = &dynamic_cast<const ClusterExperiment::Analysis&>(clex_[i][j]);
        if ( differentanas.find( ana->name() ) == differentanas.end() )
        {
          differentanas.insert(make_pair(ana->name(),ana));
        }
      }
    }
    vector<const ClusterExperiment::Analysis*> result;
    for ( map<String,const ClusterExperiment::Analysis*>::const_iterator cmit = differentanas.begin(); cmit != differentanas.end(); ++cmit )
    {
      result.push_back(cmit->second);
    }
    return result;
  }

  void ClusterExperimentView::sort()
  {
    vector<const ClusterExperiment::Analysis*> anas = getAnalysisvec_();
    if ( anas.size() )
    {
      SortDialog* bla = new SortDialog(anas,this);
      
      bla->exec();
      bla->raise();
      bla->setActiveWindow();
      clex_.sortbyResult(bla->analysis().ascii(),bla->measure().ascii());
      refresh_();
    }
    else
    {
      QMessageBox::critical(this,"OpenMS_Cat","there seems to be no Analysis in this ClusterExperiment");
    }
  }

  void ClusterExperimentView::contextMenuEvent(QContextMenuEvent* event)
  {
    if ( ! clex_.size() )
    {
      return;
    }
    QPopupMenu contextMenu(this);
    QAction* act1 = new QAction("Sort ClusterExperiment","sort by",0,this); 
    act1->addTo(&contextMenu);
    connect(act1, SIGNAL( activated() ) , this, SLOT( sort() ) );
    contextMenu.exec(event->globalPos());
  }

  void ClusterExperimentView::doubleClick(int row, int col, int /*button*/)
  {
    //check what was doubleclicked
    if ( !clex_.size() ) return;
    int crnr = 0;
    switch (rowtypes_[item(row,col)])
    {
      case DATA:
        crnr = -1;
        if (!dataexpanded_)
        {
          if (clex_.datasetname().length())
          {
            insertRows(row+1);
            setItem(row+1,0,new QTableItem(this,QTableItem::Never,"  datasetname"));
            setItem(row+1,1,new QTableItem(this,QTableItem::Never,clex_.datasetname().c_str()));
            dataexpanded_ = true;
          }
          if (clex_.datasetid() > 0 )
          {
            insertRows(row+1);
            setItem(row+1,0,new QTableItem(this,QTableItem::Never,"  datasetid"));
            setItem(row+1,1,new QTableItem(this,QTableItem::Never,QString("%1").arg(clex_.datasetid())));
            rowtypes_[item(row+1,0)] = DATA;
            dataexpanded_ = true;
          }
          if (clex_.datasetsize() > 0 )
          {
            insertRows(row+1);
            setItem(row+1,0,new QTableItem(this,QTableItem::Never,"  datasetsize"));
            setItem(row+1,1,new QTableItem(this,QTableItem::Never,QString("%1").arg(clex_.datasetsize())));
            dataexpanded_ = true;
          }
        }
        else
        {
          while (data_->row() + 1 < info_->row() )
          {
            rowtypes_.erase(rowtypes_.find(item(row+1,0)));
            removeRow(row+1);
          }
          dataexpanded_ = false;
        }
        break;
      case INFO:
        crnr = -1;
        if (!infoexpanded_)
        {
          if (clex_.infouser().length() )
          {
            insertRows(row+1);
            setItem(row+1,0,new QTableItem(this,QTableItem::Never," user"));
            setItem(row+1,1,new QTableItem(this,QTableItem::Never,clex_.infouser().c_str()));
            infoexpanded_ = true;
          }
          if (clex_.infodate().length() )
          {
            insertRows(row+1);
            setItem(row+1,0,new QTableItem(this,QTableItem::Never," date"));
            setItem(row+1,1,new QTableItem(this,QTableItem::Never,clex_.infodate().c_str()));
            infoexpanded_ = true;
          }
          if (clex_.infocomment().length() )
          {
            insertRows(row+1);
            setItem(row+1,0,new QTableItem(this,QTableItem::Never," comment"));
            setItem(row+1,1,new QTableItem(this,QTableItem::Never,clex_.infocomment().c_str()));
            setRowHeight(row+1,item(row+1,1)->sizeHint().height());
            infoexpanded_ = true;
          }
        }
        else
        {
          while ( info_->row() + 1 < crun_[0]->row() )
          {
            rowtypes_.erase(rowtypes_.find(item(row+1,0)));
            removeRow(row+1);
          }
          infoexpanded_ = false;
        }
        break;
      case CLUSTERRUN:
        {
          crnr = item(row,1)->text().toInt();
          if (crun_[crnr]->row() == row)
          {
            if (!crunexpanded_[crnr])
            {
              insertRows(row+1,5);
              setItem(row+1,0,new QTableItem(this,QTableItem::Never,"    Preprocessing"));
              setItem(row+1,1,new QTableItem(this,QTableItem::Never,QString(": %1").arg(clex_[crnr].getPreprocessqueue().size())));
              setItem(row+2,0,new QTableItem(this,QTableItem::Never,"    CompareFunctor"));
              setItem(row+2,1,new QTableItem(this,QTableItem::Never,clex_[crnr].getSimFunc()->getName().c_str()));
              setItem(row+3,0,new QTableItem(this,QTableItem::Never,"    ClusterFunctor"));
              setItem(row+3,1,new QTableItem(this,QTableItem::Never,clex_[crnr].getClusterFunc()->getName().c_str()));
              setItem(row+4,0,new QTableItem(this,QTableItem::Never,"    Clustering"));
              setItem(row+4,1,new QTableItem(this,QTableItem::Never,QString(": %1").arg(clex_[crnr].getClustering().size())));
              setItem(row+5,0,new QTableItem(this,QTableItem::Never,"    Analysis"));
              setItem(row+5,1,new QTableItem(this,QTableItem::Never,QString(": %1").arg(clex_[crnr].size())));
              rowtypes_[item(row+1,0)] = PREPROCESS;
              rowtypes_[item(row+2,0)] = SIMFUNC;
              rowtypes_[item(row+3,0)] = CLUSTERFUNC;
              rowtypes_[item(row+4,0)] = CLUSTERING;
              rowtypes_[item(row+5,0)] = ANALYSIS;
              crunexpanded_[crnr] = true;
            }
            else
            {
              uint collapsecount = 0;
              // 5 are always there + the analysis + the Filters
              collapsecount +=5;
              if (analysisexpanded_[crnr]) collapsecount += clex_[crnr].size();
              if (preprocessexpanded_[crnr]) collapsecount += clex_[crnr].getPreprocessqueue().size();
              for (uint j = 0; j <  collapsecount ;++j)
              {
                rowtypes_.erase(rowtypes_.find(item(row+1,0)));
                removeRow(row+1);
              }
              crunexpanded_[crnr] = false;
              analysisexpanded_[crnr] = false;
              preprocessexpanded_[crnr] = false;
            }
          }
        }
        break;
      case CLUSTERING:
        {
          //find the corresponding ClusterRun
          crnr = findParentClusterRun_(row,CLUSTERING);
          //emit bins(clex_[crnr].getSimFunc()->usebins());
          emit clustering(clex_[crnr].getClustering());
        }
        break;
      case SIMFUNC:
        {
          crnr = findParentClusterRun_(row,SIMFUNC);
          emit showFactoryProduct(clex_[crnr].getSimFunc());
          break;
        }
      case CLUSTERFUNC:
        {
          crnr = findParentClusterRun_(row,CLUSTERFUNC);
					// TODO 
          //emit showFactoryProduct(clex_[crnr].getClusterFunc());
          break;
        }
      case PREPROCESS:
        {
          crnr = findParentClusterRun_(row,PREPROCESS);
          if (!preprocessexpanded_[crnr])
          {
            insertRows(row+1,clex_[crnr].getPreprocessqueue().size());
            //we go backwards through the container to have it in the right order in the table
            for (int i = clex_[crnr].getPreprocessqueue().size()-1; i >= 0; --i)
            {
              setItem(row+1,0,new QTableItem(this,QTableItem::Never,"      PreprocessFunctor"));
              setItem(row+1,1,new QTableItem(this,QTableItem::Never,clex_[crnr].getPreprocessqueue()[i]->getName().c_str()));
              rowtypes_[item(row+1,0)] = PREPROCESSFUNC; 
            }
            preprocessexpanded_[crnr] = true;
          }
          else 
          {
            for (uint i = 0; i < clex_[crnr].getPreprocessqueue().size(); ++i)
            {
              rowtypes_.erase(rowtypes_.find(item(row+1,0)));
              removeRow(row+1);
            }
            preprocessexpanded_[crnr] = false; 
          }
          break;
        }
      case ANALYSIS:
        {
          crnr = findParentClusterRun_(row,ANALYSIS);
          if (!analysisexpanded_[crnr])
          {
            insertRows(row+1,clex_[crnr].size());
            //as above we want the same order as in the Object, although the order of analysis is not important, contrary to preprocessing
            for (int i = clex_[crnr].size()-1; i >= 0; --i)
            {
              setItem(row+1,0,new QTableItem(this,QTableItem::Never,("      " + clex_[crnr][i].anafuncp()->getName()).c_str()));
              setItem(row+1,1,new QTableItem(this,QTableItem::Never,QString("%1").arg(i)));
              rowtypes_[item(row+1,0)] = ANALYSISOBJECT;
            }
            analysisexpanded_[crnr] = true;
          }
          else
          {
            for (uint i = 0; i < clex_[crnr].size(); ++i)
            {
              rowtypes_.erase(rowtypes_.find(item(row+1,0)));
              removeRow(row+1);
            }
            analysisexpanded_[crnr] = false; 
          }
          break;
        }
      case ANALYSISOBJECT:
      {
        crnr = findParentClusterRun_(row,ANALYSISOBJECT);
        emit showFactoryProduct(clex_[crnr][item(row,col)->text().toInt()].anafuncp());
        emit showResult(clex_[crnr][item(row,col)->text().toInt()].results());
      }
      case DEFAULT:
        {
          //noop
          break;
        }
      default:
        {
          cerr << "unrecognized QTableItem!\n";
          //todo? throw Exception?
        }
    }
    emit(clusterRun(clex_[crnr]));
    adjustColumn(0);
    adjustColumn(1);
  }

  // to supress warning "only partially overridden"
  // uses colorgroups_.begin(), will not work at empty cologroups_
  void ClusterExperimentView::paintCell ( QPainter * p, int row, int col, const QRect & cr, bool selected )
  {
    if ( colorgroups_.size() )
    {
      paintCell(p,row,col,cr,selected,colorgroups_.begin()->second);
    }
    else
    {
      throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"cologroups_ is empty","don't know what colors to use");
    }
  }

  void ClusterExperimentView::paintCell ( QPainter * p, int row, int col, const QRect & cr, bool selected , const QColorGroup & cg )
  {
    map<RowType,QColorGroup>::iterator customcolorit = colorgroups_.find(rowtypes_[item(row,0)]);
    if ( customcolorit != colorgroups_.end())
    {
      QTable::paintCell(p,row,col,cr,selected,customcolorit->second);
    }
    else
    {
      QTable::paintCell(p,row,col,cr,selected,cg);
    }
  }

  //might make trouble if used on the rows before the clusterruns
  //tood? check?
  uint ClusterExperimentView::findParentClusterRun_(uint row, RowType type)
  {
    uint j;
    //find the minimum distance to the parent ClusterRun Row
    //todo fill in Values
    switch (type)
    {
      case DEFAULT:
        //no break
      default:
        {
          j = 0;
        }
    }
    for (; rowtypes_[item(row-j,0)] != CLUSTERRUN ;++j);
    return item(row-j,1)->text().toInt();
  }

  void ClusterExperimentView::keyPressEvent(QKeyEvent* event)
  {
    switch (event->key())
    {
    case Qt::Key_Prior:
      {
        if ( crun_.size() )
        {
          if ( currentRow() >= crun_[1]->row() )
          {
            uint crnr = findParentClusterRun_(currentRow());
            setCurrentCell(crun_[crnr-1]->row(),currentColumn());
          }
        }
        break;
      }
    case Qt::Key_Next:
      {
        if ( crun_.size() )
        {
          if( currentRow() >= crun_[0]->row() && currentRow() < crun_.back()->row())
          {
            uint crnr = findParentClusterRun_(currentRow());
            setCurrentCell(crun_[crnr+1]->row(),currentColumn());
          }
        }
        break;
      }
    default:
      {
        QTable::keyPressEvent(event);
      }
    }
  }

  const ClusterExperiment::ClusterRun& ClusterExperimentView::getCR(uint pos)
  {
    return clex_[pos];
  }
}
