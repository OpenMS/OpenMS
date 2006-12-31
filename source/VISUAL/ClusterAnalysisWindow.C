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

#include <OpenMS/VISUAL/ClusterAnalysisWindow.h>
#include <OpenMS/VISUAL/BinnedRepWidget.h>
#include <OpenMS/VISUAL/SimMatrixWidget.h>
#include <OpenMS/VISUAL/ScaleWidget.h>
#include <OpenMS/VISUAL/FactoryProductView.h>
#include <OpenMS/VISUAL/ResultView.h>
#include <OpenMS/VISUAL/ClusterExperimentView.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>
#include <OpenMS/FORMAT/DBAdapter.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/VISUAL/ClusterRunWidget.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/DIALOGS/InspectDialog.h>
#include <OpenMS/COMPARISON/CLUSTERING/SpectrumGenerator.h>

#include <qpixmap.h>
#include <qmenubar.h>
#include <qsplitter.h>
#include <qtabwidget.h>
#include <qlabel.h>
#include <qlistview.h>
#include <qstatusbar.h>
#include <qpopupmenu.h>
#include <qprogressdialog.h>
#include <qfiledialog.h>
#include <qstylesheet.h>
#include <qinputdialog.h>
#include <qmessagebox.h>
#include <qaction.h>
#include <qcursor.h>

#include <iostream>
#include <sstream>

#include <qapplication.h>

using namespace std;
namespace OpenMS
  {
  ClusterAnalysisWindow::ClusterAnalysisWindow(QWidget* parent, const char* name)
    : QVBox(parent,name),dbdialog_(0),current_matrix_(0),clusterp_(0),adapter_(0)
  {
    createLayout_();
    createWidgets_();
    doLayout_();
    init_();
    connect_();
    
  }

  ClusterAnalysisWindow::~ClusterAnalysisWindow()
  {
  }

  void ClusterAnalysisWindow::createLayout_()
  {
    menubar_ = new QMenuBar(this);
    centralwidget_ = new QWidget(this);
    mainlayout_ = new QBoxLayout(centralwidget_,QBoxLayout::LeftToRight);
    mainsplit_ = new QSplitter(centralwidget_);
    leftbox_ = new QVBox(mainsplit_);
    rightwidget_ = new QWidget(mainsplit_);
    rightsplit_ = new QSplitter(rightwidget_);
    rightsplit_->setOrientation(Qt::Vertical);
    rightlayout_ = new QBoxLayout(rightwidget_,QBoxLayout::TopToBottom);
    specbox_ = new QVBox(rightsplit_);
    infobox_ = new QHBox(rightsplit_);
  }

  void ClusterAnalysisWindow::createWidgets_()
  {
    clmatrix_ = new SimMatrixWidget(leftbox_);
    scale_ = new ScaleWidget(leftbox_);
    infotab_ = new QTabWidget(infobox_);
    specinfoview_ = new QLabel(infotab_);
    infotab_->addTab(specinfoview_,"Spectra");
    cfigview_ = new FactoryProductView(infotab_);
    infotab_->addTab(cfigview_,"Functor Details");
    resview_ = new ResultView(infotab_);
    infotab_->addTab(resview_,"Result");
    clev_ = new ClusterExperimentView(leftbox_);
    clusterlistview_ = new QListView(infobox_);
    topbin_ = new BinnedRepWidget(specbox_);
    bottombin_ = new BinnedRepWidget(specbox_);
    topspec_ = new Spectrum1DWidget(specbox_);
    bottomspec_ = new Spectrum1DWidget(specbox_); 
    statusbar_ = new QStatusBar(this);
    position_ = new QLabel("Position          ",this);
    info_ = new QLabel("Information",this);
    score_ = new QLabel("Score   ",this);
  }

  void ClusterAnalysisWindow::doLayout_()
  {
    clusterlistview_->addColumn("Nr");
    clusterlistview_->addColumn("info");
    clusterlistview_->addColumn("reference");
    clusterlistview_->setRootIsDecorated(true);
    leftbox_->setSpacing(5);
    mainlayout_->add(mainsplit_);
    mainsplit_->setResizeMode(leftbox_,QSplitter::KeepSize);
    bottombin_->draw_upside_down();
    //bottomspec_->setMirroredYAxis(1);
    topspec_->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
    bottomspec_->canvas()->setDrawMode(Spectrum1DCanvas::DM_PEAKS);
    topspec_->canvas()->showGridLines(0);
    bottomspec_->canvas()->showGridLines(0);
    rightlayout_->add(rightsplit_);
    QPopupMenu* loadmenu_ = new QPopupMenu(this);
    menubar_->insertItem("connectDB",this,SLOT(connect2DB()),0,1);
    menubar_->insertItem("&load",loadmenu_,2);
    loadmenu_->insertItem("load ClusterExperiment",this,SLOT(loadClusterExperiment()));
    menubar_->insertItem("inspect",this,SLOT(inspect()),0,3);
    menubar_->insertItem("Functors",this,SLOT(choosefunctors()),0,5);
    menubar_->setItemEnabled(1,1);
    menubar_->setItemEnabled(2,0);
    menubar_->setItemEnabled(3,0);
    menubar_->setItemEnabled(4,0);
    menubar_->setItemEnabled(5,0);
    this->setStretchFactor(centralwidget_,1);
    this->setStretchFactor(statusbar_,0);
    statusbar_->addWidget(position_);
    statusbar_->addWidget(info_,1);
    statusbar_->addWidget(score_);
    score_->setAlignment(AlignHCenter);
    score_->setMinimumSize(score_->sizeHint());
    position_->setAlignment(AlignHCenter);
    position_->setMinimumSize(position_->sizeHint());
    QValueList<int> rightsizes;
    rightsizes.append(500);
    rightsizes.append(140);
    rightsplit_->setSizes(rightsizes);
  }

  /**
  \param x index of first spectrum
  \param y index of second spectrum
   */
  void ClusterAnalysisWindow::updateStatusBar(int x, int y)
  {
    int xid, yid;
    xid = current_cspectra_[x].id();
    yid = current_cspectra_[y].id();
    QString pos = QString("X: %1 Y: %2").arg(xid).arg(yid);
    position_->setText(pos);
    QString seq1 = current_cspectra_[x].getTophit().getSequence().c_str();
    QString seq2 = current_cspectra_[y].getTophit().getSequence().c_str();
    QString scor = QString("Score %1").arg((*current_matrix_)[x][y]);
    score_->setText(scor);
    QString seqs = QString("X = %1  Y = %2").arg(seq1).arg(seq2);
    info_->setText(seqs);
  }

  void ClusterAnalysisWindow::connect_()
  {
    connect(topspec_,SIGNAL(visibleAreaChanged(double,double)),bottomspec_,SLOT(setVisibleArea(double,double)));
    connect(clmatrix_,SIGNAL(matrixEntry(int, int)),this,SLOT(displaySpectra(int,int)));
    connect(clmatrix_,SIGNAL(score(double)),scale_,SLOT(setMarker(double)));
    connect(clmatrix_,SIGNAL(info(int, int)),this,SLOT(updateStatusBar(int,int)));
    connect(topbin_,SIGNAL(boxInfo(const QString&)),score_,SLOT(setText(const QString&)));
    connect(topbin_,SIGNAL(posInfo(const QString&)),position_,SLOT(setText(const QString&)));
    connect(topbin_,SIGNAL(binInfo(const QString&)),info_,SLOT(setText(const QString&)));
    connect(bottombin_,SIGNAL(boxInfo(const QString&)),score_,SLOT(setText(const QString&)));
    connect(bottombin_,SIGNAL(posInfo(const QString&)),position_,SLOT(setText(const QString&)));
    connect(bottombin_,SIGNAL(binInfo(const QString&)),info_,SLOT(setText(const QString&)));
    connect(clusterlistview_,SIGNAL(doubleClicked(QListViewItem*)),this,SLOT(clusterListViewDoubleClick(QListViewItem*)));
    connect(clev_,SIGNAL(clustering(const std::map<int,ClusterNode*>&)),this,SLOT(showClustering(const std::map<int,ClusterNode*>&)));
    connect(clev_,SIGNAL(showFactoryProduct(const FactoryProduct* const)),cfigview_,SLOT(displayFactoryProduct(const FactoryProduct* const)));
    connect(clev_,SIGNAL(showResult(const std::map<std::string,double>&)),resview_,SLOT(displayResults(const std::map<std::string,double>&)));
    connect(clev_,SIGNAL(clusterRun(const ClusterExperiment::ClusterRun&)),this,SLOT(getClusterRun(const ClusterExperiment::ClusterRun&)));
  }

  /**
   \param cr is copied to clusterrun_
   */
  void ClusterAnalysisWindow::getClusterRun(const ClusterExperiment::ClusterRun& cr)
  {
    clusterrun_ = cr;
  }

  void ClusterAnalysisWindow::init_()
  {
  	MSExperiment<> exp;
  	exp.resize(1);
  	
  	for ( uint i = 0; i < 10; ++i )
    {
      DPeak<1> p;
      p.getPosition()[0] = i*10;
      p.getIntensity() = i;
      exp[0].push_back(p);
    }
    topbin_->loadSpectrum(exp[0]);
    bottombin_->loadSpectrum(exp[0]);
    
    topspec_->canvas()->addLayer(exp);
    bottomspec_->canvas()->addLayer(exp); 
  }


  void ClusterAnalysisWindow::createMatrix_()
  {
    uint size = current_cspectra_.size();
    QProgressDialog* progress = new QProgressDialog("computing Simliaritymatrix","cancel",size,this,0,true);
    progress->setCancelButton(new QPushButton("cancel",progress));
    vector<double> autocorr(size);
    if ( current_matrix_ ) delete current_matrix_;
    current_matrix_ = new vector<vector<double> >(size,vector<double>(size));
    const PeakSpectrumCompareFunctor* cfp = clusterrun_.getSimFunc();
    if ( !cfp )
    {
      QMessageBox::critical(this,"error while creating similarity matrix","no compare and preprocessfunctions found\nclick on a clusterrun or enter them yourself");
      return;
    }
    for ( uint i = 0; i < size; ++i )
    {
      //autocorr[i] = (*cfp)(current_cspectra_[i],current_cspectra_[i]);
    }
    for ( uint i = 0; i < size; ++i )
    {
      progress->setProgress(i+1);
      if (progress->wasCancelled()) return;
      for ( uint j = 0; j < size; ++j )
      {
        double score = clusterrun_.similarity(current_cspectra_[i],current_cspectra_[j],autocorr[i],autocorr[j]);
        (*current_matrix_)[i][j] = score; 
      }
    }
  }

  /**
  /param x index (in current_cspectra_) of first spectrum <br>
  /param y index (in current_cspectra_) of second spectrum <br>
  */
  void ClusterAnalysisWindow::displaySpectra(int x, int y)
  {
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    fillspecListView_(x,y);
    
		if ( clusterrun_.getSimFunc()/* && !clusterrun_.getSimFunc()->usebins()*/ )
    {
      topbin_->hide();
      bottombin_->hide();
      topspec_->show();
      bottomspec_->show();
      
	    //conversion TODO:remove
	    MSExperiment<> exp;
	    exp.push_back(current_cspectra_[x].getSpec());
	    topspec_->canvas()->addLayer(exp);
	    
	    exp.clear();
	    exp.push_back(current_cspectra_[y].getSpec());
	    bottomspec_->canvas()->addLayer(exp); 
      
      current_cspectra_[x].spec().updateRanges();
      current_cspectra_[y].spec().updateRanges();
      double lo_x =  min( current_cspectra_[x].getSpec().getMin()[0] , current_cspectra_[y].getSpec().getMin()[0]);
      double high_x = max( current_cspectra_[x].getSpec().getMax()[0] , current_cspectra_[y].getSpec().getMax()[0]) ;
      topspec_->canvas()->setVisibleArea(SpectrumCanvas::AreaType(lo_x,0,high_x,0));
      bottomspec_->canvas()->setVisibleArea(SpectrumCanvas::AreaType(lo_x,0,high_x,0));
    }
    else
    {
      topspec_->hide();
      bottomspec_->hide();
      topbin_->show();
      bottombin_->show();
      topbin_->loadBinnedRep(current_cspectra_[x].getBinrep());
      bottombin_->loadBinnedRep(current_cspectra_[y].getBinrep()); 
      double lo_x = min(current_cspectra_[x].getBinrep().min(),current_cspectra_[y].getBinrep().min());
      double high_x = max(current_cspectra_[x].getBinrep().max(),current_cspectra_[y].getBinrep().max());
      topbin_->setVisibleArea(lo_x,high_x);
      bottombin_->setVisibleArea(lo_x,high_x);
    }
    QApplication::restoreOverrideCursor();
  }

  void ClusterAnalysisWindow::fillspecListView_(int x, int y)
  {
    ostringstream ss;
    if ( x >= (int)current_cspectra_.size() || y >= (int)current_cspectra_.size() )
    {
      ss << "given indices don`t seem to fit into current_cspectra_\nx: " << x << 
        " y: " << y << " current_cspectra_.size() " << current_cspectra_.size();
      specinfoview_->setText(ss.str().c_str());
      return;
    }
    specinfoview_->clear();
    //get info from db
    uint id1 = current_cspectra_[x].id();
    uint id2 = current_cspectra_[y].id();
    double pz1 = current_cspectra_[x].getParentMass();
    double pz2 = current_cspectra_[y].getParentMass();
    uint charge1 = current_cspectra_[x].getParentionCharge();
    uint charge2 = current_cspectra_[y].getParentionCharge();
    double ret1 = current_cspectra_[x].getRetention();
    double ret2 = current_cspectra_[y].getRetention();
    double scor1 = current_cspectra_[x].getTophit().getScore();
    double scor2 = current_cspectra_[y].getTophit().getScore();
    QString seq1 = current_cspectra_[x].getTophit().getSequence().c_str();
    QString seq2 = current_cspectra_[y].getTophit().getSequence().c_str();
    //generate output
    ss.str("");;
    ss << "<qt>";
    ss << "<h4 > <p>" << QString("id") << " = " << QString("%1").arg(id1) << "</p> </h3>\n";
    ss << "<center>" << QString("parent m/z") << " = " << QString("%1").arg(pz1) << "</center>\n";
    ss << "<center>" << QString("precursor peak charge") << " = " << QString("%1").arg(charge1) << "</center>\n";
    ss << "<center>" << QString("retention") << "= " << QString("%1").arg(ret1) << "</center>\n";
    ss << "<center>" << QString("sequence") << " = " << QString(seq1) << "</center>\n";
    ss << "<center>" << QString("score") << " = " << QString("%1").arg(scor1) << "</center>\n";
    ss << "<h4 > <p>" << QString("id") << " = " << QString("%1").arg(id2) << "</p> </h3>\n";
    ss << "<center>" << QString("parent m/z") << " = " << QString("%1").arg(pz2) << "</center>\n";
    ss << "<center>" << QString("precursor peak charge") << " = " << QString("%1").arg(charge2) << "</center>\n";
    ss << "<center>" << QString("retention") << "= " << QString("%1").arg(ret2) << "</center>\n";
    ss << "<center>" << QString("sequence") << " = " << QString(seq2) << "</center>\n";
    ss << "<center>" << QString("score") << " = " << QString("%1").arg(scor2) << "</center>\n";
    ss << "</qt>";
    specinfoview_->setText(ss.str().c_str());
  }

  /**
  \param clustering to show
   */
  void ClusterAnalysisWindow::showClustering(const map<int,ClusterNode*>& clustering)
  {
    clusterlistview_->clear();
    clusterp_ = &clustering;
    for (std::map<int,ClusterNode*>::const_iterator it = clustering.begin(); it != clustering.end(); ++it)
    {
      new QListViewItem(clusterlistview_,QString("%1").arg(it->first),QString("size %1").arg(it->second->size(),3));
    }
  }

  void ClusterAnalysisWindow::choosefunctors()
  {
    ClusterRunWidget* pwi = new ClusterRunWidget(this);
    pwi->exec();
    pwi->raise();
    pwi->setActiveWindow();
    ClusterExperiment::ClusterRun* crp =  pwi->getClusterRun();
    clusterrun_ = *crp;
    delete crp;
    delete pwi;
  }

  void ClusterAnalysisWindow::inspect(int /*id1*/ , int /*id2*/, bool /*preprocess1*/, bool /*preprocess2*/)
  {
//TODO Persistence
//    // get input
//    if (id1 == 0)
//    {
//      InspectDialog* temp = new InspectDialog(this);
//      temp->exec();
//      temp->raise();
//      temp->setActiveWindow();
//      if ( temp->result() != QDialog::Accepted ) return;
//      
//      // if multiple ids are inserted in id1, a cluster is created
//      QStringList idsl = QStringList::split(" ",temp->text1());
//      if (idsl.size() > (uint)1) 
//      { 
//        stringstream ss;
//        vector<int> peaklists;
//        map<string,int> unknown;
//        for (uint i = 0; i < idsl.size(); ++i)
//        {
//          uint id = idsl[i].toInt();
//          string type = adapter_->type(id);
//          // check if id is a PeakList
//          if ( type == "PeakList" )
//          {
//            peaklists.push_back(id);
//          }
//          else
//          {
//            unknown.insert(make_pair(type,id));
//          }
//        }
//        if ( unknown.size() )
//        {
//          ss.str("");
//          ss << "cant cluster these ids:\n";
//          for ( map<string,int>::const_iterator cmit = unknown.begin(); cmit != unknown.end(); ++cmit )
//          {
//            ss << "\t" << cmit->first << " : " << cmit->second << "\n";
//          }
//          QMessageBox::warning(this,"OpenMS-CAt",ss.str().c_str());
//        }
//        if ( peaklists.size() )
//        {
//          showCluster(peaklists);
//        }
//        return;
//      }
//      
//      id1 = temp->selection1();
//      id2 = temp->selection2();
//     
//      preprocess1 = temp->preprocess1();
//      preprocess2 = temp->preprocess2();
//      
//      delete temp;
//    }
//
//    if ( id1 == 0 )
//    {
//      return;
//    }
//    
//    if ( adapter_->type(id1) != "PeakList" || ( id2 > 0 && adapter_->type(id2) != "PeakList" ))
//    {
//      uint wrong = 0;
//      stringstream ss;
//      if ( adapter_->type(id1) != "PeakList" )
//      {
//        ++wrong;
//        ss << id1;
//      }
//      if ( id2 > 0 && adapter_->type(id2) != "PeakList" )
//      {
//        if ( wrong ) 
//        {
//          ss << " and ";
//        }
//        ss << id2;
//      }
//      if ( wrong == 2 ) 
//      {
//        ss << " are not PeakLists";
//      }
//      else
//      {
//        ss << " is not a PeakList";
//      }
//      QMessageBox::warning(this,"OpenMS-CAt",ss.str().c_str());
//      return;
//    }
//   
//    current_cspectra_.clear();
//    
//    addClusterSpectrum(ClusterSpectrum(id1,adapter_,clusterrun_.getBinSize(),clusterrun_.getBinSpread()),preprocess1);
//    if ( id2 > 0 )
//    {
//      addClusterSpectrum(ClusterSpectrum(id2,adapter_,clusterrun_.getBinSize(),clusterrun_.getBinSpread()),preprocess2);
//    }
//    // theoretical
//    else if ( id2 == -2 )
//    {
//      string sequence = "";
//      sequence = current_cspectra_.begin()->getTophit().getSequence();
//      addClusterSpectrum(ClusterSpectrum(SpectrumGenerator::instance()->getspectrum(sequence),0,clusterrun_.getBinSize(),clusterrun_.getBinSpread()),preprocess2);
//    }
//    // self
//    else
//    {
//      addClusterSpectrum(ClusterSpectrum(id1,adapter_,clusterrun_.getBinSize(),clusterrun_.getBinSpread()),preprocess2);
//    }
  }

  void ClusterAnalysisWindow::clusterListViewDoubleClick(QListViewItem* clitem)
  {
    if (clitem->listView() == clusterlistview_)
    {
      if (clitem->depth() == 0) // cluster
      {
        ClusterNode* cluster = (*const_cast<map<int,ClusterNode*>* >(clusterp_))[clitem->text(0).toInt()];
        vector<int> clusterids;
        for ( list<SignedInt>::const_iterator cit = cluster->children().begin(); cit != cluster->children().end(); ++cit )
        {
          clusterids.push_back(*cit);
        }
        showCluster(clusterids,clitem);
      }
      else //peptide
      {
        int id = clitem->text(0).toInt();
        for (uint i = 0; i < ids_.size(); ++i)
        {
          if ( ids_[i] == id ) 
          {
            clmatrix_->highlight(i);
            break;
          }
        }
        uint index = 0;
        for ( uint i = 0; i < current_cspectra_.size() ; ++i)
        {
          if (current_cspectra_[i].id() == id )
          {
            index = i;
          }
        }
        displaySpectra(index,index);
      }
    }
  }

  /**
  \param cspec new ClusterSpectrum
  \param preprocess whether to use preprocessing
  */
  void ClusterAnalysisWindow::addClusterSpectrum(const ClusterSpectrum& cspec,bool preprocess)
  {
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    current_cspectra_.push_back(cspec);
    if ( preprocess ) clusterrun_.preprocess(current_cspectra_.rbegin()->spec());
    ids_.push_back(current_cspectra_.rbegin()->id());
    createMatrix_();
    clmatrix_->setMatrix(current_matrix_);
    QApplication::restoreOverrideCursor();
  }

  /**
  spectra are cached in current_cspectra_
  \param ids of spectra in clster
  \param clitem QListViewItem to expand
   */
  void ClusterAnalysisWindow::showCluster(const vector<int>& ids, QListViewItem* clitem)
  {
    current_cspectra_.clear();
    ids_.clear();
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    for ( uint i = 0; i < ids.size(); ++i )
    {
      current_cspectra_.push_back(ClusterSpectrum(ids[i],adapter_,clusterrun_.getBinSize(),clusterrun_.getBinSpread()));
      try
      {
        clusterrun_.preprocess(current_cspectra_.rbegin()->spec());
        ids_.push_back(current_cspectra_.rbegin()->id());
        if ( clitem != 0 && clitem->childCount() < (int)ids.size() )
        {
          QString id = QString("%1").arg(current_cspectra_.rbegin()->id());
          QString sequence = current_cspectra_.rbegin()->getTophit().getSequence().c_str();
          new QListViewItem(clitem,id,sequence);
        }
      }
      catch (Exception::Base& e)
      {
        cerr << "id " << ids[i] << " is no PeakList!\n";
        current_cspectra_.pop_back();
      }
    }
    if ( current_cspectra_.size() )
    {
      createMatrix_();
      clmatrix_->setMatrix(current_matrix_);
    }
    QApplication::restoreOverrideCursor();
  }

  void ClusterAnalysisWindow::connect2DB()
  {
//TODO Persistence
//    if ( !dbdialog_ ) dbdialog_ = new DBDialog(this);
//    dbdialog_->exec();
//    dbdialog_->raise();
//    dbdialog_->setActiveWindow();
//    if (dbdialog_->result() == QDialog::Accepted )
//    {
//      adapter_ = dbdialog_->adapter();
//      menubar_->setItemEnabled(2,1);
//      menubar_->setItemEnabled(3,1);
//      menubar_->setItemEnabled(4,1);
//      menubar_->setItemEnabled(5,1);
//      clev_->setDBAdapter(adapter_);
//    }
  }

  void ClusterAnalysisWindow::loadClusterExperiment()
  {
    QFileDialog* filedialog = new QFileDialog(this);
    //debug
    filedialog->exec();
    filedialog->raise();
    filedialog->setActiveWindow();
    if (filedialog->result() == QDialog::Accepted)
    {
      QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
      string filename = filedialog->selectedFile().ascii();
      clev_->load(filename);
      QApplication::restoreOverrideCursor();
    }
  }

}
