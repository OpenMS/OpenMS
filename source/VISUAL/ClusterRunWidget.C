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
#include <OpenMS/VISUAL/ClusterRunWidget.h>

#include <OpenMS/VISUAL/DIALOGS/FactoryProductDialog.h>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>

#include <qinputdialog.h>

using namespace std;

namespace OpenMS
{

  ClusterRunWidget::ClusterRunWidget( QWidget* parent, const char* name )
    : QDialog(parent,name),cfp_(0),mowers_(),binsize_(1),binspread_(0)
  {
    gridlayout_ = new QGridLayout(this,5,3);

    cflabel_ = new QLabel("CompareFunctor",this);
    pplabel_ = new QLabel("Mowers",this);
    
    cfbox_ = new QComboBox(this);
    ppbox_ = new QComboBox(this);
    ppnames_ = new QLabel(this);

    addpp_ = new QPushButton("add",this);
    clearpp_ = new QPushButton("clear",this);
    usecf_ = new QPushButton("use",this);
    
    ok_ = new QPushButton("ok",this);
    
    gridlayout_->addWidget(cflabel_,0,0);
    gridlayout_->addWidget(cfbox_,0,1);
    gridlayout_->addWidget(usecf_,0,2);
    gridlayout_->addWidget(pplabel_,1,0);
    gridlayout_->addWidget(ppbox_,1,1);
    gridlayout_->addWidget(addpp_,1,2);
    gridlayout_->addWidget(ppnames_,2,1);
    gridlayout_->addWidget(clearpp_,2,0);
    gridlayout_->addWidget(ok_,4,2);
    
    fillbox_("CompareFunctor",cfbox_);
    fillbox_("PreprocessingFunctor",ppbox_);
    connect(addpp_,SIGNAL(clicked()),this,SLOT(addpp()));
    connect(clearpp_,SIGNAL(clicked()),this,SLOT(clearpp()));
    connect(usecf_,SIGNAL(clicked()),this,SLOT(usecf()));
    connect(ok_,SIGNAL(clicked()),this,SLOT(ok()));
  }

  void ClusterRunWidget::configure_(FactoryProduct* cf)
  {
    if ( cf->getParam().size() )
    {
      FactoryProductDialog* cfd = new FactoryProductDialog();
      cfd->setFactoryProduct(cf);
      cfd->exec();
      cfd->raise();
      cfd->setActiveWindow();
      delete cfd;
    }
  }

  void ClusterRunWidget::usecf()
  {
    delete cfp_;
    cfp_ = dynamic_cast<CompareFunctor*>(ClusterFactory::instance()->create(cfbox_->currentText().ascii()));
    configure_(cfp_);
  }

  void ClusterRunWidget::addpp()
  {
    PreprocessingFunctor* mfp = dynamic_cast<PreprocessingFunctor*>(ClusterFactory::instance()->create(ppbox_->currentText().ascii()));
    configure_(mfp);
    mowers_.push_back(mfp);
    ppnames_->setText(ppnames_->text() + "\n" + ppbox_->currentText());
  }

  void ClusterRunWidget::clearpp()
  {
    for ( uint i = 0; i < mowers_.size(); ++i )
    {
      delete mowers_[i];
    }
    mowers_.clear();
    ppnames_->clear();
  }

  void ClusterRunWidget::fillbox_(String type, QComboBox* box)
  {
    vector<String> candidates = ClusterFactory::instance()->catalogue(type);
    box->clear();
    for ( uint i = 0; i < candidates.size(); ++i )
    {
      box->insertItem(candidates[i].c_str());
    }
  }

  void ClusterRunWidget::ok()
  {
    if ( !cfp_ )
    {
      usecf();
    }
    if ( cfp_->usebins() )
    {
      binsize_ = QInputDialog::getDouble("OpenMS-CAt","please insert binsize",1);
      binspread_ = QInputDialog::getInteger("OpenMS-CAt","please insert binspread");
    }
    accept();
  }

  ClusterExperiment::ClusterRun* ClusterRunWidget::getClusterRun()
  {
    ClusterExperiment::ClusterRun* result = new ClusterExperiment::ClusterRun();
		/*
		TODO
    if (cfp_) 
    {
      CompareFunctor* cfp = dynamic_cast<CompareFunctor*>(ClusterFactory::instance()->duplicate(cfp_));
      result->setSimFunc(cfp);
    }
    for ( vector<PreprocessingFunctor*>::const_iterator cvit = mowers_.begin();
        cvit != mowers_.end(); ++cvit )
    {
      PreprocessingFunctor* mfp = dynamic_cast<PreprocessingFunctor*>(ClusterFactory::instance()->duplicate(*cvit));
      result->addMower(mfp);
    }
    result->setBinSize(binsize_);
    result->setBinSpread(binspread_);
		*/
    return result;
  }

}
