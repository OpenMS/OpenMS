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
#include <OpenMS/VISUAL/DIALOGS/SortDialog.h>

using namespace std;

namespace OpenMS
{

  SortDialog::SortDialog(vector<const ClusterExperiment::Analysis*> anafuncs, QWidget* parent, const char* name)
    :QDialog(parent,name),anafuncs_(anafuncs)
  {
    analysislabel_ = new QLabel("Analysis",this);
    measurelabel_ = new QLabel("Measure",this);

    analysisbox_ = new QComboBox(this);
    measurebox_ = new QComboBox(this);

    ok_ = new QPushButton("ok",this);

    gridlayout_ = new QGridLayout(this,2,3);
    gridlayout_->addWidget(analysislabel_,0,0);
    gridlayout_->addWidget(measurelabel_,1,0);
    gridlayout_->addWidget(analysisbox_,0,1);
    gridlayout_->addWidget(measurebox_,1,1);
    gridlayout_->addWidget(ok_,2,2);

    connect(ok_,SIGNAL(clicked()),this,SLOT(ok()));
    connect(analysisbox_,SIGNAL(textChanged(const QString &)),this,SLOT(fillmeasure()));

    fillanalysis_();
    fillmeasure();
  }

  void SortDialog::fillanalysis_()
  {
    for ( vector<const ClusterExperiment::Analysis*>::const_iterator cvit = anafuncs_.begin(); cvit != anafuncs_.end(); ++cvit )
    {
      analysisbox_->insertItem((*cvit)->name().c_str());
    }
  }

  void SortDialog::fillmeasure()
  {
    for ( vector<const ClusterExperiment::Analysis*>::const_iterator cvit = anafuncs_.begin(); cvit != anafuncs_.end(); ++cvit )
    {
      if ( analysisbox_->currentText().ascii() == (*cvit)->name() )
      {
        measurebox_->clear();
        for ( map<String,double>::const_iterator cmit = (*cvit)->results().begin(); cmit != (*cvit)->results().end(); ++cmit )
        {
          measurebox_->insertItem(cmit->first.c_str());
        }
      }
    }
  }

  void SortDialog::ok()
  {
    accept();
  }

  QString SortDialog::analysis()
  {
    return analysisbox_->currentText();
  }

  QString SortDialog::measure()
  {
    return measurebox_->currentText();
  }

}
