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
// $Id: FactoryProductDialog.C,v 1.2 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#include <OpenMS/VISUAL/DIALOGS/FactoryProductDialog.h>

using namespace std;

namespace OpenMS
{

  FactoryProductDialog::FactoryProductDialog(QWidget* parent, const char* name)
    :QDialog(parent,name),paramwidgets_()
  {
    ok_ = new QPushButton("ok",this);
    connect(ok_,SIGNAL(released()),this,SLOT(ok()));
  }

  void FactoryProductDialog::setFactoryProduct(FactoryProduct* configurable)
  {
    configurable_ = configurable;
    Param mfparams = configurable_->getParam();
    gridlayout_ = new QGridLayout(this,2,mfparams.size() +1 );
    int i = 0;
    for (Param::ConstIterator it = mfparams.begin(); it != mfparams.end(); ++it)
    {
      paramwidgets_.insert(make_pair(it->first,make_pair(new QLabel(it->first.c_str(),this),new QLineEdit(QString("%1").arg((String)it->second),this))));
      gridlayout_->addWidget(paramwidgets_[it->first].first,i,0);
      gridlayout_->addWidget(paramwidgets_[it->first].second,i,1);
      ++i;
    }
    gridlayout_->addWidget(ok_,i,1);
  }

  void FactoryProductDialog::ok()
  {
    for (map<string,pair<QLabel*,QLineEdit*> > ::iterator it = paramwidgets_.begin(); it != paramwidgets_.end();++it)
    {
      configurable_->getParam().setValue(it->first,it->second.second->text().toDouble()); 
    }
    accept();
  }

}
