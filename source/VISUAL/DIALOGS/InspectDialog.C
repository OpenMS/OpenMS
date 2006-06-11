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
// $Id: InspectDialog.C,v 1.2 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#include <OpenMS/VISUAL/DIALOGS/InspectDialog.h>

namespace OpenMS
{

  InspectDialog::InspectDialog(QWidget* parent, const char* name)
    : QDialog(parent,name)
  {
    label1_ = new QLabel("id1",this);
    label2_ = new QLabel("id2",this);

    edit1_ = new QLineEdit(this);
    box2_ = new QComboBox(this);
    box2_->setEditable(1);
    box2_->insertItem("enter id");
    box2_->insertItem("self");
    box2_->insertItem("theoretical");

    pp1_ = new QCheckBox("preprocess",this);
    pp2_ = new QCheckBox("preprocess",this);
    
    ok_ = new QPushButton("ok",this);

    gridlayout_ = new QGridLayout(this,2,3);

    gridlayout_->addWidget(label1_,0,0);
    gridlayout_->addWidget(label2_,1,0);
    gridlayout_->addWidget(edit1_,0,1);
    gridlayout_->addWidget(box2_,1,1);
    gridlayout_->addWidget(pp1_,0,2);
    gridlayout_->addWidget(pp2_,1,2);

    gridlayout_->addWidget(ok_,2,1);

    connect(ok_,SIGNAL(clicked()),this,SLOT(ok()));
  }

  void InspectDialog::ok()
  {
    accept();
  }

  QString InspectDialog::text1()
  {
    return edit1_->text();
  }

  int InspectDialog::selection1()
  {
    return edit1_->text().toInt();
  }

  // result = id of spectrum or 0 if nothing was given, 
  //   -1 for unpreprocessed
  //   -2 for theoretical
  int InspectDialog::selection2()
  {
    int result;
    if ( box2_->currentText() == "enter id" )
    {
      result = 0;
    }
    else if ( box2_->currentText() == "self" )
    {
      result = -1;
    }
    else if ( box2_->currentText() == "theoretical" )
    {
      result = -2;
    }
    else 
    {
      result = box2_->currentText().toInt();
    }
    return result;
  }

  bool InspectDialog::preprocess1()
  {
    return pp1_->isChecked();
  }

  bool InspectDialog::preprocess2()
  {
    return pp2_->isChecked();
  }
}
