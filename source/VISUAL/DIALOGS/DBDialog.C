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
// $Id: DBDialog.C,v 1.4 2006/03/28 08:03:39 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#include <OpenMS/VISUAL/DIALOGS/DBDialog.h>

#include <OpenMS/FORMAT/DBConnection.h>

#include <qmessagebox.h>

#include <fstream>

using namespace std;

namespace OpenMS
{

  void DBDialog::createWidgets_()
  {
    hostlabel_ = new QLabel(tr("Host"),this);
    dblabel_ = new QLabel(tr("DB"),this);
    loginlabel_ = new QLabel(tr("login"),this);
    pwlabel_ = new QLabel(tr("password"),this);

    hostedit_ = new QLineEdit(this);
    dbedit_ = new QLineEdit(this);
    loginedit_ = new QLineEdit(this);
    pwedit_ = new QLineEdit(this);
    pwedit_->setEchoMode(QLineEdit::Password);

    cblogin_ = new QCheckBox("save logininfo",this);
    cbpw_ = new QCheckBox("save password",this);
    cbpw_->setEnabled(0);
    
    okbutton_ = new QPushButton(tr("Ok"),this);
  }

  void DBDialog::doLayout_()
  {
    layout_ = new QGridLayout(this,8,4,10,5);
    layout_->addWidget(hostlabel_,0,0);
    layout_->addWidget(hostedit_,0,1);
    layout_->addWidget(dblabel_,0,2);
    layout_->addWidget(dbedit_,0,3);
    layout_->addWidget(loginlabel_,1,0);
    layout_->addWidget(loginedit_,1,1);
    layout_->addWidget(pwlabel_,1,2);
    layout_->addWidget(pwedit_,1,3);
    layout_->addMultiCellWidget(cblogin_,2,2,0,1);
    layout_->addMultiCellWidget(cbpw_,3,3,0,1);
    layout_->addWidget(okbutton_,3,3);
  }

  DBDialog::DBDialog(QWidget* parent, const char* name, int modal)
    :QDialog(parent,name,modal),adapter_(0)
  {
    createWidgets_();
    doLayout_();
    connect_();
    getContents();
  }

  void DBDialog::saveContents()
  {
    ofstream pref("preferences");
    pref << hostedit_->text() << " " << dbedit_->text() << " " << loginedit_->text();
    if ( cbpw_->isChecked() ) pref << " " << pwedit_->text();
  }

  void DBDialog::getContents()
  {
    ifstream pref("preferences");
    string host,db,login,pw;
    pref >> host ;
    pref >> db ;
    pref >> login ;
    pref >> pw;
    hostedit_->setText(host.c_str());
    dbedit_->setText(db.c_str());
    loginedit_->setText(login.c_str());
    pwedit_->setText(pw.c_str());
  }

  void DBDialog::connect_()
  {
    connect(okbutton_,SIGNAL(released()),this,SLOT(ok()));
    connect(cblogin_,SIGNAL(released()),this,SLOT(savepressed()));
  }

  void DBDialog::savepressed()
  {
    if ( cblogin_->isChecked() )
    {
      QMessageBox::warning( this, "OpenMS-CAt","All information\nincluding the password is saved unencrypted!\n");
      cbpw_->setEnabled(1);
    }
    else
    {
      cbpw_->setEnabled(0);
    }
  }

  void DBDialog::ok()
  {
    //todo: validate qlineedits
    if ( ! adapter_ ) 
    adapter_ = new DBConnection();
    try
    {
      adapter_->connect(dbedit_->text().ascii(),loginedit_->text().ascii(),pwedit_->text().ascii(),hostedit_->text().ascii());
    }
    catch (DBConnection::InvalidQuery e)
    {
      setCaption(tr("could not connect to host"));
      std::cerr << e.what() << std::endl;
      return;
    }
    if ( cblogin_->isChecked() )
    {
      saveContents();
    }
    accept();
  }

}

