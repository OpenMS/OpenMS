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
// $Id: DBDialog.h,v 1.4 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------
//

#ifndef OPENMS_VISUAL_DIALOGS_DBDIALOG_H
#define OPENMS_VISUAL_DIALOGS_DBDIALOG_H

#include <qdialog.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qpushbutton.h>
#include <qlayout.h>
#include <qcheckbox.h>

#include <OpenMS/FORMAT/DBConnection.h>

namespace OpenMS
{
  // lets the user connect to a OpenMS DB
  // the settings can be saved to HDD (filename preferences)
  // but the settings are saved in plaintext
  class DBDialog : public QDialog
  {
    Q_OBJECT
  public:

    DBDialog(QWidget* = 0, const char* = 0 , int = 1);

    //accessors
    DBConnection* adapter() {return adapter_;}

  public slots:

    void ok();
    void getContents();
    void saveContents();
    void savepressed();

  protected:

  private:

    void createWidgets_();
    void doLayout_();
    void connect_();

    DBConnection* adapter_;

    QLabel* hostlabel_;
    QLabel* dblabel_;
    QLabel* loginlabel_;
    QLabel* pwlabel_;

    QLineEdit* hostedit_;
    QLineEdit* dbedit_;
    QLineEdit* loginedit_;
    QLineEdit* pwedit_;

    QPushButton* okbutton_;
    QCheckBox* cblogin_;
    QCheckBox* cbpw_;

    QGridLayout* layout_;

  };

}
#endif // OPENMS_VISUAL_DIALOGS_DBDIALOG_H

