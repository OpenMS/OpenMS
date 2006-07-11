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
#ifndef OPENMS_VISUAL_DIALOGS_INSPECTDIALOG_H
#define OPENMS_VISUAL_DIALOGS_INSPECTDIALOG_H

#include <qdialog.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qcombobox.h>
#include <qcheckbox.h>
#include <qpushbutton.h>
#include <qlayout.h>

#include <string>

namespace OpenMS
{

  class InspectDialog : public QDialog
  {
    Q_OBJECT
  public:
    InspectDialog(QWidget* = 0, const char* = 0);
    
    QString text1();
    
    int selection1();
    int selection2();
    
    bool preprocess1();
    bool preprocess2();
  public slots:
    void ok();
  private:
    QGridLayout* gridlayout_;
    QLabel* label1_;
    QLabel* label2_;
    QLineEdit* edit1_;
    QComboBox* box2_;
    QCheckBox* pp1_;
    QCheckBox* pp2_;
    QPushButton* ok_;
  };
}
#endif //OPENMS_VISUAL_DIALOGS_INSPECTDIALOG_H
