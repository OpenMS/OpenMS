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
#ifndef OPENMS_VISUAL_DIALOGS_FACTORYPRODUCTDIALOG_H
#define OPENMS_VISUAL_DIALOGS_FACTORYPRODUCTDIALOG_H

#include <qdialog.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qpushbutton.h>

#include <map>
#include <string>

#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{

  /**
   allows to change the paramaters for a FactoryProduct <br>
   */
  class FactoryProductDialog : public QDialog
  {
    Q_OBJECT
  public:
    FactoryProductDialog(QWidget* = 0, const char* = 0);
    void setFactoryProduct(FactoryProduct*);
    uint size() const { return paramwidgets_.size();}
  public slots:
    void ok();
  private:
    FactoryProduct* configurable_;
    QGridLayout* gridlayout_;
    std::map<std::string,std::pair<QLabel*,QLineEdit*> > paramwidgets_;
    QPushButton* ok_;
  };
}
#endif //OPENMS_VISUAL_DIALOGS_FACTORYPRODUCTDIALOG
