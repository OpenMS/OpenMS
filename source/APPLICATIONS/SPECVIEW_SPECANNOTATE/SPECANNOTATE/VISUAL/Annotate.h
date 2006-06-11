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
//
// --------------------------------------------------------------------------
// $Id: Annotate.h,v 1.3 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef ANNOTATE_H
#define ANNOTATE_H

#include <qvariant.h>
#include <qdialog.h>
#include <qstring.h>
#include <qmap.h>
#include <string>
#include <vector>
#include <ext/hash_map>
#include <qdatetime.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include "../FUNCTION/string_hash_stl_fixes.h"
#include "AnnotateThread.h"
#include <qsqldatabase.h>


class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QLabel;
class QTextBrowser;
class QLCDNumber;
class QPushButton;


namespace OpenMS
  {

  class Annotate : public QDialog
    {
      Q_OBJECT

    public:
      Annotate( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
      ~Annotate();

      QLabel* textLabel2;
      QTextBrowser* textBrowser1;
      QLabel* textLabel1;
      QLCDNumber* lCDNumber1;
      QLabel* textLabel1_2;
      QLCDNumber* lCDNumber2;
      QLabel* textLabel2_2;
      QLCDNumber* lCDNumber3;
      QLabel* textLabel3;
      QLCDNumber* lCDNumber4;
      QLabel* textLabel1_3;
      QLCDNumber* lCDNumber5;
      QLabel* textLabel2_3;
      QLCDNumber* lCDNumber6;
      QPushButton* pushButton1;
      QPushButton* pushButton3;

    public slots:
      virtual void run( __gnu_cxx::hash_map<std::string, std::string> sample_data, std::vector<OpenMS::Spectrum1DWidget::Spectrum1D::iterator> & peaklist, std::vector<std::string> ov_mods, QMap<QString, QString> * settings );
      virtual void addOutput( std::string s );
      virtual void ready();
      virtual void abort();
      virtual void closeWindow();
      virtual void timerEvent( QTimerEvent * e );
      virtual void customEvent( QCustomEvent * e );

    protected:
      QVBoxLayout* AnnotateLayout;
      QVBoxLayout* layout3;
      QHBoxLayout* layout5;
      QHBoxLayout* layout4;
      QHBoxLayout* layout5_2;
      QHBoxLayout* layout6;
      QHBoxLayout* layout8;
      QHBoxLayout* layout9;
      QHBoxLayout* layout9_2;

    protected slots:
      virtual void languageChange();

      virtual void updateDBDisplay();
      void dbConnect();

    private:
      int db_display_update_timer;
      QTime t;
      int timer_ID;
      QMap<QString, QString>* settings_;
      AnnotateThread* qathread;
      QSqlDatabase* defaultDB;

    };
}
#endif // ANNOTATE_H
