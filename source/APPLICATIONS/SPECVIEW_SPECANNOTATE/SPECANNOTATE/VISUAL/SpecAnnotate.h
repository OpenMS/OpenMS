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
// $Id: SpecAnnotate.h,v 1.3 2006/03/28 08:03:28 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------



#ifndef SPECANNOTATE_H
#define SPECANNOTATE_H

//QT includes
#include <qmainwindow.h>

//OpenMS includes
#include <OpenMS/FORMAT/Param.h>


//forward declaration
class QString;
template<class Key, class T>
class QMap;



namespace OpenMS
  {

  class SpecAnnotate: public QMainWindow
    {
      Q_OBJECT

    public:
      SpecAnnotate();
      ~SpecAnnotate();

      QMap<QString,QString>* getSettings();

    protected:
      bool connectToDatabase();

    private:
      Param main_param_;
      QMap<QString,QString>* settings_;

    public slots:
      void loadSettings();
      void quit();

    private slots:

      void newSample();
      void loadSample();
      void saveSample();
      void openSettingsDialog();

      void dbMod();
      void dbEnz();
      void dbAmi();
      void dbProt();
      void dbSeq();
      void dbDigFrag();
      void dbModComb();
      void dbModCombPosless();
      void dbProtModScen();
      void dbRealMod();
      void dbRealModPosless();
      void dbSample();
      void dbAnnot();

      void xmlEnzyme();
      void xmlModification();

      void annotate();
      void resetDB();
      void setupDB();
      void insertProtIntoDB();
      void insertProtIntoXML();

      void about();
      void aboutQt();
      void usageWithoutDB();
      void threeMethods();

    };
}

#endif
