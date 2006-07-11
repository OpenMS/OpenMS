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
// $Maintainer:  $
// --------------------------------------------------------------------------

#include "AnnotateThread.h"
#include "Annotate.h"
#include "CustomEvents.h"
#include "../config_specannotate.h"

#include <qapplication.h>


using namespace OpenMS;


AnnotateThread::AnnotateThread(__gnu_cxx::hash_map<std::string, std::string> sample_data,
                               std::vector<OpenMS::DPeakArray<1, OpenMS::DPeak<1> >::iterator>& peaklist,
                               std::vector<std::string> ov_mods,
                               std::string db_username, std::string db_password, std::string db_host,
                               OpenMS::Annotate* qannotate) : QThread()
{
  sample_data_ = sample_data;
  peaklist_ = peaklist;
  ov_mods_ = ov_mods;
  db_username_ = db_username;
  db_password_ = db_password;
  db_host_ = db_host;
  qannotate_ = qannotate;

}



void AnnotateThread::run()
{
  Sample sample(sample_data_, peaklist_, ov_mods_, db_username_, db_password_, db_host_, qannotate_);
  sample.annotate();
  if (sample_data_["outputdir"] != "")
    {
      sample.printAnnotations();
    }
  else
    {

      sample.storeAnnotations();


      //std::pair<OpenMS::DPeakArray<1, OpenMS::DPeak<1> >*, std::vector<std::vector<Annotation> > > annot_peaks = sample.getAnnotations();
    }

  //clean up database connections
  // while(QSqlDatabase::contains(QSqlDatabase::defaultConnection))
  //     {
  //       QSqlDatabase::removeDatabase(QSqlDatabase::defaultConnection);
  //     }
#ifndef ANNOTATE_XML
  while(QSqlDatabase::contains("db_handle_"))
    {
      QSqlDatabase::removeDatabase("db_handle_");
    }
#endif

  //send the finish signal!!
  FinishEvent* fe = new FinishEvent();
  QApplication::postEvent(qannotate_, fe);
}


