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


#ifndef __ANNOTATETHREAD_H__
#define __ANNOTATETHREAD_H__


//STL includes
#include <string>

//QT includes
#include <qthread.h>
#include <qevent.h>

//OpenMS includes
#include <OpenMS/VISUAL/Spectrum1DWidget.h>


#include "../FUNCTION/Sample.h"


namespace OpenMS
  {

  class Annotate;

  //Thread for running the part of SpecAnnotate containing the actual functionality
  class AnnotateThread : public QThread
    {
    private:
      __gnu_cxx::hash_map<std::string, std::string> sample_data_;
      std::vector<std::string> ov_mods_;
      std::vector<Spectrum1DWidget::Spectrum1D::iterator> peaklist_;
      std::string db_username_;
      std::string db_password_;
      std::string db_host_;
      Annotate* qannotate_;

      friend class Annotate;

    public:

      AnnotateThread() : QThread()
      {}
      ;
      AnnotateThread(__gnu_cxx::hash_map<std::string, std::string> sample_data,
                     std::vector<OpenMS::Spectrum1DWidget::Spectrum1D::iterator>& peaklist,
                     std::vector<std::string> ov_mods,
                     std::string db_username, std::string db_password, std::string db_host, Annotate* qannotate);

      void run();


    };
}
#endif
