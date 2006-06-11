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
// $Id: SequestOutfile.h,v 1.4 2006/06/10 06:40:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FORMAT_SEQUESTOUTFILE_H
#define OPENMS_FORMAT_SEQUESTOUTFILE_H


#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <fstream>
#include <iostream>
#include <string>
//#include <OpenMS/KERNEL/PeptideHit.h>

namespace OpenMS
{
  class SequestOutfile
  {
    public:
      SequestOutfile(const std::string& filename = "");
      SequestOutfile(const SequestOutfile& outfile);
      ~SequestOutfile();
      SequestOutfile& operator>>(Identification& identification);
      SequestOutfile& operator>>(PeptideHit& peptidehit);
    protected:
      static int length_(const char* charp);
      static int containlength_(const char* charp);
      void replaceapostrophe_(std::string& str);
      std::string filename_;
      bool phfinished_;
      int phstartingpos_;
  };
} //namespace OpenMS

#endif // OPENMS_FORMAT_SEQUESTOUTFILE_H

