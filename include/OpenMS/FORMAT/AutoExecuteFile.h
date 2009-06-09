// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    AutoExecute file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz
// $Authors: Guillaume Belz
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_AUTOEXECUTEFILE_H
#define OPENMS_FORMAT_AUTOEXECUTEFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <fstream>
using namespace std;

namespace OpenMS
{
 	/**
 		@brief File adapter for autoExecute file or list file.
  	
    For example, to merge fid files "find source | grep fid > destination"
  
  	@ingroup FileIO
  */
  
  class OPENMS_DLLAPI AutoExecuteFile
  {
// 	  private:
//		  PeakFileOptions options_;

    public:
      /// Default constructor
      AutoExecuteFile();
      
      /// Import file list
      StringList getFileList(
	      const String & filename,  
	      const bool isAutoExecute = false,
	      const unsigned int begin = 0,
	      const unsigned int end = 0,
	      const String src_dir = "");
	      
	  private:
	    unsigned int PosOnScout_; // position on chip, format [A-P]:[1-24]
	    unsigned int SpectrumDirectory_; // spectrum directory
	    unsigned int SpectrumFilename_; // spectrum filename
	    unsigned int ChipOnScout_; // normal or calibrant position, "0" or "1"

      String autoExecuteToFilename(const String & line, const String & src_dir);
      void readAutoExecuteHeader(std::ifstream & is_);
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_AUTOEXECUTEFILE_H

