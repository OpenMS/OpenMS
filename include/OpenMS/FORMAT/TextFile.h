// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_TEXTFILE_H
#define OPENMS_FORMAT_TEXTFILE_H

#include <OpenMS/DATASTRUCTURES/StringList.h>

namespace OpenMS
{
	/**
		@brief This class provides some basic file handling methods for text files.
  
  	@ingroup FileIO
	*/
  class OPENMS_DLLAPI TextFile
  	: public StringList
  {
  	
    public:

    	///Default constructor
			TextFile();

			/// destructor
			virtual ~TextFile();
			
    	/**
    		@brief Constructor with filename
    	
    		@param filename The input file name.
    		@param trim_lines Whether or not the lines are trimmed when reading them from file.
    		@param first_n If set, only @p first_n lines the lines from the beginning of the file are read.

				@exception Exception::FileNotFound is thrown if the file could not be opened.
    	*/
			TextFile(const String& filename, bool trim_lines=false, Int first_n=-1);

    	/**
    		@brief Loads data from a text file.
    	
    		@param filename The input file name.
    		@param trim_lines Whether or not the lines are trimmed when reading them from file.
    		@param first_n If set, only @p first_n lines the lines from the beginning of the file are read.

				@exception Exception::FileNotFound is thrown if the file could not be opened.
    	*/
			void load(const String& filename, bool trim_lines=false, Int first_n=-1);

    	/**
    		@brief Writes the data to a file
    		
    		@note this function uses unix-style linebreaks

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
    	*/
			void store(const String& filename);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_TEXTFILE_H
