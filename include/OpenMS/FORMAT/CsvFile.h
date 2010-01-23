// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CSVFILE_H
#define OPENMS_FORMAT_CSVFILE_H

#include <OpenMS/FORMAT/TextFile.h>

namespace OpenMS
{
	class StringList;
	/**
		@brief This class handles csv files. Currently only loading is implemented.

		@note items are allowed to be enclosed by only one character e.g. "item" where " is enclosing character
	
  	@ingroup FileIO
	*/
  class OPENMS_DLLAPI CsvFile
  	: public TextFile
  {
    public:

    	///Default constructor
			CsvFile();

			/// destructor
			virtual ~CsvFile();
			
    	/**
    		@brief Constructor with filename
    	
    		@param  filename The input file name.
    		@param  is character which seperates the items.
    		@param  ie Whether or not every item is enclosed.
    		@param  first_n  Only the given number of lines are read, starting from the beginning of the file.


				@exception Exception::FileNotFound is thrown if the file could not be opened.
    	*/
			CsvFile(const String& filename, char is = ',',bool ie = false, Int first_n = -1);

    	/**
    		@brief Loads data from a text file.
    	
    		@param  filename The input file name.
    		@param  is character which seperates the items.
    		@param  ie Whether or not every item is enclosed.
    		@param  first_n  Only the given number of lines are read, starting from the beginning of the file.

				@exception Exception::FileNotFound is thrown if the file could not be opened.
    	*/
			void fload(const String& filename, char is = ',', bool ie = false, Int first_n = -1);
			
    	/**
    		@brief writes all items from a row to list
    		
    		@param row the row which will be read
    		@param list StringList which will contain all items of the row
    		
    		@exception Exception::InvalidIterator is thrown if the row is not existing
    		
    		@return  returns false if the given row could not be seperated into items
    		
    	*/
    	bool getRow(Size row, StringList& list);
			
		private:
			char itemseperator_;
			bool itemenclosed_;

			
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_CSVFILE_H
