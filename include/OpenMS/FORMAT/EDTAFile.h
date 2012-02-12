// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_EDTAFILE_H
#define OPENMS_FORMAT_EDTAFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>
#include <vector>

namespace OpenMS
{
 	/**
 		@brief File adapter for Enhanced DTA files.
 		
    Input text file containing tab, space or comma separated columns.
    The separator between columns is checked in the first line in this order.

    It supports three variants of this format.

    - Columns are: RT, MZ, Intensity
      A header is optional.

    - Columns are: RT, MZ, Intensity, Charge, <Meta-Data> columns{0,}
      A header is mandatory.

      Example:
      @code
      RT m/z Intensity charge mymeta1 mymeta2
      321 405.233 24543534 2 lala  lili
      321 406.207 4343344  2 blubb blabb
      @endcode

    - Columns are: (RT, MZ, Intensity, Charge){1,}, <Meta-Data> columns{0,}
      Header is mandatory.
      First quadruplet is the consensus. All following quadruplets describe the sub-features.
      This variant is discerned from variant #2 by the name of the fifth column, which is required to be RT1 (or rt1).
      All other column names for sub-features are faithfully ignored.

      Example:
      @code
      RT MZ INT CHARGE RT1 MZ1 INT1 CHARGE1 RT2 MZ2 INT2 CHARGE2
      321 405 100 2 321 405 100 2 321 406 50 2
      323 406 200 2 323 406 200 2 323 407 100 2 323 407 50 2
      @endcode

  	@ingroup FileIO
  */
  class OPENMS_DLLAPI EDTAFile
  {
    public:
      /// Default constructor
      EDTAFile();
			/// Destructor
      virtual ~EDTAFile();
      
    private:
      /**
       * Check if column exists and convert String into DoubleReal.
       */
      DoubleReal checkedToDouble_(const std::vector<String> &parts, Size index, DoubleReal def = -1);

      /**
       * Check if column exists and convert String into Int.
       */
      Int checkedToInt_(const std::vector<String> &parts, Size index, Int def = -1);

    public:
      /**
        @brief Loads a EDTA file into a consensusXML.
 				
 				The content of the file is stored in @p features.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
      */
      void load(const String& filename, ConsensusMap& consensus_map);

      /**
      	@brief Stores a ConsensusMap as an enhanced DTA file.
      	
        NOT IMPLEMENTED

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
      */
      void store(const String& filename, const ConsensusMap& map) const;
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_EDTAFILE_H

