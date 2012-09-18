// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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


      /**
      	@brief Stores a FeatureMap as an enhanced DTA file.
      	
        Creating columns: RT, m/z, intensity, charge

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
      */
      void store(const String& filename, const FeatureMap<>& map) const;
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_EDTAFILE_H

