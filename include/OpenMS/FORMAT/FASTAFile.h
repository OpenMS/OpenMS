// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FASTAFILE_H
#define OPENMS_FORMAT_FASTAFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
		 @brief This class serves for reading in FASTA files

  */
  class FASTAFile
  {
	 public:

		/**@brief
			 FASTA entry type (first is comment, second is sequence)

			 The first String in the pair corresponds to the comment that is
			 written after the > in the FASTA file and the second String in the
			 pair corresponds to the sequence.

		*/
		typedef std::pair<String, String> FASTAEntry;

		/// FASTA type (a sequence of FASTA entries)
		typedef std::vector< FASTAEntry > FASTAType;

		/// Copy constructor
		FASTAFile();

		/// Destructor
		virtual ~FASTAFile();

		/**
			 @brief loads a FASTA file given by 'filename' and stores the information in 'data'
		*/
		void load(const String& filename, FASTAType& data) throw (Exception::FileNotFound,Exception::ParseError);

		/**
			 @brief stores the data given by 'data' at the file 'filename'
		*/
		void store(const String& filename, const FASTAType& data) const throw (Exception::UnableToCreateFile);

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_FASTAFILE_H
