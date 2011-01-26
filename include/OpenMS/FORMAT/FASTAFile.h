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
// $Maintainer: Sandro Andreotti $
// $Authors: Nico Pfeifer $
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
  class OPENMS_DLLAPI FASTAFile
  {
	 public:

		/**@brief
			 FASTA entry type (identifier, description and sequence)

			 The first String corresponds to the identifier that is
			 written after the > in the FASTA file. The part after the 
			 first whitespace is stored in description and the text 
			 from the next line until the next > (exclusive) is stored
			 in sequence. 

		*/
		struct FASTAEntry
		{
			String identifier;
			String description;
			String sequence;

			FASTAEntry()
				: identifier(""),
					description(""),
					sequence("")
			{
      }

			FASTAEntry(String id, String desc, String seq)
				: identifier(id),
					description(desc),
					sequence(seq)
			{
      }

			bool operator == (const FASTAEntry& rhs) const
			{
				return identifier == rhs.identifier
					&& description == rhs.description
					&& sequence == rhs.sequence;
			}
				
		};

		/// Copy constructor
		FASTAFile();

		/// Destructor
		virtual ~FASTAFile();

		/**
			 @brief loads a FASTA file given by 'filename' and stores the information in 'data'

			 @exception Exception::FileNotFound is thrown if the file does not exists.
			 @exception Exception::ParseError is thrown if the file does not suit to the standard.
		*/
		void load(const String& filename, std::vector<FASTAEntry>& data);

		/**
			 @brief stores the data given by 'data' at the file 'filename'

			 @exception Exception::UnableToCreateFile is thrown if the process is not able to write the file.
		*/
		void store(const String& filename, const std::vector<FASTAEntry>& data) const;

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_FASTAFILE_H
