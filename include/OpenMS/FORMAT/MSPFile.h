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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MSPFILE_H
#define OPENMS_FORMAT_MSPFILE_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief File adapter for MSP files (NIST spectra library)


		@htmlinclude OpenMS_MSPFile.parameters

		@ingroup FileIO
	*/
	class OPENMS_DLLAPI MSPFile : public DefaultParamHandler
	{
		public:

			/** Constructors and destructors
			*/
			//@{
			///Default constructor
			MSPFile();

			/// Copy constructor
			MSPFile(const MSPFile& rhs);

			///Destructor
			virtual ~MSPFile();
			//@}

			/// assignment operator
			MSPFile& operator = (const MSPFile& rhs);

			/**
				@brief Loads a map from a MSPFile file.

				@param exp RichPeakMap which contains the spectra after reading
				@param filename the filename of the experiment
				@param ids output parameter which contains the peptide identifications from the spectra anntations

				@throw FileNotFound is thrown if the file could not be found
				@throw ParseError is thrown if the given file could not be parsed
				@throw ElementNotFound is thrown if a annotated modification cannot be found in ModificationsDB (PSI-MOD definitions)
			*/
			void load(const String& filename, std::vector<PeptideIdentification>& ids, RichPeakMap& exp);

			/**
				@brief Stores a map in a MSPFile file.

				@throw UnableToCreateFile is thrown if the given file could not be created
			*/
			void store(const String& filename, const RichPeakMap& exp) const;

			protected:

				/// reads the header information and stores it as metainfo in the spectrum
				void parseHeader_(const String& header, RichPeakSpectrum& spec);

	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_MSPFILE_H
