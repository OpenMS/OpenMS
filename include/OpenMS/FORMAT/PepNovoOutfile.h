// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_PEPNOVOOUTFILE_H
#define OPENMS_FORMAT_PEPNOVOOUTFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <map>

namespace OpenMS
{
	class ProteinIdentification;

	/**
		@brief Representation of a PepNovo output file
	
		This class serves to read in a PepNovo outfile. The information can be
		retrieved via the load function.

		@ingroup FileIO
	*/
	class OPENMS_DLLAPI PepNovoOutfile
	{
		public:
			/// Constructor
			PepNovoOutfile();

			/// copy constructor
			PepNovoOutfile(const PepNovoOutfile& pepnovo_outfile);

			/// destructor
			virtual ~PepNovoOutfile();

			/// assignment operator
			PepNovoOutfile& operator=(const PepNovoOutfile& pepnovo_outfile);

			/// equality operator
			bool operator==(const PepNovoOutfile& pepnovo_outfile) const;

			 /**
				@brief loads data from a Inspect outfile

				@param result_filename the file to be loaded
				@param peptide_identifications the peptide identification
				@param protein_identification the protein identifications
				@param score_threshold cutoff threshold for the PepNovo score (PnvScr)
				@param dta_filenames_and_precursor_retention_times retention times

				@throw Exception::FileNotFound is thrown if the result file could not be found
				@throw Exception::ParseErro is thrown if the results file could not be parsed

				This class serves to read in a PepNovo outfile. The information can be
				retrieved via the load function.
			*/
			void load(const std::string& result_filename, std::vector< PeptideIdentification >& peptide_identifications, ProteinIdentification& protein_identification, const Real& score_threshold, const std::map< String, Real >& dta_filenames_and_precursor_retention_times);

			/** get the search engine and it's version from a file that is the output of PepNovo run without parameters

				@param pepnovo_output_without_parameters_filename
				@param protein_identification 

				@throw Exception::FileNotFound is thrown if the results file could not be found
			*/
			void getSearchEngineAndVersion(const String& pepnovo_output_without_parameters_filename, ProteinIdentification& protein_identification);

	};

} //namespace OpenMS

#endif // OPENMS_FORMAT_PEPNOVOOUTFILE_H

