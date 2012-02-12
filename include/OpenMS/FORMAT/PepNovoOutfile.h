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
// $Maintainer: Sandro Andreotti $
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
				@brief loads data from a PepNovo outfile

				@param result_filename the file to be loaded
				@param peptide_identifications the peptide identifications
				@param protein_identification the protein identification
				@param score_threshold cutoff threshold for the PepNovo score (PnvScr)
				@param id_rt_mz map the spectrum identifiers returned by PepNovo
				to the rt and mz values of the spectrum (used to map the identifications back to the spectra). key= &lt;PepNovo Id&gt;, value= &lt;pair&lt;rt,mz&gt; &gt;.
				For spectra not present in this map identifications cannot be mapped back.
				@param mod_id_map map the OpenMS id for modifications (FullId) to the ids returned by PepNovo key= &lt;PepNovo_key&gt;, value= &lt;OpenMS FullId&gt;
			*/
			void load(const std::string& result_filename, std::vector< PeptideIdentification >& peptide_identifications,
			    ProteinIdentification& protein_identification,
			    const DoubleReal& score_threshold,
			    const std::map< String, std::pair<DoubleReal,DoubleReal> >& id_rt_mz,
			    const std::map<String, String> &mod_id_map);

			/** @brief get the search engine version and search parameters from a PepNovo output file
			 *
			 * search parameters (precursor tolerance, peak mass tolerance, allowed modifications)are stored in the protein_identification.

				@param pepnovo_output_without_parameters_filename
				@param protein_identification
			*/
			void getSearchEngineAndVersion(const String& pepnovo_output_without_parameters_filename, ProteinIdentification& protein_identification);

	};

} //namespace OpenMS

#endif // OPENMS_FORMAT_PEPNOVOOUTFILE_H

