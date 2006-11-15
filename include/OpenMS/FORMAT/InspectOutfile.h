// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_INSPECTOUTFILE_H
#define OPENMS_FORMAT_INSPECTOUTFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace OpenMS
{
	/**
		@brief Representation of an Inspect outfile
		
		This class serves to read in an Inspect outfile and write an AnalysisXML file
		
		@todo write test(Martin)
		
		@ingroup FileIO
	*/
	class InspectOutfile
	{
		public:
			/// Constructor
			InspectOutfile();
			
			/// load the results of an InsPecT search
			std::vector< UnsignedInt >
			load(
				const std::string& result_filename,
				std::vector< Identification >& identifications,
				ProteinIdentification& protein_identification,
				std::vector< Real >& precursor_retention_times,
				std::vector< Real >& precursor_mz_values,
				Real p_value_threshold)
//				const std::string& database_filename)
			throw (
				Exception::FileNotFound,
				Exception::ParseError,
				Exception::IllegalArgument);
			
			std::vector< UnsignedInt >
			getWantedRecords(
				const std::string& result_filename,
				Real p_value_threshold)
			throw (
				Exception::FileNotFound,
				Exception::ParseError,
				Exception::IllegalArgument);
			
			void
			compressTrieDB(
				const std::string& database_filename,
				const std::string& index_filename,
				const std::vector< UnsignedInt >& wanted_records,
				const std::string& snd_database_filename,
				const std::string& snd_index_filename,
				bool append = false)
			throw (
				Exception::FileNotFound,
				Exception::ParseError,
				Exception::UnableToCreateFile);
			
			void
			generateTrieDB(
				const std::string& source_database_filename,
				const std::string& database_filename,
				const std::string& index_filename,
				bool append = false,
				const std::string species = "")
			throw (
				Exception::FileNotFound,
				Exception::ParseError,
				Exception::UnableToCreateFile);
			
		protected:
			/// get the accession and accession type of a protein
			void getACAndACType(String line, std::string& accession, std::string& accession_type) throw (Exception::ParseError);
			
			/// given a vector of peptide hits, either insert the new peptide hit or update its ProteinHits, returns whether an update took place
			bool updatePeptideHits(PeptideHit& peptide_hit, std::vector< PeptideHit >& peptide_hits);
			
			void getPrecursorRTandMZ(
				const std::vector< std::pair< String, std::vector< UnsignedInt > > >& files_and_scan_numbers,
				std::vector< Real >& precursor_retention_times,
				std::vector< Real >& precursor_mz_values,
				UnsignedInt scans)
			throw(
				Exception::ParseError);
			
			void
			getLabels(
				const std::string& source_database_filename,
				std::string& ac_label,
				std::string& sequence_start_label,
				std::string& sequence_end_label,
				std::string& comment_label,
				std::string& species_label)
			throw (
				Exception::FileNotFound,
				Exception::ParseError);
			
			std::vector< UnsignedInt >
			getSequences(
				const std::string& database_filename,
				const std::map< UnsignedInt, UnsignedInt >& wanted_records, // < record number, number of protein in a vector >
				std::vector< String >& sequences)
			throw (
				Exception::FileNotFound);
			
			void
			getLabels(
				const std::string& source_database_filename,
				const std::string& ac_label,
				const std::string& sequence_start_label,
				const std::string& sequence_end_label,
				const std::string& comment_label,
				const std::string& species_label)
			throw (
				Exception::FileNotFound,
				Exception::ParseError);
			
			static const UnsignedInt db_pos_length_;
			static const UnsignedInt trie_db_pos_length_;
			static const UnsignedInt protein_name_length_;
			static const UnsignedInt record_length_;
			static const char trie_delimiter_;
			static const std::string score_type_;
	};
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTOUTFILE_H
