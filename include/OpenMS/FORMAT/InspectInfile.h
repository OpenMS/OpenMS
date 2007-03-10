// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_INSPECTINFILE_H
#define OPENMS_FORMAT_INSPECTINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>
#include <fstream>


namespace OpenMS
{
	/**
		@brief Inspect input file adapter.
		
		Creates a file that can be used for Inspect search from a peak list.
  	
  	@ingroup FileIO
	*/
  class InspectInfile
  {
		public:
			
			/// default constructor
			InspectInfile();

			/// copy constructor
			InspectInfile(const InspectInfile& inspect_infile);

			/// destructor
			virtual ~InspectInfile();

			/// assignment operator
			InspectInfile& operator=(const InspectInfile& inspect_infile);

			/// stores the experiment data in an Inspect input file that can be used as input for Inspect shell execution
			void store(const String& filename) throw (Exception::UnableToCreateFile);
	    
			/**
				@brief Specifies a spectrum file to search.
				
				You can specify the name of a directory to search every file in that directory (non-recursively). Supported spectra file formats are .mzXML, .mzData, .ms2, dta, and .pkl. Multiple spectra in one .dta file are not supported.
			*/
			const std::string& getSpectra() const;
			void setSpectra(const std::string& spectra);
			
			/**
				@brief Specifies the name of a database (.trie file) to search.
				
				The .trie file contains one or more protein sequences delimited by asterisks, with no whitespace or other data. Use PrepDB.py (see above) to prepare a .trie file. Most .trie files have a corresponding .index file giving the names of the proteins. You can specify at most one database.
			*/
			const String& getDb() const;
			void setDb(const String& db);
			
			/// Specifies the name of a protease. "Trypsin", "None", and "Chymotrypsin" are the available values.
			const String& getProtease() const;
			void setProtease(const String& protease);
			
			/**
				@brief Specifies an amino acid modification.
				
				The delta mass (in daltons) and affected amino acids are required. The first four characters of the name should be unique. Valid values for "type" are "fix", "cterminal", "nterminal", and "opt" (the default). Examples:
				mod,+57,C,fix - Most searches should include this line; it reflects CAM (carbamidomethylation) reaction which prevents cysteines from forming disulfide bonds.
				mod,80,STY,opt,phosphorylation
				mod,16,M (Oxidation of methionine, seen in many samples)
				mod,43,*,nterminal (N-terminal carbamylation, common if sample is treated with urea)
				Important note: When searching for phosphorylation sites, use a modification with the name "phosphorylation". This lets Inspect know that it should use its model of phosphopeptide fragmentation when generating tags and scoring matches. (Phosphorylation of serine dramatically affects fragmentation, so modeling it as simply an 80Da offset is typically not sufficient to detect sites with high sensitivity)
			*/
			const std::vector< std::vector< String > >& getMod() const;
			std::vector< std::vector< String > >& getMod();
			void setMod(const std::vector< std::vector< String > >& mod);
			void addMod(const std::vector< String >& mod);
			
			/// Number of PTMs permitted in a single peptide. Defaults to 2.
			const int getMods() const;
			void setMods(int mods);
			
			/**
				@brief run Inspect in a blind mode
				
				If true, use the MS-Alignment algorithm to perform a blind search (allowing arbitrary modification masses). Running a blind search with one mod per peptide is slower than the normal (tag-based) search; running time is approximately 1 second per spectra per megabyte of database. Running a blind search with two mods is significantly slower. We recommend performing "blind" searches against a small database, containing proteins output by an earlier search.
			*/
			const unsigned int getBlind() const;
			void setBlind(unsigned int blind);
			
			/**
				@brief the maximum modification size (in Da) to consider in a blind search
				
				Defaults to 200. Larger values require more time to search.
			*/
			const DoubleReal getMaxPTMsize() const;
			void setMaxPTMsize(DoubleReal maxptmsize);
			
			/**
				@brief Specifies the parent mass tolerance, in Daltons.
				
				A candidate's flanking mass can differ from the tag's flanking mass by no more than this amount.
			*/
			const DoubleReal getPMTolerance() const;
			void setPMTolerance(DoubleReal PM_tolerance);
			
			/**
				@brief How far b and y peaks can be shifted from their expected masses.
				
				Default is 0.5. Higher values produce a more sensitive but much slower search.
			*/
			const DoubleReal getIonTolerance() const;
			void setIonTolerance(DoubleReal ion_tolerance);
			
			/// If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible.
			const unsigned int getMulticharge() const;
			void setMulticharge(unsigned int multicharge);
			
			/// If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.
			const String& getInstrument() const;
			void setInstrument(const String& instrument);
			
			/// Number of tags to generate.
			const int getTagCount() const;
			void setTagCount(int TagCount);
			
		private:

			std::string spectra_; ///< Specifies a spectrum file to search. You can specify the name of a directory to search every file in that directory (non-recursively). Supported spectra file formats are .mzXML, .mzData, .ms2, dta, and .pkl. Multiple spectra in one .dta file are not supported.

    	String db_; ///< Specifies the name of a database (.trie file) to search. The .trie file contains one or more protein sequences delimited by asterisks, with no whitespace or other data. Use PrepDB.py (see above) to prepare a .trie file. Most .trie files have a corresponding .index file giving the names of the proteins. You can specify at most one database.

			String protease_; ///< Specifies the name of a protease. "Trypsin", "None", and "Chymotrypsin" are the available values.

			std::vector< std::vector< String > > mod_; ///< Specifies an amino acid modification. The delta mass (in daltons) and affected amino acids are required. The first four characters of the name should be unique. Valid values for "type" are "fix", "cterminal", "nterminal", and "opt" (the default). Examples:

			int mods_; ///< Number of PTMs permitted in a single peptide. Defaults to 2. -1 is not set

			unsigned int blind_; ///< If true, use the MS-Alignment algorithm to perform a blind search (allowing arbitrary modification masses). Running a blind search with one mod per peptide is slower than the normal (tag-based) search; running time is approximately 1 second per spectra per megabyte of database. Running a blind search with two mods is significantly slower. We recommend performing "blind" searches against a small database, containing proteins output by an earlier search. (The "Summary.py" script can be used to generate a second-pass database from initial search results)
			/// 0 - false, 1 - true, 2 - not set

		  	DoubleReal maxptmsize_; ///< For blind search, specifies the maximum modification size (in Da) to consider. Defaults to 200. Larger values require more time to search. <0 is not set

			DoubleReal PM_tolerance_; ///< Specifies the parent mass tolerance, in Daltons. A candidate's flanking mass can differ from the tag's flanking mass by no more than ths amount. <0 is not set

			DoubleReal ion_tolerance_; ///< How far b and y peaks can be shifted from their expected masses. Default is 0.5. Higher values produce a more sensitive but much slower search. <0 is not set
			
			unsigned int multicharge_; ///< If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible.
			/// 0 - false, 1 - true, 2 - not set

			String instrument_; ///< If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.

			int tag_count_; ///< Number of tags to generate. <0 is not set
		
  };

} /// namespace OpenMS

#endif /// OPENMS_FORMAT_INSPECTINFILE_H
