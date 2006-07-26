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
// $Id: InspectInfile.h,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
// $Author: martinlangwisch $
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_INSPECTINFILE_H
#define OPENMS_FORMAT_INSPECTINFILE_H

#include <OpenMS/FORMAT/InspectFile.h>

#include <map>
#include <set>
#include <stdlib.h>


namespace OpenMS
{
	/**
		@brief Inspect input file adapter.
		
		Creates a file that can be used for Inspect search from a peak list.
  
  	@ingroup FileIO
	*/
  class InspectInfile:
		public InspectFile
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

			/// generate a database from an Inspect result file; this new database can be used for a blind search
			/// either a trie database is given as input, or a corresponding database is generated from the input
			void generateSecondDatabase(const std::string& result_filename_, const std::string& result_path, const std::string& database_path, const std::string& database_filename_, double cutoff_p_value, int min_annotated_spectra_per_protein, std::string second_database_filename_, std::string second_index_filename_, std::string index_filename_, std::string second_database_path, std::string species = "None") throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument);

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
			
			/**
				@brief Specifies the name of a FASTA-format protein database to search.
				
				If you plan to search a large database, it is more efficient to preprocess it using PrepDB.py and use the "db" command instead. You can specify at most one sequence file, but you can include both one "db" and one "sequence_file".
			*/
			const String& getSequenceFile() const;
			void setSequenceFile(const String& sequence_file);
			
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
			const double getMaxPTMsize() const;
			void setMaxPTMsize(double maxptmsize);
			
			/**
				@brief Specifies the parent mass tolerance, in Daltons.
				
				A candidate's flanking mass can differ from the tag's flanking mass by no more than ths amount.
			*/
			const double getPMTolerance() const;
			void setPMTolerance(double PM_tolerance);
			
			/**
				@brief How far b and y peaks can be shifted from their expected masses.
				
				Default is 0.5. Higher values produce a more sensitive but much slower search.
			*/
			const double getIonTolerance() const;
			void setIonTolerance(double ion_tolerance);
			
			/**
				@brief Use the file to specify PTM frequencies, for use in tag generation.
				
				This is more accurate tagging than than the default behavior (where tags can contain any PTM), but requires the creation of the jump frequency file.
			*/
			const String& getJumpscores() const;
			void setJumpscores(const String& jumpscores);
			
			/// If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible.
			const unsigned int getMulticharge() const;
			void setMulticharge(unsigned int multicharge);
			
			/// If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.
			const String& getInstrument() const;
			void setInstrument(const String& instrument);
			
			/// Number of tags to generate for the first pass of a two-pass search.
			const int getTagCountA() const;
			void setTagCountA(int TagCountA);
			
			/// Number of tags to generate for the second pass of a two-pass search, OR the number of tags to use in a one-pass search.
			const int getTagCountB() const;
			void setTagCountB(int TagCountB);
			
			/// Use two-pass search. The first pass uses fewer tags, and produces a list of proteins which are re-searched in the second pass.
			const unsigned int getTwopass() const;
			void setTwopass(unsigned int twopass);
			
		private:

			std::string spectra_; ///< Specifies a spectrum file to search. You can specify the name of a directory to search every file in that directory (non-recursively). Supported spectra file formats are .mzXML, .mzData, .ms2, dta, and .pkl. Multiple spectra in one .dta file are not supported.

    	String db_; ///< Specifies the name of a database (.trie file) to search. The .trie file contains one or more protein sequences delimited by asterisks, with no whitespace or other data. Use PrepDB.py (see above) to prepare a .trie file. Most .trie files have a corresponding .index file giving the names of the proteins. You can specify at most one database.

    	String sequence_file_; ///< Specifies the name of a FASTA-format protein database to search. If you plan to search a large database, it is more efficient to preprocess it using PrepDB.py and use the "db" command instead. You can specify at most one sequence file, but you can include both one "db" and one "sequence_file".

			String protease_; ///< Specifies the name of a protease. "Trypsin", "None", and "Chymotrypsin" are the available values.

			std::vector< std::vector< String > > mod_; ///< Specifies an amino acid modification. The delta mass (in daltons) and affected amino acids are required. The first four characters of the name should be unique. Valid values for "type" are "fix", "cterminal", "nterminal", and "opt" (the default). Examples:

			int mods_; ///< Number of PTMs permitted in a single peptide. Defaults to 2. -1 is not set

			unsigned int blind_; ///< If true, use the MS-Alignment algorithm to perform a blind search (allowing arbitrary modification masses). Running a blind search with one mod per peptide is slower than the normal (tag-based) search; running time is approximately 1 second per spectra per megabyte of database. Running a blind search with two mods is significantly slower. We recommend performing "blind" searches against a small database, containing proteins output by an earlier search. (The "Summary.py" script can be used to generate a second-pass database from initial search results)
			/// 0 - false, 1 - true, 2 - not set

		  	double maxptmsize_; ///< For blind search, specifies the maximum modification size (in Da) to consider. Defaults to 200. Larger values require more time to search. <0 is not set

			double PM_tolerance_; ///< Specifies the parent mass tolerance, in Daltons. A candidate's flanking mass can differ from the tag's flanking mass by no more than ths amount. <0 is not set

			double ion_tolerance_; ///< How far b and y peaks can be shifted from their expected masses. Default is 0.5. Higher values produce a more sensitive but much slower search. <0 is not set

			String jumpscores_; ///< Use the file to specify PTM frequencies, for use in tag generation. This is more accurate tagging than than the default behavior (where tags can contain any PTM), but requires the creation of the jump frequency file.

			unsigned int multicharge_; ///< If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible.
			/// 0 - false, 1 - true, 2 - not set

			String instrument_; ///< If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass.

			int tag_count_a_; ///< Number of tags to generate for the first pass of a two-pass search. <0 is not set

			int tag_count_b_; ///< Number of tags to generate for the second pass of a two-pass search, OR the number of tags to use in a one-pass search. <0 is not set

			unsigned int twopass_; ///< Use two-pass search. The first pass uses fewer tags, and produces a list of proteins which are re-searched in the second pass.
			/// 0 - false, 1 - true, 2 - not set
		
  };

} /// namespace OpenMS

#endif /// OPENMS_FORMAT_INSPECTINFILE_H
