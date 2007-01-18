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

#ifndef OPENMS_FORMAT_SEQUESTINFILE_H
#define OPENMS_FORMAT_SEQUESTINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>


namespace OpenMS
{
	/**
		@brief Sequest input file adapter.
		
		Creates a sequest.params file for Sequest search from a peak list.
  	  	
  	@ingroup FileIO
	*/
  class SequestInfile
  {
		public:
			/// default constructor
			SequestInfile();

			/// copy constructor
			SequestInfile(const SequestInfile& sequest_infile);

			/// destructor
			virtual ~SequestInfile();

			/// assignment operator
			SequestInfile& operator=(const SequestInfile& sequest_infile);

			/// stores the experiment data in a Sequest input file that can be used as input for Sequest shell execution
			void store(const String& filename) throw (Exception::UnableToCreateFile);

			void setDynMod(char amino_acid, Real mass);

			/// returns the enzyme used for cleavage
			const String& getEnzymeInfo() const;
			/// sets the enzyme used for cleavage
			void setEnzymeInfo(String& value);

			/// returns the used database
			const String& getDatabase() const;
			/// sets the used database
			void setDatabase(const String& value);

			/// returns the snd database used
			const String& getSndDatabase() const;
			/// sets the second database used
			void setSndDatabase(const String& value);

			/// returns whether neutral losses are considered for the a-, b- and y-ions
			const String& getNeutralLossesForIons() const;
			/// sets whether neutral losses are considered for the a-, b- and y-ions
			void setNeutralLossesForIons(const String& neutral_losses_for_ions);

			/// returns the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
			const String& getIonSeriesWeights() const;
			/// sets the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
			void setIonSeriesWeights(const String& ion_series_weights);

			/// returns the dynamic modifications
			const String& getDynMods() const;
			/// sets the dynamic modifications
			void setDynMods(const String& dyn_mods);

			/// returns the partial sequences (space delimited) that have to occur in the theortical spectra
			const String& getPartialSequence() const;
			/// sets the partial sequences (space delimited) that have to occur in the theortical spectra
			void setPartialSequence(const String& partial_sequence);

			/// returns the sequences (space delimited) that have to occur, or be absent (preceeded by a tilde) in the header of a protein to be considered
			const String& getSequenceHeaderFilter() const;
			/// sets the sequences (space delimited) that have to occur, or be absent (preceeded by a tilde) in the header of a protein to be considered
			void setSequenceHeaderFilter(const String& sequence_header_filter);

			/// returns the protein mass filter (either min and max mass, or mass and tolerance value in percent)
			const String& getProteinMassFilter() const;
			/// sets the protein mass filter (either min and max mass, or mass and tolerance value in percent)
			void setProteinMassFilter(const String& protein_mass_filter);


			/// returns the peptide mass tolerance
			Real getPeptideMassTolerance() const;
			/// sets the peptide mass tolerance
			void setPeptideMassTolerance(Real value);

			/// returns the fragment ion tolerance
			Real getFragmentIonTolerance() const;
			/// sets the fragment ion tolerance
			void setFragmentIonTolerance(Real value);

			/// returns the match peak tolerance
			Real getMatchPeakTolerance() const;
			/// sets the match peak tolerance
			void setMatchPeakTolerance(Real value);

			/// returns the the cutoff of the ratio matching theoretical peaks/theoretical peaks
			Real getIonCutoffPercentage() const;
			/// sets the ion cutoff of the ratio matching theoretical peaks/theoretical peaks
			void setIonCutoffPercentage(Real value);

			/// returns the dynamic modification for the N-terminal of the peptide
			Real getDynNTermMod() const;
			/// sets the dynamic modification for the N-terminal of the peptide
			void setDynNTermMod(Real mass);

			/// returns the dynamic modification for the C-terminal of the peptide
			Real getDynCTermMod() const;
			/// sets the dynamic modification for the C-terminal of the peptide
			void setDynCTermMod(Real mass);

			/// returns the static modification for the N-terminal of the peptide
			Real getStatNTermMod() const;
			/// sets the static modification for the N-terminal of the peptide
			void setStatNTermMod(Real mass);

			/// returns the static modification for the C-terminal of the peptide
			Real getStatCTermMod() const;
			/// sets the static modification for the C-terminal of the peptide
			void setStatCTermMod(Real mass);

			/// returns the static modification for the N-terminal of the protein
			Real getStatNTermProtMod() const;
			/// sets the static modification for the N-terminal of the protein
			void setStatNTermProtMod(Real mass);

			/// returns the static modification for the C-terminal of the protein
			Real getStatCTermProtMod() const;
			/// sets the static modification for the C-terminal of the protein
			void setStatCTermProtMod(Real mass);


			/// returns the peptide mass unit
			SignedInt getPeptideMassUnit() const;
			/// sets the peptide mass unit
			void setPeptideMassUnit(SignedInt value);

			/// return the number of peptides to be displayed
			SignedInt getOutputLines() const;
			/// sets the number of peptides to be displayed
			void setOutputLines(SignedInt value);

			/// returns the enzyme used for cleavage (by means of the number from a list of enzymes)
			SignedInt getEnzymeNumber() const;
			/// sets the enzyme used for cleavage (by means of the number from a list of enzymes)
			SignedInt setEnzymeNumber(SignedInt value);

			/// returns the maximum number of amino acids containing the same modification in a peptide
			SignedInt getMaxAAPerModPerPeptide() const;
			/// sets the maximum number of amino acids containing the same modification in a peptide
			void setMaxAAPerModPerPeptide(SignedInt value);

			/// returns the maximum number of modifications that are allowed in a peptide
			SignedInt getMaxModsPerPeptide() const;
			/// set the maximum number of modifications that are allowed in a peptide
			void setMaxModsPerPeptide(SignedInt value);

			/// returns the nucleotide reading frame
			SignedInt getNucleotideReadingFrame() const;
			/// sets the nucleotide reading frame:
			///		0 	The FASTA file contains amino acid codes. No translation is needed. This is the best and fastest case.
			///		1 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the first DNA code.
			///		2 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the second DNA code.
			///		3 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the third DNA code.
			///		4 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the first DNA code.
			///		5 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the second DNA code.
			///		6 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the third DNA code.
			///		7 	Use each of the DNA translations of the codes 1, 2, 3.
			///		8 	Use each of the DNA translations of the codes 4, 5, 6.
			///		9 	Use each of the DNA translations of the codes 1, 2, 3, 4, 5, 6.
			void setNucleotideReadingFrame(SignedInt value);

			/// returns the maximum number of internal cleavage sites
			SignedInt getMaxInternalCleavageSites() const;
			/// sets the maximum number of internal cleavage sites
			void setMaxInternalCleavageSites(SignedInt value);

			/// returns the number of top abundant peaks to match with theoretical ones
			SignedInt getMatchPeakCount() const;
			/// sets the number of top abundant peaks to with theoretical ones
			void setMatchPeakCount(SignedInt value);

			/// returns the number of top abundant peaks that are allowed not to match with a theoretical peak
			SignedInt getMatchPeakAllowedError() const;
			/// sets the number of top abundant peaks that are allowed not to match with a theoretical peak
			void setMatchPeakAllowedError(SignedInt value);


			/// returns whether fragment ions shall be displayed
			bool getShowFragmentIons() const;
			/// sets whether fragment ions shall be displayed
			void setShowFragmentIons(bool value);

			/// returns whether all proteins containing a found peptide should be displayed
			bool getPrintDuplicateReferences() const;
			/// sets whether all proteins containing a found peptide should be displayed
			void setPrintDuplicateReferences(bool value);

// 			bool getUsePhosphoFragmentation() const;
// 			void setUsePhosphoFragmentation(bool value);

			/// return whether peaks near (15 amu) the precursor peak are removed
			bool getRemovePrecursorNearPeaks() const;
			/// sets whether peaks near (15 amu) the precursor peak are removed
			void setRemovePrecursorNearPeaks(bool value);

			/// return the mass type of the parent (0 - monoisotopic, 1 - average mass)
			bool getMassTypeParent() const;
			/// sets the mass type of the parent (0 - monoisotopic, 1 - average mass)
			void setMassTypeParent(bool value);

			/// return the mass type of the fragments (0 - monoisotopic, 1 - average mass)
			bool getMassTypeFragment() const;
			/// sets the mass type of the fragments (0 - monoisotopic, 1 - average mass)
			void setMassTypeFragment(bool value);

			/// returns whether normalized xcorr values are displayed
			bool getNormalizeXcorr() const;
			/// sets whether normalized xcorr values are displayed
			void setNormalizeXcorr(bool value);

			/// returns whether residues are in upper case
			bool getResiduesInUpperCase() const;
			/// sets whether residues are in upper case
			void setResiduesInUpperCase(bool value);


			/// adds an enzyme to the list and sets is as used
			/// the vector constists of four strings:
			/// name, cut direction: 0 (N to C) / 1, cuts after (list of aa), doesn't cut before (list of aa)
			void addEnzymeInfo(std::vector< String >& enzyme_info);

			/// returns the static modifications (map of amino acids and corresponding modification)
			const std::map< char, Real >& getStatMods() const;
			/// set the static modification for an amino acid
			char setStatMod(String amino_acid, Real mass);

		protected:
			/// returns some standard enzymes (used to initialize the enzyme list)
			String getStandardEnzymeInfo();

			/// the amino acids in one-letter-code
			static const String aas_single_letter_;// = "GASPVTCLIXNOBDQKZEMHFRYW";

			/// the static modifications (map of amino acids and corresponding modification)
			std::map< char, Real > stat_mods_;

			String enzyme_info_; ///< an endline-delimited list of enzymes; each with cutting direction 0 (N to C) /1; cuts after (list of aa); doesn't cut before (list of aa); the attributes are tab-delimited
			String database_; ///< database used
			String snd_database_; ///< second database used
			String neutral_losses_for_ions_; ///< whether neutral losses are considered for the a-; b- and y-ions (e.g. 011 for b- and y-ions)
			String ion_series_weights_;///< weights for the a-; b-; c-; d-; v-; w-; x-; y- and z-ion series; space delimited
			String dyn_mods_; ///< space-delimited list of dynamic modifications; each with weight and aas (space delimited)
			String partial_sequence_; ///< space-delimited list of sequence parts that have to occur in the theortical spectra
			String sequence_header_filter_;///< space-delimited list of sequences that have to occur or be absend (preceeded by a tilde) in a protein header; to be considered
			String protein_mass_filter_;
			
			
			Real peptide_mass_tolerance_;///< tolerance for matching a theoretical to an experimental peptide
			Real fragment_ion_tolerance_;///< tolerance for matching a theoretical to an experimental peak
			Real match_peak_tolerance_;///< minimum distance between two experimental peaks
			Real ion_cutoff_percentage_;///< cutoff of the ratio matching theoretical peaks/theoretical peaks
			Real dyn_n_term_mod_;///< dynamic modifications for the N-terminal of a peptide
			Real dyn_c_term_mod_;///< dynamic modifications for the C-terminal of a peptide
			Real stat_n_term_mod_;///< static modifications for the N-terminal of a peptide
			Real stat_c_term_mod_;///< static modifications for the C-terminal of a peptide
			Real stat_n_term_prot_mod_;///< static modifications for the N-terminal of a protein
			Real stat_c_term_prot_mod_;///< static modifications for the C-terminal of a protein
			
			
			SignedInt peptide_mass_unit_;///< peptide mass unit (0 = amu; 1 = mmu; 2 = ppm)
			SignedInt output_lines_;///< number of peptides to be displayed
			SignedInt enzyme_number_;///< number of the enzyme used for cleavage
			SignedInt highest_enzyme_number_;///< highest enzyme number
			SignedInt max_AA_per_mod_per_peptide_;///< maximum number of amino acids containing the same modification in a peptide
			SignedInt max_mods_per_peptide_;///< maximum number of modifications per peptide
			SignedInt nucleotide_reading_frame_;///< nucleotide reading frame:
					/// 0 	The FASTA file contains amino acid codes. No translation is needed. This is the best and fastest case.
					/// 1 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the first DNA code.
					/// 2 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the second DNA code.
					/// 3 	The DNA sequence is scanned left to right (forward direction). The amino acid code starts with the third DNA code.
					/// 4 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the first DNA code.
					/// 5 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the second DNA code.
					/// 6 	The DNA sequence is scanned right to left (backward direction for the complementary strand). The amino acid code starts with the third DNA code.
					/// 7 	Use each of the DNA translations of the codes 1; 2; 3.
					/// 8 	Use each of the DNA translations of the codes 4; 5; 6.
					/// 9 	Use each of the DNA translations of the codes 1; 2; 3; 4; 5; 6.
			SignedInt max_internal_cleavage_sites_;///< maximum number of internal cleavage sites
			SignedInt match_peak_count_;///< number of the top abundant peaks to match with theoretical one
			SignedInt match_peak_allowed_error_;///< number of peaks that may lack this test
			
			
			bool show_fragment_ions_;///< wether to display fragment ions
			bool print_duplicate_references_;///< whether all proteins containing a found peptide should be displayed
//		bool use_phospho_fragmentation_;///< 
			bool remove_precursor_near_peaks_;///< whether peaks near (15 amu) the precursor peak are removed
			bool mass_type_parent_;///< mass type of the parent peak (0 - monoisotopic; 1 - average)
			bool mass_type_fragment_;///< mass type of fragment peaks (0 - monoisotopic; 1 - average)
			bool normalize_xcorr_;///< whether to display normalized xcorr values
			bool residues_in_upper_case_;///< whether residues are in upper case
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_SEQUESTINFILE_H
