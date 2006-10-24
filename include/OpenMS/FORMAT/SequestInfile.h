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

#ifndef OPENMS_FORMAT_SEQUESTINFILE_H
#define OPENMS_FORMAT_SEQUESTINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>


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
			
			// default constructor
			SequestInfile();

			// copy constructor
			SequestInfile(const SequestInfile& sequest_infile);

			// destructor
			//virtual ~SequestInfile();

			// assignment operator
			SequestInfile& operator=(const SequestInfile& sequest_infile);

			// stores the experiment data in an Sequest input file that can be used as input for Sequest shell execution
			void store(const String& filename) throw (Exception::UnableToCreateFile);
			
			void setDynMod(char amino_acid, float mass);
			
			
			const String& getEnzymeInfo() const;
			void setEnzymeInfo(String& value);
			
			const String& getDatabase() const;
			void setDatabase(String& value);
			
			const String& getSndDatabase() const;
			void setSndDatabase(String& value);
			
			const String& getNeutralLossesForIons() const;
			void setNeutralLossesForIons(String& neutral_losses_for_ions);
			
			const String& getIonSeriesWeights() const;
			void setIonSeriesWeights(String& ion_series_weights);
			
			const String& getDynMods() const;
			void setDynMods(String& dyn_mods);
			
			const String& getPartialSequence() const;
			void setPartialSequence(String& partial_sequence);
			
			const String& getSequenceHeaderFilter() const;
			void setSequenceHeaderFilter(String& sequence_header_filter);
			
			
			float getPeptideMassTolerance() const;
			void setPeptideMassTolerance(float value);
			
			float getFragmentIonTolerance() const;
			void setFragmentIonTolerance(float value);
			
			float getMatchPeakTolerance() const;
			void setMatchPeakTolerance(float value);
			
			float getIonCutoffPercentage() const;
			void setIonCutoffPercentage(float value);
			
			float getMinimumProteinMass() const;
			void setMinimumProteinMass(float value);
			
			float getMaximumProteinMass() const;
			void setMaximumProteinMass(float value);
			
			float getDynNTermMod() const;
			void setDynNTermMod(float mass);
			
			float getDynCTermMod() const;
			void setDynCTermMod(float mass);
			
			float getStatNTermMod() const;
			void setStatNTermMod(float mass);
			
			float getStatCTermMod() const;
			void setStatCTermMod(float mass);
			
			float getStatNTermProtMod() const;
			void setStatNTermProtMod(float mass);
			
			float getStatCTermProtMod() const;
			void setStatCTermProtMod(float mass);
			
			
			int getPeptideMassUnits() const;
			void setPeptideMassUnits(int value);
			
			int getNumOutputLines() const;
			void setNumOutputLines(int value);
			
			int getEnzymeNumber() const;
			void setEnzymeNumber(int value);
			
			int getMaxNumDifAAPerMod() const;
			void setMaxNumDifAAPerMod(int value);
			
			int getMaxNumModsPerPeptide() const;
			void setMaxNumModsPerPeptide(int value);
			
			int getNucleotideReadingFrame() const;
			void setNucleotideReadingFrame(int value);
			
			int getMaxNumInternalCleavageSites() const;
			void setMaxNumInternalCleavageSites(int value);
			
			int getMatchPeakCount() const;
			void setMatchPeakCount(int value);
			
			int getMatchPeakAllowedError() const;
			void setMatchPeakAllowedError(int value);
			
			
			bool getShowFragmentIons() const;
			void setShowFragmentIons(bool value);
			
			bool getPrintDuplicateReferences() const;
			void setPrintDuplicateReferences(bool value);
			
			bool getUsePhosphoFragmentation() const;
			void setUsePhosphoFragmentation(bool value);
			
			bool getRemovePrecursorPeak() const;
			void setRemovePrecursorPeak(bool value);
			
			bool getMassTypeParent() const;
			void setMassTypeParent(bool value);
			
			bool getMassTypeFragment() const;
			void setMassTypeFragment(bool value);
			
			bool getNormalizeXcorr() const;
			void setNormalizeXcorr(bool value);
			
			bool getResiduesInUpperCase() const;
			void setResiduesInUpperCase(bool value);
			
			
			void addEnzymeInfo(std::vector< String >& enzyme_info);
			
			const std::map< char, float >& getStatMods() const;
			void setStatMods(std::map< char, float >& stat_mods);
			
			void setStatMod(char amino_acid, float mass);
		
		protected:
		
			String getStandardEnzymeInfo();
			
			const static String aas_single_letter_;// = "GASPVTCLIXNOBDQKZEMHFRYW";
			
			std::map< char, float > stat_mods_;
			
			String enzyme_info_,
						 database_,
						 snd_database_,
						 neutral_losses_for_ions_,
						 ion_series_weights_,
						 dyn_mods_,
						 partial_sequence_,
						 sequence_header_filter_;
			
			float peptide_mass_tolerance_,
						fragment_ion_tolerance_,
						dyn_n_term_mod_,
						dyn_c_term_mod_,
						ion_cutoff_percentage_,
						min_protein_mass_,
						max_protein_mass_,
						match_peak_tolerance_,
						stat_n_term_mod_,
						stat_c_term_mod_,
						stat_n_term_prot_mod_,
						stat_c_term_prot_mod_;
			
			int peptide_mass_units_,
					num_output_lines_,
					enzyme_number_,
					highest_enzyme_number_,
					max_num_differential_AA_per_mod_,
					max_num_mods_per_peptide_,
					nucleotide_reading_frame_,
					max_num_internal_cleavage_sites_,
					match_peak_count_,
					match_peak_allowed_error_;
			
			bool show_fragment_ions_,
					 print_duplicate_references_,
					 use_phospho_fragmentation_,
					 remove_precursor_peak_,
					 mass_type_parent_,
					 mass_type_fragment_,
					 normalize_xcorr_,
					 residues_in_upper_case_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_SEQUESTINFILE_H
