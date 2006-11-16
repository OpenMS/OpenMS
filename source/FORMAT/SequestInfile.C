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

#include <OpenMS/FORMAT/SequestInfile.h>

namespace OpenMS
{

	SequestInfile::SequestInfile():
		neutral_losses_for_ions_("0 1 1"),
		ion_series_weights_("0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0"),
		dyn_mods_("0 X"),
		peptide_mass_tolerance_(0),
		fragment_ion_tolerance_(0),
		dyn_n_term_mod_(0),
		dyn_c_term_mod_(0),
		ion_cutoff_percentage_(0),
		min_protein_mass_(0),
		max_protein_mass_(0),
		match_peak_tolerance_(0),
		stat_n_term_mod_(0),
		stat_c_term_mod_(0),
		stat_n_term_prot_mod_(0),
		stat_c_term_prot_mod_(0),
		peptide_mass_units_(0),
		enzyme_number_(0),
		max_num_differential_AA_per_mod_(0),
		max_num_mods_per_peptide_(0),
		nucleotide_reading_frame_(0),
		max_num_internal_cleavage_sites_(0),
		match_peak_count_(0),
		match_peak_allowed_error_(0),
		show_fragment_ions_(1),
		print_duplicate_references_(1),
		use_phospho_fragmentation_(0),
		remove_precursor_peak_(0),
		mass_type_parent_(0),
		mass_type_fragment_(0),
		normalize_xcorr_(0),
		residues_in_upper_case_(1)
	{
		for ( String::const_iterator aa_i = aas_single_letter_.begin(); aa_i != aas_single_letter_.end(); ++ aa_i )	stat_mods_[*aa_i] = 0.0;
		enzyme_info_ = getStandardEnzymeInfo();
	}
	
	SequestInfile::SequestInfile(const SequestInfile& sequest_infile)
	{
		stat_mods_ = sequest_infile.getStatMods();
		enzyme_info_ = sequest_infile.getEnzymeInfo(),
		database_ = sequest_infile.getDatabase(),
		snd_database_ = sequest_infile.getSndDatabase(),
		neutral_losses_for_ions_ = sequest_infile.getNeutralLossesForIons(),
		ion_series_weights_ = sequest_infile.getIonSeriesWeights(),
		dyn_mods_ = sequest_infile.getDynMods(),
		partial_sequence_ = sequest_infile.getPartialSequence(),
		sequence_header_filter_ = sequest_infile.getSequenceHeaderFilter();
		peptide_mass_tolerance_ = sequest_infile.getPeptideMassTolerance(),
		fragment_ion_tolerance_ = sequest_infile.getFragmentIonTolerance(),
		dyn_n_term_mod_ = sequest_infile.getDynNTermMod(),
		dyn_c_term_mod_ = sequest_infile.getDynCTermMod(),
		ion_cutoff_percentage_ = sequest_infile.getIonCutoffPercentage(),
		min_protein_mass_ = sequest_infile.getMinimumProteinMass(),
		max_protein_mass_ = sequest_infile.getMaximumProteinMass(),
		match_peak_tolerance_ = sequest_infile.getMatchPeakTolerance(),
		stat_n_term_mod_ = sequest_infile.getStatNTermMod(),
		stat_c_term_mod_ = sequest_infile.getStatCTermMod(),
		stat_n_term_prot_mod_ = sequest_infile.getStatNTermProtMod(),
		stat_c_term_prot_mod_ = sequest_infile.getStatCTermProtMod();
		peptide_mass_units_ = sequest_infile.getPeptideMassUnits(),
		num_output_lines_ = sequest_infile.getNumOutputLines(),
		enzyme_number_ = sequest_infile.getEnzymeNumber(),
		max_num_differential_AA_per_mod_ = sequest_infile.getMaxNumDifAAPerMod(),
		max_num_mods_per_peptide_ = sequest_infile.getMaxNumModsPerPeptide(),
		nucleotide_reading_frame_ = sequest_infile.getNucleotideReadingFrame(),
		max_num_internal_cleavage_sites_ = sequest_infile.getMaxNumInternalCleavageSites(),
		match_peak_count_ = sequest_infile.getMatchPeakCount(),
		match_peak_allowed_error_ = sequest_infile.getMatchPeakAllowedError();
		show_fragment_ions_ = sequest_infile.getShowFragmentIons(),
		print_duplicate_references_ = sequest_infile.getPrintDuplicateReferences(),
		use_phospho_fragmentation_ = sequest_infile.getUsePhosphoFragmentation(),
		remove_precursor_peak_ = sequest_infile.getRemovePrecursorPeak(),
		mass_type_parent_ = sequest_infile.getMassTypeParent(),
		mass_type_fragment_ = sequest_infile.getMassTypeFragment(),
		normalize_xcorr_ = sequest_infile.getNormalizeXcorr(),
		residues_in_upper_case_ = sequest_infile.getResiduesInUpperCase();
	}
	
	SequestInfile& SequestInfile::operator=(const SequestInfile& sequest_infile)
	{
		if ( this != &sequest_infile )
		{
			stat_mods_ = sequest_infile.getStatMods();
			enzyme_info_ = sequest_infile.getEnzymeInfo(),
			database_ = sequest_infile.getDatabase(),
			snd_database_ = sequest_infile.getSndDatabase(),
			neutral_losses_for_ions_ = sequest_infile.getNeutralLossesForIons(),
			ion_series_weights_ = sequest_infile.getIonSeriesWeights(),
			dyn_mods_ = sequest_infile.getDynMods(),
			partial_sequence_ = sequest_infile.getPartialSequence(),
			sequence_header_filter_ = sequest_infile.getSequenceHeaderFilter();
			peptide_mass_tolerance_ = sequest_infile.getPeptideMassTolerance(),
			fragment_ion_tolerance_ = sequest_infile.getFragmentIonTolerance(),
			dyn_n_term_mod_ = sequest_infile.getDynNTermMod(),
			dyn_c_term_mod_ = sequest_infile.getDynCTermMod(),
			ion_cutoff_percentage_ = sequest_infile.getIonCutoffPercentage(),
			min_protein_mass_ = sequest_infile.getMinimumProteinMass(),
			max_protein_mass_ = sequest_infile.getMaximumProteinMass(),
			match_peak_tolerance_ = sequest_infile.getMatchPeakTolerance(),
			stat_n_term_mod_ = sequest_infile.getStatNTermMod(),
			stat_c_term_mod_ = sequest_infile.getStatCTermMod(),
			stat_n_term_prot_mod_ = sequest_infile.getStatNTermProtMod(),
			stat_c_term_prot_mod_ = sequest_infile.getStatCTermProtMod();
			peptide_mass_units_ = sequest_infile.getPeptideMassUnits(),
			num_output_lines_ = sequest_infile.getNumOutputLines(),
			enzyme_number_ = sequest_infile.getEnzymeNumber(),
			max_num_differential_AA_per_mod_ = sequest_infile.getMaxNumDifAAPerMod(),
			max_num_mods_per_peptide_ = sequest_infile.getMaxNumModsPerPeptide(),
			nucleotide_reading_frame_ = sequest_infile.getNucleotideReadingFrame(),
			max_num_internal_cleavage_sites_ = sequest_infile.getMaxNumInternalCleavageSites(),
			match_peak_count_ = sequest_infile.getMatchPeakCount(),
			match_peak_allowed_error_ = sequest_infile.getMatchPeakAllowedError();
			show_fragment_ions_ = sequest_infile.getShowFragmentIons(),
			print_duplicate_references_ = sequest_infile.getPrintDuplicateReferences(),
			use_phospho_fragmentation_ = sequest_infile.getUsePhosphoFragmentation(),
			remove_precursor_peak_ = sequest_infile.getRemovePrecursorPeak(),
			mass_type_parent_ = sequest_infile.getMassTypeParent(),
			mass_type_fragment_ = sequest_infile.getMassTypeFragment(),
			normalize_xcorr_ = sequest_infile.getNormalizeXcorr(),
			residues_in_upper_case_ = sequest_infile.getResiduesInUpperCase();
		}
		return *this;
	}
	
	void
	SequestInfile::store(const String& filename) throw (Exception::UnableToCreateFile)
	{
		std::ofstream ofs(filename.c_str());
		if ( !ofs ) throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		
		// the header
		ofs << "[SEQUEST]" << std::endl;
		
		ofs << "database_name = " << database_ << std::endl;
		
		if ( !snd_database_.empty() ) ofs << "first_database_name = " << database_ << std::endl << "second_database_name = " << snd_database_ << std::endl;
		
		ofs << "peptide_mass_tolerance = " << peptide_mass_tolerance_ << std::endl;
		
		ofs << "peptide_mass_units = " << peptide_mass_units_ << "; 0=amu, 1=mmu, 2=ppm" << std::endl;
		
		ofs << "ion_series = " << neutral_losses_for_ions_ << " " << ion_series_weights_  << ";nABY ABCDVWXYZ" << std::endl;
		
		ofs << "fragment_ion_tolerance = " << fragment_ion_tolerance_ << std::endl;
		
		ofs << "num_output_lines = " << num_output_lines_ << std::endl;
		
		ofs << "num_results = " << num_output_lines_ << std::endl;
		
		ofs << "num_description_lines = 0" << std::endl;
		
		ofs << "show_fragment_ions = " << show_fragment_ions_ << std::endl;
		
		ofs << "print_duplicate_references = " << print_duplicate_references_ << std::endl;
		
		ofs << "enzyme_number = " << enzyme_number_ << std::endl;
		
		ofs << "diff_search_options = " << dyn_mods_ << std::endl;
		
		ofs << "term_diff_search_options = " << dyn_n_term_mod_ << " " << dyn_c_term_mod_ << std::endl;
		
		ofs << "use_phospho_fragmentation = " << use_phospho_fragmentation_ << std::endl;
		
		ofs << "remove_precursor_peak = " << remove_precursor_peak_ << std::endl;
		
		ofs << "ion_cutoff_percentage = " << ion_cutoff_percentage_ << std::endl;
		
		ofs << "protein_mass_filter = " << min_protein_mass_ << " " << max_protein_mass_ << std::endl;
		
		ofs << "max_num_differential_AA_per_mod = " << max_num_differential_AA_per_mod_ << std::endl;
		
		ofs << "max_num_differential_per_peptide = " << max_num_mods_per_peptide_ << std::endl;
		
		ofs << "nucleotide_reading_frame = " << nucleotide_reading_frame_ << ";\t\t0=protein db, 1-6, 7 = forward three, 8-reverse three, 9=all six" << std::endl;
		
		ofs << "mass_type_parent = " << mass_type_parent_ << ";\t\t0=average masses, 1=monoisotopic masses" << std::endl;
		
		ofs << "mass_type_fragment = " << mass_type_fragment_ << ";\t\t0=average masses, 1=monoisotopic masses" << std::endl;
		
		ofs << "normalize_xcorr = " << normalize_xcorr_ << std::endl;
		
		ofs << "max_num_internal_cleavage_sites = " << max_num_internal_cleavage_sites_ << std::endl;
		
		ofs << "create_output_files = 1" << std::endl;
		
		ofs << "partial_sequence = " << partial_sequence_ << std::endl;
		
		ofs << "sequence_header_filter = " << sequence_header_filter_ << std::endl;
		
		ofs << "match_peak_count = " << match_peak_count_ << ";\t\tnumber of auto-detected peaks to try matching (max 5)" << std::endl;
		
		ofs << "match_peak_allowed_error = " << match_peak_allowed_error_ << std::endl;
		
		ofs << "match_peak_tolerance = " << match_peak_tolerance_ << std::endl;
		
		ofs << "residues_in_upper_case = " << residues_in_upper_case_ << std::endl << std::endl << std::endl;
		
		ofs << "add_Nterm_peptide = " << stat_n_term_mod_ << std::endl;
		
		ofs << "add_Cterm_peptide = " << stat_c_term_mod_ << std::endl;
		
		ofs << "add_Nterm_protein = " << stat_n_term_prot_mod_ << std::endl;
		
		ofs << "add_Cterm_protein = " << stat_c_term_prot_mod_ << std::endl << std::endl;
		
		ofs << "add_G_Glycine = " << stat_mods_['G'] << ";\t\tadded to G - avg.  57.0519, mono.  57.02146" << std::endl;
		ofs << "add_A_Alanine = " << stat_mods_['A'] << ";\t\tadded to A - avg.  71.0788, mono.  71.03711" << std::endl;
		ofs << "add_S_Serine = " << stat_mods_['S'] << ";\t\tadded to S - avg.  87.0782, mono.  87.03203" << std::endl;
		ofs << "add_P_Proline = " << stat_mods_['P'] << ";\t\tadded to P - avg.  97.1167, mono.  97.05276" << std::endl;
		ofs << "add_V_Valine = " << stat_mods_['V'] << ";\t\tadded to V - avg.  99.1326, mono.  99.06841" << std::endl;
		ofs << "add_T_Threonine = " << stat_mods_['T'] << ";\t\tadded to T - avg. 101.1051, mono. 101.04768" << std::endl;
		ofs << "add_C_Cysteine = " << stat_mods_['C'] << ";\t\tadded to C - avg. 103.1388, mono. 103.00919" << std::endl;
		ofs << "add_L_Leucine = " << stat_mods_['L'] << ";\t\tadded to L - avg. 113.1594, mono. 113.08406" << std::endl;
		ofs << "add_I_Isoleucine = " << stat_mods_['I'] << ";\t\tadded to I - avg. 113.1594, mono. 113.08406" << std::endl;
		ofs << "add_X_LorI = " << stat_mods_['X'] << ";\t\tadded to X - avg. 113.1594, mono. 113.08406" << std::endl;
		ofs << "add_N_Asparagine = " << stat_mods_['N'] << ";\t\tadded to N - avg. 114.1038, mono. 114.04293" << std::endl;
		ofs << "add_O_Ornithine = " << stat_mods_['O'] << ";\t\tadded to O - avg. 114.1472, mono  114.07931" << std::endl;
		ofs << "add_B_avg_NandD = " << stat_mods_['B'] << ";\t\tadded to B - avg. 114.5962, mono. 114.53494" << std::endl;
		ofs << "add_D_Aspartic_Acid = " << stat_mods_['D'] << ";\t\tadded to D - avg. 115.0886, mono. 115.02694" << std::endl;
		ofs << "add_Q_Glutamine = " << stat_mods_['Q'] << ";\t\tadded to Q - avg. 128.1307, mono. 128.05858" << std::endl;
		ofs << "add_K_Lysine = " << stat_mods_['K'] << ";\t\tadded to K - avg. 128.1741, mono. 128.09496" << std::endl;
		ofs << "add_Z_avg_QandE = " << stat_mods_['Z'] << ";\t\tadded to Z - avg. 128.6231, mono. 128.55059" << std::endl;
		ofs << "add_E_Glutamic_Acid = " << stat_mods_['E'] << ";\t\tadded to E - avg. 129.1155, mono. 129.04259" << std::endl;
		ofs << "add_M_Methionine = " << stat_mods_['M'] << ";\t\tadded to M - avg. 131.1926, mono. 131.04049" << std::endl;
		ofs << "add_H_Histidine = " << stat_mods_['H'] << ";\t\tadded to H - avg. 137.1411, mono. 137.05891" << std::endl;
		ofs << "add_F_Phenylalanine = " << stat_mods_['F'] << ";\t\tadded to F - avg. 147.1766, mono. 147.06841" << std::endl;
		ofs << "add_R_Arginine = " << stat_mods_['R'] << ";\t\tadded to R - avg. 156.1875, mono. 156.10111" << std::endl;
		ofs << "add_Y_Tyrosine = " << stat_mods_['Y'] << ";\t\tadded to Y - avg. 163.1760, mono. 163.06333" << std::endl;
		ofs << "add_W_Tryptophan = " << stat_mods_['W'] << ";\t\tadded to W - avg. 186.2132, mono. 186.07931" << std::endl << std::endl;
		
		ofs << getStandardEnzymeInfo();
	}
	
	String SequestInfile::getStandardEnzymeInfo()
	{
		std::stringstream ss;
		ss << "[SEQUEST_ENZYME_INFO]" << std::endl;
		ss << "0.  No_Enzyme              0      -           -" << std::endl;
		ss << "1.  Trypsin_Strict         1      KR          -" << std::endl;
		ss << "2.  Trypsin                1      KRLNH       -" << std::endl;
		ss << "3.  Chymotrypsin           1      FWYL        -" << std::endl;
		ss << "4.  Chymotrypsin_WYF       1      FWY         -" << std::endl;
		ss << "5.  Clostripain            1      R           -" << std::endl;
		ss << "6.  Cyanogen_Bromide       1      M           -" << std::endl;
		ss << "7.  IodosoBenzoate         1      W           -" << std::endl;
		ss << "8.  Proline_Endopept       1      P           -" << std::endl;
		ss << "9.  GluC                   1      E           -" << std::endl;
		ss << "10. GluC_ED                1      ED          -" << std::endl;
		ss << "11. LysC                   1      K           -" << std::endl;
		ss << "12. AspN                   0      D           -" << std::endl;
		ss << "13. AspN_DE                0      DE          -" << std::endl;
		ss << "14. Elastase               1      ALIV        P" << std::endl;
		ss << "15. Elastase/Tryp/Chymo    1      ALIVKRWFY   P" << std::endl;
		ss << "16. Trypsin/Chymo          1      KRLFWYN     -" << std::endl;

		highest_enzyme_number_ = 16;
		
		return ss.str();
	}
	
	const String& SequestInfile::getEnzymeInfo() const {return enzyme_info_;}
	void SequestInfile::setEnzymeInfo(String& value){enzyme_info_ = value;}
	
	const String& SequestInfile::getDatabase() const {return database_;}
	void SequestInfile::setDatabase(String& value){database_ = value;}
	
	const String& SequestInfile::getSndDatabase() const {return snd_database_;}
	void SequestInfile::setSndDatabase(String& value){snd_database_ = value;}
	
	const String& SequestInfile::getNeutralLossesForIons() const {return neutral_losses_for_ions_;}
	void SequestInfile::setNeutralLossesForIons(String& neutral_losses_for_ions){neutral_losses_for_ions_ = neutral_losses_for_ions;}
	
	const String& SequestInfile::getIonSeriesWeights() const {return ion_series_weights_;}
	void SequestInfile::setIonSeriesWeights(String& ion_series_weights){ion_series_weights_ = ion_series_weights;}
	
	const String& SequestInfile::getDynMods() const {return dyn_mods_;}
	void SequestInfile::setDynMods(String& dyn_mods){dyn_mods_ = dyn_mods;}
	
	const String& SequestInfile::getPartialSequence() const {return partial_sequence_;}
	void SequestInfile::setPartialSequence(String& partial_sequence){partial_sequence_ = partial_sequence;}
	
	const String& SequestInfile::getSequenceHeaderFilter() const {return sequence_header_filter_;}
	void SequestInfile::setSequenceHeaderFilter(String& sequence_header_filter){sequence_header_filter_ = sequence_header_filter;}
	
	
	float SequestInfile::getPeptideMassTolerance() const {return peptide_mass_tolerance_;}
	void SequestInfile::setPeptideMassTolerance(float value){peptide_mass_tolerance_ = value;}
	
	float SequestInfile::getFragmentIonTolerance() const {return fragment_ion_tolerance_;}
	void SequestInfile::setFragmentIonTolerance(float value){fragment_ion_tolerance_ = value;}
	
	float SequestInfile::getMatchPeakTolerance() const {return match_peak_tolerance_;}
	void SequestInfile::setMatchPeakTolerance(float value){match_peak_tolerance_ = value;}
	
	float SequestInfile::getIonCutoffPercentage() const {return ion_cutoff_percentage_;}
	void SequestInfile::setIonCutoffPercentage(float value){ion_cutoff_percentage_ = value;}
	
	float SequestInfile::getMinimumProteinMass() const {return min_protein_mass_;}
	void SequestInfile::setMinimumProteinMass(float value){min_protein_mass_ = value;}
	
	float SequestInfile::getMaximumProteinMass() const {return max_protein_mass_;}
	void SequestInfile::setMaximumProteinMass(float value){max_protein_mass_ = value;}
	
	float SequestInfile::getDynNTermMod() const {return dyn_n_term_mod_;}
	void SequestInfile::setDynNTermMod(float mass){dyn_n_term_mod_ = mass;}
	
	float SequestInfile::getDynCTermMod() const {return dyn_c_term_mod_;}
	void SequestInfile::setDynCTermMod(float mass){dyn_c_term_mod_ = mass;}
	
	float SequestInfile::getStatNTermMod() const {return stat_n_term_mod_;}
	void SequestInfile::setStatNTermMod(float mass){stat_n_term_mod_ = mass;}
	
	float SequestInfile::getStatCTermMod() const {return stat_c_term_mod_;}
	void SequestInfile::setStatCTermMod(float mass){stat_c_term_mod_ = mass;}
	
	float SequestInfile::getStatNTermProtMod() const {return stat_n_term_prot_mod_;}
	void SequestInfile::setStatNTermProtMod(float mass){stat_n_term_prot_mod_ = mass;}
	
	float SequestInfile::getStatCTermProtMod() const {return stat_c_term_prot_mod_;}
	void SequestInfile::setStatCTermProtMod(float mass){stat_c_term_prot_mod_ = mass;}
	
	
	int SequestInfile::getPeptideMassUnits() const {return peptide_mass_units_;}
	void SequestInfile::setPeptideMassUnits(int value){peptide_mass_units_ = value;}
	
	int SequestInfile::getNumOutputLines() const {return num_output_lines_;}
	void SequestInfile::setNumOutputLines(int value){num_output_lines_ = value;}
	
	int SequestInfile::getEnzymeNumber() const {return enzyme_number_;}
	void SequestInfile::setEnzymeNumber(int value){enzyme_number_ = value;}
	
	int SequestInfile::getMaxNumDifAAPerMod() const {return max_num_differential_AA_per_mod_;}
	void SequestInfile::setMaxNumDifAAPerMod(int value){max_num_differential_AA_per_mod_ = value;}
	
	int SequestInfile::getMaxNumModsPerPeptide() const {return max_num_mods_per_peptide_;}
	void SequestInfile::setMaxNumModsPerPeptide(int value){max_num_mods_per_peptide_ = value;}
	
	int SequestInfile::getNucleotideReadingFrame() const {return nucleotide_reading_frame_;}
	void SequestInfile::setNucleotideReadingFrame(int value){nucleotide_reading_frame_ = value;}
	
	int SequestInfile::getMaxNumInternalCleavageSites() const {return max_num_internal_cleavage_sites_;}
	void SequestInfile::setMaxNumInternalCleavageSites(int value){max_num_internal_cleavage_sites_ = value;}
	
	int SequestInfile::getMatchPeakCount() const {return match_peak_count_;}
	void SequestInfile::setMatchPeakCount(int value){match_peak_count_ = value;}
	
	int SequestInfile::getMatchPeakAllowedError() const {return match_peak_allowed_error_;}
	void SequestInfile::setMatchPeakAllowedError(int value){match_peak_allowed_error_ = value;}
	
	
	bool SequestInfile::getShowFragmentIons() const {return show_fragment_ions_;}
	void SequestInfile::setShowFragmentIons(bool value){show_fragment_ions_ = value;}
	
	bool SequestInfile::getPrintDuplicateReferences() const {return print_duplicate_references_;}
	void SequestInfile::setPrintDuplicateReferences(bool value){print_duplicate_references_ = value;}
	
	bool SequestInfile::getUsePhosphoFragmentation() const {return use_phospho_fragmentation_;}
	void SequestInfile::setUsePhosphoFragmentation(bool value){use_phospho_fragmentation_ = value;}
	
	bool SequestInfile::getRemovePrecursorPeak() const {return remove_precursor_peak_;}
	void SequestInfile::setRemovePrecursorPeak(bool value){remove_precursor_peak_ = value;}
	
	bool SequestInfile::getMassTypeParent() const {return mass_type_parent_;}
	void SequestInfile::setMassTypeParent(bool value){mass_type_parent_ = value;}
	
	bool SequestInfile::getMassTypeFragment() const {return mass_type_fragment_;}
	void SequestInfile::setMassTypeFragment(bool value){mass_type_fragment_ = value;}
	
	bool SequestInfile::getNormalizeXcorr() const {return normalize_xcorr_;}
	void SequestInfile::setNormalizeXcorr(bool value){normalize_xcorr_ = value;}
	
	bool SequestInfile::getResiduesInUpperCase() const {return residues_in_upper_case_;}
	void SequestInfile::setResiduesInUpperCase(bool value){residues_in_upper_case_ = value;}
	
	
	void SequestInfile::addEnzymeInfo(std::vector< String >& enzyme_info)
	{
		enzyme_number_ = ++highest_enzyme_number_;
		
		enzyme_info_.append(String(enzyme_number_) + ". " + enzyme_info[0].substr(0, 22) + String(std::max(22 - (int) enzyme_info[0].length(), 0), ' ') + " " + enzyme_info[1] + String(6, ' ') + enzyme_info[2] + " " + enzyme_info[3] + "\n");
	}
	
	const std::map< char, float >& SequestInfile::getStatMods() const {return stat_mods_;}
	void SequestInfile::setStatMods(std::map< char, float >& stat_mods){stat_mods_ = stat_mods;}
	
	void SequestInfile::setStatMod(char amino_acid, float mass)
	{
		if ( amino_acid > 96 ) amino_acid -= 32;
		stat_mods_[amino_acid] = mass;
	}
	
	const String SequestInfile::aas_single_letter_ = "GASPVTCLIXNOBDQKZEMHFRYW";
}
