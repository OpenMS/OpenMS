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

#include <OpenMS/FORMAT/SequestInfile.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{

	SequestInfile::SequestInfile():
		neutral_losses_for_ions_("0 1 1"),
		ion_series_weights_("0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0"),
		dyn_mods_("0 X"),
		protein_mass_filter_("0 0"),
		peptide_mass_tolerance_(0),
		fragment_ion_tolerance_(0),
		match_peak_tolerance_(0),
		ion_cutoff_percentage_(0),
		dyn_n_term_mod_(0),
		dyn_c_term_mod_(0),
		stat_n_term_mod_(0),
		stat_c_term_mod_(0),
		stat_n_term_prot_mod_(0),
		stat_c_term_prot_mod_(0),
		peptide_mass_unit_(0),
		enzyme_number_(0),
		max_AA_per_mod_per_peptide_(0),
		max_mods_per_peptide_(0),
		nucleotide_reading_frame_(0),
		max_internal_cleavage_sites_(0),
		match_peak_count_(0),
		match_peak_allowed_error_(0),
		show_fragment_ions_(1),
		print_duplicate_references_(1),
//		use_phospho_fragmentation_(0),
		remove_precursor_near_peaks_(0),
		mass_type_parent_(0),
		mass_type_fragment_(0),
		normalize_xcorr_(0),
		residues_in_upper_case_(1)
	{
		for ( String::const_iterator aa_i = aas_single_letter_.begin(); aa_i != aas_single_letter_.end(); ++ aa_i )	stat_mods_[*aa_i] = 0.0;
		setStandardEnzymeInfo();
	}
	
	SequestInfile::SequestInfile(const SequestInfile& sequest_infile)
	{
		stat_mods_ = sequest_infile.getStatMods();
		enzyme_info_ = sequest_infile.getEnzymeInfo(),
		database_ = sequest_infile.getDatabase(),
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
		protein_mass_filter_ = sequest_infile.getProteinMassFilter(),
		match_peak_tolerance_ = sequest_infile.getMatchPeakTolerance(),
		stat_n_term_mod_ = sequest_infile.getStatNTermMod(),
		stat_c_term_mod_ = sequest_infile.getStatCTermMod(),
		stat_n_term_prot_mod_ = sequest_infile.getStatNTermProtMod(),
		stat_c_term_prot_mod_ = sequest_infile.getStatCTermProtMod();
		peptide_mass_unit_ = sequest_infile.getPeptideMassUnit(),
		output_lines_ = sequest_infile.getOutputLines(),
		enzyme_number_ = sequest_infile.getEnzymeNumber(),
		max_AA_per_mod_per_peptide_ = sequest_infile.getMaxAAPerModPerPeptide(),
		max_mods_per_peptide_ = sequest_infile.getMaxModsPerPeptide(),
		nucleotide_reading_frame_ = sequest_infile.getNucleotideReadingFrame(),
		max_internal_cleavage_sites_ = sequest_infile.getMaxInternalCleavageSites(),
		match_peak_count_ = sequest_infile.getMatchPeakCount(),
		match_peak_allowed_error_ = sequest_infile.getMatchPeakAllowedError();
		show_fragment_ions_ = sequest_infile.getShowFragmentIons(),
		print_duplicate_references_ = sequest_infile.getPrintDuplicateReferences(),
//		use_phospho_fragmentation_ = sequest_infile.getUsePhosphoFragmentation(),
		remove_precursor_near_peaks_ = sequest_infile.getRemovePrecursorNearPeaks(),
		mass_type_parent_ = sequest_infile.getMassTypeParent(),
		mass_type_fragment_ = sequest_infile.getMassTypeFragment(),
		normalize_xcorr_ = sequest_infile.getNormalizeXcorr(),
		residues_in_upper_case_ = sequest_infile.getResiduesInUpperCase();
	}

	SequestInfile::~SequestInfile()
	{
		stat_mods_.clear();
	}
	
	SequestInfile& SequestInfile::operator=(const SequestInfile& sequest_infile)
	{
		if ( this != &sequest_infile )
		{
			stat_mods_ = sequest_infile.getStatMods();
			enzyme_info_ = sequest_infile.getEnzymeInfo(),
			database_ = sequest_infile.getDatabase(),
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
			protein_mass_filter_ = sequest_infile.getProteinMassFilter(),
			match_peak_tolerance_ = sequest_infile.getMatchPeakTolerance(),
			stat_n_term_mod_ = sequest_infile.getStatNTermMod(),
			stat_c_term_mod_ = sequest_infile.getStatCTermMod(),
			stat_n_term_prot_mod_ = sequest_infile.getStatNTermProtMod(),
			stat_c_term_prot_mod_ = sequest_infile.getStatCTermProtMod();
			peptide_mass_unit_ = sequest_infile.getPeptideMassUnit(),
			output_lines_ = sequest_infile.getOutputLines(),
			enzyme_number_ = sequest_infile.getEnzymeNumber(),
			max_AA_per_mod_per_peptide_ = sequest_infile.getMaxAAPerModPerPeptide(),
			max_mods_per_peptide_ = sequest_infile.getMaxModsPerPeptide(),
			nucleotide_reading_frame_ = sequest_infile.getNucleotideReadingFrame(),
			max_internal_cleavage_sites_ = sequest_infile.getMaxInternalCleavageSites(),
			match_peak_count_ = sequest_infile.getMatchPeakCount(),
			match_peak_allowed_error_ = sequest_infile.getMatchPeakAllowedError();
			show_fragment_ions_ = sequest_infile.getShowFragmentIons(),
			print_duplicate_references_ = sequest_infile.getPrintDuplicateReferences(),
//			use_phospho_fragmentation_ = sequest_infile.getUsePhosphoFragmentation(),
			remove_precursor_near_peaks_ = sequest_infile.getRemovePrecursorNearPeaks(),
			mass_type_parent_ = sequest_infile.getMassTypeParent(),
			mass_type_fragment_ = sequest_infile.getMassTypeFragment(),
			normalize_xcorr_ = sequest_infile.getNormalizeXcorr(),
			residues_in_upper_case_ = sequest_infile.getResiduesInUpperCase();
		}
		return *this;
	}
	
	void
	SequestInfile::store(
		const String& filename)
	throw (
		Exception::UnableToCreateFile)
	{
		ofstream ofs(filename.c_str());
		if ( !ofs ) throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		
		// the header
		ofs << "[SEQUEST]" << endl;
		
		ofs << "database_name = " << database_ << endl;
		
		ofs << "peptide_mass_tolerance = " << peptide_mass_tolerance_ << endl;
		
		ofs << "peptide_mass_units = " << peptide_mass_unit_ << "; 0=amu, 1=mmu, 2=ppm" << endl;
		
		ofs << "ion_series = " << neutral_losses_for_ions_ << " " << ion_series_weights_  << ";nABY ABCDVWXYZ" << endl;
		
		ofs << "fragment_ion_tolerance = " << fragment_ion_tolerance_ << endl;
		
		ofs << "num_output_lines = " << output_lines_ << endl;
		
		ofs << "num_results = " << output_lines_ << endl;
		
		ofs << "num_description_lines = 0" << endl;
		
		ofs << "show_fragment_ions = " << show_fragment_ions_ << endl;
		
		ofs << "print_duplicate_references = " << print_duplicate_references_ << endl;
		
		ofs << "enzyme_number = " << enzyme_number_ << endl;
		
		ofs << "diff_search_options = " << dyn_mods_ << endl;
		
		ofs << "term_diff_search_options = " << dyn_n_term_mod_ << " " << dyn_c_term_mod_ << endl;
		
//		ofs << "use_phospho_fragmentation = " << use_phospho_fragmentation_ << endl;
		
		ofs << "remove_precursor_peak = " << remove_precursor_near_peaks_ << endl;
		
		ofs << "ion_cutoff_percentage = " << ion_cutoff_percentage_ << endl;
		
		ofs << "protein_mass_filter = " << protein_mass_filter_ << endl;
		
		ofs << "max_differential_AA_per_mod = " << max_AA_per_mod_per_peptide_ << endl;
		
		ofs << "max_differential_per_peptide = " << max_mods_per_peptide_ << endl;
		
		ofs << "nucleotide_reading_frame = " << nucleotide_reading_frame_ << "; 0=protein db, 1-6, 7 = forward three, 8-reverse three, 9=all six" << endl;
		
		ofs << "mass_type_parent = " << mass_type_parent_ << "; 0=average masses, 1=monoisotopic masses" << endl;
		
		ofs << "mass_type_fragment = " << mass_type_fragment_ << "; 0=average masses, 1=monoisotopic masses" << endl;
		
		ofs << "normalize_xcorr = " << normalize_xcorr_ << endl;
		
		ofs << "max_internal_cleavage_sites = " << max_internal_cleavage_sites_ << endl;
		
		ofs << "create_output_files = 1" << endl;
		
		ofs << "partial_sequence = " << partial_sequence_ << endl;
		
		ofs << "sequence_header_filter = " << sequence_header_filter_ << endl;
		
		ofs << "match_peak_count = " << match_peak_count_ << "; number of auto-detected peaks to try matching (max 5)" << endl;
		
		ofs << "match_peak_allowed_error = " << match_peak_allowed_error_ << endl;
		
		ofs << "match_peak_tolerance = " << match_peak_tolerance_ << endl;
		
		ofs << "residues_in_upper_case = " << residues_in_upper_case_ << endl << endl << endl;
		
		ofs << "add_Nterm_peptide = " << stat_n_term_mod_ << endl;
		
		ofs << "add_Cterm_peptide = " << stat_c_term_mod_ << endl;
		
		ofs << "add_Nterm_protein = " << stat_n_term_prot_mod_ << endl;
		
		ofs << "add_Cterm_protein = " << stat_c_term_prot_mod_ << endl << endl;
		
		ofs << "add_G_Glycine = " << stat_mods_['G'] << "; added to G - avg.  57.0519, mono.  57.02146" << endl;
		ofs << "add_A_Alanine = " << stat_mods_['A'] << "; added to A - avg.  71.0788, mono.  71.03711" << endl;
		ofs << "add_S_Serine = " << stat_mods_['S'] << "; added to S - avg.  87.0782, mono.  87.03203" << endl;
		ofs << "add_P_Proline = " << stat_mods_['P'] << "; added to P - avg.  97.1167, mono.  97.05276" << endl;
		ofs << "add_V_Valine = " << stat_mods_['V'] << "; added to V - avg.  99.1326, mono.  99.06841" << endl;
		ofs << "add_T_Threonine = " << stat_mods_['T'] << "; added to T - avg. 101.1051, mono. 101.04768" << endl;
		ofs << "add_C_Cysteine = " << stat_mods_['C'] << "; added to C - avg. 103.1388, mono. 103.00919" << endl;
		ofs << "add_L_Leucine = " << stat_mods_['L'] << "; added to L - avg. 113.1594, mono. 113.08406" << endl;
		ofs << "add_I_Isoleucine = " << stat_mods_['I'] << "; added to I - avg. 113.1594, mono. 113.08406" << endl;
		ofs << "add_X_LorI = " << stat_mods_['X'] << "; added to X - avg. 113.1594, mono. 113.08406" << endl;
		ofs << "add_N_Asparagine = " << stat_mods_['N'] << "; added to N - avg. 114.1038, mono. 114.04293" << endl;
		ofs << "add_O_Ornithine = " << stat_mods_['O'] << "; added to O - avg. 114.1472, mono  114.07931" << endl;
		ofs << "add_B_avg_NandD = " << stat_mods_['B'] << "; added to B - avg. 114.5962, mono. 114.53494" << endl;
		ofs << "add_D_Aspartic_Acid = " << stat_mods_['D'] << "; added to D - avg. 115.0886, mono. 115.02694" << endl;
		ofs << "add_Q_Glutamine = " << stat_mods_['Q'] << "; added to Q - avg. 128.1307, mono. 128.05858" << endl;
		ofs << "add_K_Lysine = " << stat_mods_['K'] << "; added to K - avg. 128.1741, mono. 128.09496" << endl;
		ofs << "add_Z_avg_QandE = " << stat_mods_['Z'] << "; added to Z - avg. 128.6231, mono. 128.55059" << endl;
		ofs << "add_E_Glutamic_Acid = " << stat_mods_['E'] << "; added to E - avg. 129.1155, mono. 129.04259" << endl;
		ofs << "add_M_Methionine = " << stat_mods_['M'] << "; added to M - avg. 131.1926, mono. 131.04049" << endl;
		ofs << "add_H_Histidine = " << stat_mods_['H'] << "; added to H - avg. 137.1411, mono. 137.05891" << endl;
		ofs << "add_F_Phenylalanine = " << stat_mods_['F'] << "; added to F - avg. 147.1766, mono. 147.06841" << endl;
		ofs << "add_R_Arginine = " << stat_mods_['R'] << "; added to R - avg. 156.1875, mono. 156.10111" << endl;
		ofs << "add_Y_Tyrosine = " << stat_mods_['Y'] << "; added to Y - avg. 163.1760, mono. 163.06333" << endl;
		ofs << "add_W_Tryptophan = " << stat_mods_['W'] << "; added to W - avg. 186.2132, mono. 186.07931" << endl << endl;
		
		ofs << getEnzymeInfoAsString();
	}

	const map< String, vector< String > >& SequestInfile::getEnzymeInfo() const {return enzyme_info_;}
	
	const String SequestInfile::getEnzymeInfoAsString() const
	{
		stringstream ss;
		UnsignedInt i = 0;
		String::size_type max_name_length = 0;
		String::size_type max_cut_before_length = 0;
		String::size_type max_doesnt_cut_after_length = 0;
		ss << "[SEQUEST_ENZYME_INFO]" << endl;
		for ( map< String, vector< String > >::const_iterator einfo_i = enzyme_info_.begin(); einfo_i != enzyme_info_.end(); ++einfo_i )
		{
			max_name_length = max(max_name_length, einfo_i->first.length());
			max_cut_before_length = max(max_cut_before_length, einfo_i->second[1].length());
			max_doesnt_cut_after_length = max(max_doesnt_cut_after_length, einfo_i->second[2].length());
		}
		for ( map< String, vector< String > >::const_iterator einfo_i = enzyme_info_.begin(); einfo_i != enzyme_info_.end(); ++einfo_i, ++i )
		{
			ss << i << ".  " << einfo_i->first << String(max_name_length + 5 - einfo_i->first.length(), ' ') << einfo_i->second[0] << "     " << einfo_i->second[1] << String(max_cut_before_length + 5 - einfo_i->second[1].length(), ' ') << einfo_i->second[2] << endl;
		}
		return String(ss.str());
	}
	
	void SequestInfile::addEnzymeInfo(vector< String >& enzyme_info)
	{
		// remove duplicates from the concerned amino acids
		set< char > aas;
		for ( String::const_iterator s_i = enzyme_info[2].begin(); s_i != enzyme_info[2].end(); ++s_i )
		{
			aas.insert(*s_i);
		}
		if ( enzyme_info[2].length() != aas.size() )
		{
			enzyme_info[2].clear();
			enzyme_info[2].reserve(aas.size());
			for ( set< char >::const_iterator aa_i = aas.begin(); aa_i != aas.end(); ++aa_i )
			{
				enzyme_info[2].append(1, *aa_i);
			}
		}
		
		String enzyme_name = enzyme_info[0];
		enzyme_info.erase(enzyme_info.begin());
		enzyme_info_[enzyme_name] = enzyme_info;
		enzyme_number_ = 0;
		for ( std::map< String, std::vector< String > >::const_iterator einfo_i = enzyme_info_.begin(); einfo_i != enzyme_info_.end(); ++einfo_i, ++enzyme_number_ )
		{
			if ( einfo_i->first == enzyme_name ) break;
		}
	}
	
	const String& SequestInfile::getDatabase() const {return database_;}
	void SequestInfile::setDatabase(const String& value){database_ = value;}
	
	const String& SequestInfile::getNeutralLossesForIons() const {return neutral_losses_for_ions_;}
	void SequestInfile::setNeutralLossesForIons(const String& neutral_losses_for_ions){neutral_losses_for_ions_ = neutral_losses_for_ions;}
	
	const String& SequestInfile::getIonSeriesWeights() const {return ion_series_weights_;}
	void SequestInfile::setIonSeriesWeights(const String& ion_series_weights){ion_series_weights_ = ion_series_weights;}
	
	const String& SequestInfile::getDynMods() const {return dyn_mods_;}
	void SequestInfile::setDynMods(const String& dyn_mods){dyn_mods_ = dyn_mods;}
	
	const String& SequestInfile::getPartialSequence() const {return partial_sequence_;}
	void SequestInfile::setPartialSequence(const String& partial_sequence){partial_sequence_ = partial_sequence;}
	
	const String& SequestInfile::getSequenceHeaderFilter() const {return sequence_header_filter_;}
	void SequestInfile::setSequenceHeaderFilter(const String& sequence_header_filter){sequence_header_filter_ = sequence_header_filter;}

	const String& SequestInfile::getProteinMassFilter() const {return protein_mass_filter_;}
	void SequestInfile::setProteinMassFilter(const String& protein_mass_filter){protein_mass_filter_ = protein_mass_filter;}
	
	
	Real SequestInfile::getPeptideMassTolerance() const {return peptide_mass_tolerance_;}
	void SequestInfile::setPeptideMassTolerance(Real value){peptide_mass_tolerance_ = value;}
	
	Real SequestInfile::getFragmentIonTolerance() const {return fragment_ion_tolerance_;}
	void SequestInfile::setFragmentIonTolerance(Real value){fragment_ion_tolerance_ = value;}
	
	Real SequestInfile::getMatchPeakTolerance() const {return match_peak_tolerance_;}
	void SequestInfile::setMatchPeakTolerance(Real value){match_peak_tolerance_ = value;}
	
	Real SequestInfile::getIonCutoffPercentage() const {return ion_cutoff_percentage_;}
	void SequestInfile::setIonCutoffPercentage(Real value){ion_cutoff_percentage_ = value;}
	
	Real SequestInfile::getDynNTermMod() const {return dyn_n_term_mod_;}
	void SequestInfile::setDynNTermMod(Real mass){dyn_n_term_mod_ = mass;}
	
	Real SequestInfile::getDynCTermMod() const {return dyn_c_term_mod_;}
	void SequestInfile::setDynCTermMod(Real mass){dyn_c_term_mod_ = mass;}
	
	Real SequestInfile::getStatNTermMod() const {return stat_n_term_mod_;}
	void SequestInfile::setStatNTermMod(Real mass){stat_n_term_mod_ = mass;}
	
	Real SequestInfile::getStatCTermMod() const {return stat_c_term_mod_;}
	void SequestInfile::setStatCTermMod(Real mass){stat_c_term_mod_ = mass;}
	
	Real SequestInfile::getStatNTermProtMod() const {return stat_n_term_prot_mod_;}
	void SequestInfile::setStatNTermProtMod(Real mass){stat_n_term_prot_mod_ = mass;}
	
	Real SequestInfile::getStatCTermProtMod() const {return stat_c_term_prot_mod_;}
	void SequestInfile::setStatCTermProtMod(Real mass){stat_c_term_prot_mod_ = mass;}
	
	
	SignedInt SequestInfile::getPeptideMassUnit() const {return peptide_mass_unit_;}
	void SequestInfile::setPeptideMassUnit(SignedInt value){peptide_mass_unit_ = value;}
	
	SignedInt SequestInfile::getOutputLines() const {return output_lines_;}
	void SequestInfile::setOutputLines(SignedInt value){output_lines_ = value;}
	
	SignedInt SequestInfile::getEnzymeNumber() const {return enzyme_number_;}
	SignedInt SequestInfile::setEnzyme(String enzyme_name)
	{
		enzyme_number_ = 0;
		std::map< String, std::vector< String > >::const_iterator einfo_i;
		for ( einfo_i = enzyme_info_.begin(); einfo_i != enzyme_info_.end(); ++einfo_i, ++enzyme_number_ )
		{
			if ( einfo_i->first == enzyme_name ) break;
		}
		return ( einfo_i == enzyme_info_.end() ) ? enzyme_info_.size() : 0;
	}
	
	SignedInt SequestInfile::getMaxAAPerModPerPeptide() const {return max_AA_per_mod_per_peptide_;}
	void SequestInfile::setMaxAAPerModPerPeptide(SignedInt value){max_AA_per_mod_per_peptide_ = value;}
	
	SignedInt SequestInfile::getMaxModsPerPeptide() const {return max_mods_per_peptide_;}
	void SequestInfile::setMaxModsPerPeptide(SignedInt value){max_mods_per_peptide_ = value;}
	
	SignedInt SequestInfile::getNucleotideReadingFrame() const {return nucleotide_reading_frame_;}
	void SequestInfile::setNucleotideReadingFrame(SignedInt value){nucleotide_reading_frame_ = value;}
	
	SignedInt SequestInfile::getMaxInternalCleavageSites() const {return max_internal_cleavage_sites_;}
	void SequestInfile::setMaxInternalCleavageSites(SignedInt value){max_internal_cleavage_sites_ = value;}
	
	SignedInt SequestInfile::getMatchPeakCount() const {return match_peak_count_;}
	void SequestInfile::setMatchPeakCount(SignedInt value){match_peak_count_ = value;}
	
	SignedInt SequestInfile::getMatchPeakAllowedError() const {return match_peak_allowed_error_;}
	void SequestInfile::setMatchPeakAllowedError(SignedInt value){match_peak_allowed_error_ = value;}
	
	
	bool SequestInfile::getShowFragmentIons() const {return show_fragment_ions_;}
	void SequestInfile::setShowFragmentIons(bool value){show_fragment_ions_ = value;}
	
	bool SequestInfile::getPrintDuplicateReferences() const {return print_duplicate_references_;}
	void SequestInfile::setPrintDuplicateReferences(bool value){print_duplicate_references_ = value;}
	
// 	bool SequestInfile::getUsePhosphoFragmentation() const {return use_phospho_fragmentation_;}
// 	void SequestInfile::setUsePhosphoFragmentation(bool value){use_phospho_fragmentation_ = value;}
	
	bool SequestInfile::getRemovePrecursorNearPeaks() const {return remove_precursor_near_peaks_;}
	void SequestInfile::setRemovePrecursorNearPeaks(bool value){remove_precursor_near_peaks_ = value;}
	
	bool SequestInfile::getMassTypeParent() const {return mass_type_parent_;}
	void SequestInfile::setMassTypeParent(bool value){mass_type_parent_ = value;}
	
	bool SequestInfile::getMassTypeFragment() const {return mass_type_fragment_;}
	void SequestInfile::setMassTypeFragment(bool value){mass_type_fragment_ = value;}
	
	bool SequestInfile::getNormalizeXcorr() const {return normalize_xcorr_;}
	void SequestInfile::setNormalizeXcorr(bool value){normalize_xcorr_ = value;}
	
	bool SequestInfile::getResiduesInUpperCase() const {return residues_in_upper_case_;}
	void SequestInfile::setResiduesInUpperCase(bool value){residues_in_upper_case_ = value;}
	
	
	const map< char, Real >& SequestInfile::getStatMods() const {return stat_mods_;}
	
	char SequestInfile::setStatMod(String amino_acids, Real mass)
	{
		amino_acids.toUpper();
		for ( String::const_iterator s_i = amino_acids.begin(); s_i != amino_acids.end(); ++s_i )
		{
			if ( aas_single_letter_.find(*s_i) == string::npos ) return *s_i; // if an unknown amino acid is used, report it
			stat_mods_[*s_i] = mass;
		}
		return 0;
	}
	
	void SequestInfile::setStandardEnzymeInfo()
	{
		vector< String > info;
		//		 cuts n to c?							 cuts before							doesn't cut after
		info.push_back("0");info.push_back("-");info.push_back("-");enzyme_info_["No_Enzyme"] = info; info.clear();
		info.push_back("1");info.push_back("KR");info.push_back("-");enzyme_info_["Trypsin_Strict"] = info; info.clear();
		info.push_back("1");info.push_back("KRLNH");info.push_back("-");enzyme_info_["Trypsin"] = info; info.clear();
		info.push_back("1");info.push_back("FWYL");info.push_back("-");enzyme_info_["Chymotrypsin"] = info; info.clear();
		info.push_back("1");info.push_back("FWY");info.push_back("-");enzyme_info_["Chymotrypsin_WYF"] = info; info.clear();
		info.push_back("1");info.push_back("R");info.push_back("-");enzyme_info_["Clostripain"] = info; info.clear();
		info.push_back("1");info.push_back("M");info.push_back("-");enzyme_info_["Cyanogen_Bromide"] = info; info.clear();
		info.push_back("1");info.push_back("W");info.push_back("-");enzyme_info_["IodosoBenzoate"] = info; info.clear();
		info.push_back("1");info.push_back("P");info.push_back("-");enzyme_info_["Proline_Endopept"] = info; info.clear();
		info.push_back("1");info.push_back("E");info.push_back("-");enzyme_info_["GluC"] = info; info.clear();
		info.push_back("1");info.push_back("ED");info.push_back("-");enzyme_info_["GluC_ED"] = info; info.clear();
		info.push_back("1");info.push_back("K");info.push_back("-");enzyme_info_["LysC"] = info; info.clear();
		info.push_back("0");info.push_back("D");info.push_back("-");enzyme_info_["AspN"] = info; info.clear();
		info.push_back("0");info.push_back("DE");info.push_back("-");enzyme_info_["AspN_DE"] = info; info.clear();
		info.push_back("1");info.push_back("ALIV");info.push_back("P");enzyme_info_["Elastase"] = info; info.clear();
		info.push_back("1");info.push_back("ALIVKRWFY");info.push_back("P");enzyme_info_["Elastase/Tryp/Chymo"] = info; info.clear();
		info.push_back("1");info.push_back("KRLFWYN");info.push_back("-");enzyme_info_["Trypsin/Chymo"] = info; info.clear();
	}
	
	const String SequestInfile::aas_single_letter_ = "GASPVTCLIXNOBDQKZEMHFRYW";
}
