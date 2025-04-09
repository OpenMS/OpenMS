// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SequestInfile.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>

#include <fstream>
#include <sstream>

using namespace std;

namespace OpenMS
{

  SequestInfile::SequestInfile() :
    neutral_losses_for_ions_("0 1 1"),
    ion_series_weights_("0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0"),
    protein_mass_filter_("0 0"),
    precursor_mass_tolerance_(0),
    peak_mass_tolerance_(0),
    match_peak_tolerance_(0),
    ion_cutoff_percentage_(0),
    peptide_mass_unit_(0),
    output_lines_(0),
    enzyme_number_(0),
    max_AA_per_mod_per_peptide_(0),
    max_mods_per_peptide_(0),
    nucleotide_reading_frame_(0),
    max_internal_cleavage_sites_(0),
    match_peak_count_(0),
    match_peak_allowed_error_(0),
    show_fragment_ions_(true),
    print_duplicate_references_(true),
    remove_precursor_near_peaks_(false),
    mass_type_parent_(false),
    mass_type_fragment_(false),
    normalize_xcorr_(false),
    residues_in_upper_case_(true)
  {
    setStandardEnzymeInfo_();
  }

  SequestInfile::SequestInfile(const SequestInfile & sequest_infile)
  {
    enzyme_info_ = sequest_infile.getEnzymeInfo_(),
    database_ = sequest_infile.getDatabase(),
    neutral_losses_for_ions_ = sequest_infile.getNeutralLossesForIons(),
    ion_series_weights_ = sequest_infile.getIonSeriesWeights(),
    partial_sequence_ = sequest_infile.getPartialSequence(),
    sequence_header_filter_ = sequest_infile.getSequenceHeaderFilter();
    precursor_mass_tolerance_ = sequest_infile.getPrecursorMassTolerance(),
    peak_mass_tolerance_ = sequest_infile.getPeakMassTolerance(),
    ion_cutoff_percentage_ = sequest_infile.getIonCutoffPercentage(),
    protein_mass_filter_ = sequest_infile.getProteinMassFilter(),
    match_peak_tolerance_ = sequest_infile.getMatchPeakTolerance(),
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
    remove_precursor_near_peaks_ = sequest_infile.getRemovePrecursorNearPeaks(),
    mass_type_parent_ = sequest_infile.getMassTypeParent(),
    mass_type_fragment_ = sequest_infile.getMassTypeFragment(),
    normalize_xcorr_ = sequest_infile.getNormalizeXcorr(),
    residues_in_upper_case_ = sequest_infile.getResiduesInUpperCase();
    PTMname_residues_mass_type_ = sequest_infile.getModifications();
  }

  SequestInfile::~SequestInfile()
  {
    PTMname_residues_mass_type_.clear();
  }

  SequestInfile & SequestInfile::operator=(const SequestInfile & sequest_infile)
  {
    if (this == &sequest_infile)
      return *this;

    enzyme_info_ = sequest_infile.getEnzymeInfo_();
    database_ = sequest_infile.getDatabase();
    neutral_losses_for_ions_ = sequest_infile.getNeutralLossesForIons();
    ion_series_weights_ = sequest_infile.getIonSeriesWeights();
    partial_sequence_ = sequest_infile.getPartialSequence();
    sequence_header_filter_ = sequest_infile.getSequenceHeaderFilter();
    precursor_mass_tolerance_ = sequest_infile.getPrecursorMassTolerance();
    peak_mass_tolerance_ = sequest_infile.getPeakMassTolerance();
    ion_cutoff_percentage_ = sequest_infile.getIonCutoffPercentage();
    protein_mass_filter_ = sequest_infile.getProteinMassFilter();
    match_peak_tolerance_ = sequest_infile.getMatchPeakTolerance();
    peptide_mass_unit_ = sequest_infile.getPeptideMassUnit();
    output_lines_ = sequest_infile.getOutputLines();
    enzyme_number_ = sequest_infile.getEnzymeNumber();
    max_AA_per_mod_per_peptide_ = sequest_infile.getMaxAAPerModPerPeptide();
    max_mods_per_peptide_ = sequest_infile.getMaxModsPerPeptide();
    nucleotide_reading_frame_ = sequest_infile.getNucleotideReadingFrame();
    max_internal_cleavage_sites_ = sequest_infile.getMaxInternalCleavageSites();
    match_peak_count_ = sequest_infile.getMatchPeakCount();
    match_peak_allowed_error_ = sequest_infile.getMatchPeakAllowedError();
    show_fragment_ions_ = sequest_infile.getShowFragmentIons();
    print_duplicate_references_ = sequest_infile.getPrintDuplicateReferences();
    remove_precursor_near_peaks_ = sequest_infile.getRemovePrecursorNearPeaks();
    mass_type_parent_ = sequest_infile.getMassTypeParent();
    mass_type_fragment_ = sequest_infile.getMassTypeFragment();
    normalize_xcorr_ = sequest_infile.getNormalizeXcorr();
    residues_in_upper_case_ = sequest_infile.getResiduesInUpperCase();
    PTMname_residues_mass_type_ = sequest_infile.getModifications();
    return *this;
  }

  bool SequestInfile::operator==(const SequestInfile & sequest_infile) const
  {
    bool equal = true;

    equal &= (enzyme_info_ == sequest_infile.getEnzymeInfo_());
    equal &= (database_ == sequest_infile.getDatabase());
    equal &= (neutral_losses_for_ions_ == sequest_infile.getNeutralLossesForIons());
    equal &= (ion_series_weights_ == sequest_infile.getIonSeriesWeights());
    equal &= (partial_sequence_ == sequest_infile.getPartialSequence());
    equal &= (sequence_header_filter_ == sequest_infile.getSequenceHeaderFilter());
    equal &= (precursor_mass_tolerance_ == sequest_infile.getPrecursorMassTolerance());
    equal &= (peak_mass_tolerance_ == sequest_infile.getPeakMassTolerance());
    equal &= (ion_cutoff_percentage_ == sequest_infile.getIonCutoffPercentage());
    equal &= (protein_mass_filter_ == sequest_infile.getProteinMassFilter());
    equal &= (match_peak_tolerance_ == sequest_infile.getMatchPeakTolerance());
    equal &= (peptide_mass_unit_ == sequest_infile.getPeptideMassUnit());
    equal &= (output_lines_ == sequest_infile.getOutputLines());
    equal &= (enzyme_number_ == sequest_infile.getEnzymeNumber());
    equal &= (max_AA_per_mod_per_peptide_ == sequest_infile.getMaxAAPerModPerPeptide());
    equal &= (max_mods_per_peptide_ == sequest_infile.getMaxModsPerPeptide());
    equal &= (nucleotide_reading_frame_ == sequest_infile.getNucleotideReadingFrame());
    equal &= (max_internal_cleavage_sites_ == sequest_infile.getMaxInternalCleavageSites());
    equal &= (match_peak_count_ == sequest_infile.getMatchPeakCount());
    equal &= (match_peak_allowed_error_ == sequest_infile.getMatchPeakAllowedError());
    equal &= (show_fragment_ions_ == sequest_infile.getShowFragmentIons());
    equal &= (print_duplicate_references_ == sequest_infile.getPrintDuplicateReferences());
    equal &= (remove_precursor_near_peaks_ == sequest_infile.getRemovePrecursorNearPeaks());
    equal &= (mass_type_parent_ == sequest_infile.getMassTypeParent());
    equal &= (mass_type_fragment_ == sequest_infile.getMassTypeFragment());
    equal &= (normalize_xcorr_ == sequest_infile.getNormalizeXcorr());
    equal &= (residues_in_upper_case_ == sequest_infile.getResiduesInUpperCase());
    equal &= (PTMname_residues_mass_type_ == sequest_infile.getModifications());

    return equal;
  }

  void
  SequestInfile::store(
    const String & filename)
  {
    ofstream ofs(filename.c_str());
    if (!ofs)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    stringstream file_content;

    float dyn_n_term_mod(0.0), dyn_c_term_mod(0.0), stat_n_term_mod(0.0), stat_c_term_mod(0.0), stat_n_term_prot_mod(0.0), stat_c_term_prot_mod(0.0);

    map<char, float> stat_mods, dyn_mods;
    map<char, float> * mods_p = nullptr;

    // compute the masses for the amino acids, divided into fixed and optional modifications
    float mass(0.0);
    String residues, dyn_mods_string;
    for (map<String, vector<String> >::const_iterator mods_i = PTMname_residues_mass_type_.begin(); mods_i != PTMname_residues_mass_type_.end(); ++mods_i)
    {
      if (mods_i->second[0] == "CTERM")
      {
        if (mods_i->second[2] == "OPT")
        {
          dyn_c_term_mod += mods_i->second[1].toFloat();
        }
        if (mods_i->second[2] == "FIX")
        {
          stat_c_term_mod += mods_i->second[1].toFloat();
        }
      }
      else if (mods_i->second[0] == "NTERM")
      {
        if (mods_i->second[2] == "OPT")
        {
          dyn_n_term_mod += mods_i->second[1].toFloat();
        }
        if (mods_i->second[2] == "FIX")
        {
          stat_n_term_mod += mods_i->second[1].toFloat();
        }
      }
      else if (mods_i->second[0] == "CTERM_PROT")
      {
        stat_c_term_prot_mod += mods_i->second[1].toFloat();
      }
      else if (mods_i->second[0] == "NTERM_PROT")
      {
        stat_n_term_prot_mod += mods_i->second[1].toFloat();
      }
      else
      {
        if (mods_i->second[2] == "FIX")
        {
          mods_p = &stat_mods;
        }
        else
        {
          mods_p = &dyn_mods;
        }
        mass = mods_i->second[1].toFloat();
        residues = mods_i->second[0];
        for (String::const_iterator residue_i = residues.begin(); residue_i != residues.end(); ++residue_i)
          (*mods_p)[*residue_i] += mass;
      }
    }

    // now put together all optional modifications with the same mass change
    map<float, String> dyn_mods_masses;
    for (map<char, float>::const_iterator dyn_mod_i = dyn_mods.begin(); dyn_mod_i != dyn_mods.end(); ++dyn_mod_i)
    {
      dyn_mods_masses[dyn_mod_i->second].append(1, dyn_mod_i->first);
    }
    // and write them down
    if (dyn_mods_masses.empty())
    {
      dyn_mods_string = "0 X";
    }
    else
    {
      for (map<float, String>::const_iterator dyn_mod_i = dyn_mods_masses.begin(); dyn_mod_i != dyn_mods_masses.end(); ++dyn_mod_i)
      {
        dyn_mods_string.append(String(dyn_mod_i->first) + " " + dyn_mod_i->second + " ");
      }
      dyn_mods_string.erase(dyn_mods_string.length() - 1);
    }

    // the header
    file_content << "[SEQUEST]" << "\n";

    file_content << "database_name = " << database_ << "\n";

    file_content << "peptide_mass_tolerance = " << precursor_mass_tolerance_ << "\n";

    file_content << "peptide_mass_units = " << peptide_mass_unit_ << "; 0=amu, 1=mmu, 2=ppm" << "\n";

    file_content << "ion_series = " << neutral_losses_for_ions_ << " " << ion_series_weights_  << ";nABY ABCDVWXYZ" << "\n";

    file_content << "fragment_ion_tolerance = " << peak_mass_tolerance_ << "\n";

    file_content << "num_output_lines = " << output_lines_ << "\n";

    file_content << "num_results = " << output_lines_ << "\n";

    file_content << "num_description_lines = 0" << "\n";

    file_content << "show_fragment_ions = " << show_fragment_ions_ << "\n";

    file_content << "print_duplicate_references = " << print_duplicate_references_ << "\n";

    file_content << "enzyme_number = " << enzyme_number_ << "\n";

    file_content << "diff_search_options = " << dyn_mods_string << "\n";

    file_content << "term_diff_search_options = " << dyn_n_term_mod << " " << dyn_c_term_mod << "\n";

    file_content << "remove_precursor_peak = " << remove_precursor_near_peaks_ << "\n";

    file_content << "ion_cutoff_percentage = " << ion_cutoff_percentage_ << "\n";

    file_content << "protein_mass_filter = " << protein_mass_filter_ << "\n";

    file_content << "max_differential_AA_per_mod = " << max_AA_per_mod_per_peptide_ << "\n";

    file_content << "max_differential_per_peptide = " << max_mods_per_peptide_ << "\n";

    file_content << "nucleotide_reading_frame = " << nucleotide_reading_frame_ << "; 0=protein db, 1-6, 7 = forward three, 8-reverse three, 9=all six" << "\n";

    file_content << "mass_type_parent = " << mass_type_parent_ << "; 0=average masses, 1=monoisotopic masses" << "\n";

    file_content << "mass_type_fragment = " << mass_type_fragment_ << "; 0=average masses, 1=monoisotopic masses" << "\n";

    file_content << "normalize_xcorr = " << normalize_xcorr_ << "\n";

    file_content << "max_internal_cleavage_sites = " << max_internal_cleavage_sites_ << "\n";

    file_content << "create_output_files = 1" << "\n";

    file_content << "partial_sequence = " << partial_sequence_ << "\n";

    file_content << "sequence_header_filter = " << sequence_header_filter_ << "\n";

    file_content << "match_peak_count = " << match_peak_count_ << "; number of auto-detected peaks to try matching (max 5)" << "\n";

    file_content << "match_peak_allowed_error = " << match_peak_allowed_error_ << "\n";

    file_content << "match_peak_tolerance = " << match_peak_tolerance_ << "\n";

    file_content << "residues_in_upper_case = " << residues_in_upper_case_ << "\n" << "\n" << "\n";

    file_content << "add_Nterm_peptide = " << stat_n_term_mod << "\n";

    file_content << "add_Cterm_peptide = " << stat_c_term_mod << "\n";

    file_content << "add_Nterm_protein = " << stat_n_term_prot_mod << "\n";

    file_content << "add_Cterm_protein = " << stat_c_term_prot_mod << "\n" << "\n";

    file_content << "add_G_Glycine = " << stat_mods['G'] << "; added to G - avg.  57.0519, mono.  57.02146" << "\n";
    file_content << "add_A_Alanine = " << stat_mods['A'] << "; added to A - avg.  71.0788, mono.  71.03711" << "\n";
    file_content << "add_S_Serine = " << stat_mods['S'] << "; added to S - avg.  87.0782, mono.  87.03203" << "\n";
    file_content << "add_P_Proline = " << stat_mods['P'] << "; added to P - avg.  97.1167, mono.  97.05276" << "\n";
    file_content << "add_V_Valine = " << stat_mods['V'] << "; added to V - avg.  99.1326, mono.  99.06841" << "\n";
    file_content << "add_T_Threonine = " << stat_mods['T'] << "; added to T - avg. 101.1051, mono. 101.04768" << "\n";
    file_content << "add_C_Cysteine = " << stat_mods['C'] << "; added to C - avg. 103.1388, mono. 103.00919" << "\n";
    file_content << "add_L_Leucine = " << stat_mods['L'] << "; added to L - avg. 113.1594, mono. 113.08406" << "\n";
    file_content << "add_I_Isoleucine = " << stat_mods['I'] << "; added to I - avg. 113.1594, mono. 113.08406" << "\n";
    file_content << "add_X_LorI = " << stat_mods['X'] << "; added to X - avg. 113.1594, mono. 113.08406" << "\n";
    file_content << "add_N_Asparagine = " << stat_mods['N'] << "; added to N - avg. 114.1038, mono. 114.04293" << "\n";
    file_content << "add_O_Ornithine = " << stat_mods['O'] << "; added to O - avg. 114.1472, mono  114.07931" << "\n";
    file_content << "add_B_avg_NandD = " << stat_mods['B'] << "; added to B - avg. 114.5962, mono. 114.53494" << "\n";
    file_content << "add_D_Aspartic_Acid = " << stat_mods['D'] << "; added to D - avg. 115.0886, mono. 115.02694" << "\n";
    file_content << "add_Q_Glutamine = " << stat_mods['Q'] << "; added to Q - avg. 128.1307, mono. 128.05858" << "\n";
    file_content << "add_K_Lysine = " << stat_mods['K'] << "; added to K - avg. 128.1741, mono. 128.09496" << "\n";
    file_content << "add_Z_avg_QandE = " << stat_mods['Z'] << "; added to Z - avg. 128.6231, mono. 128.55059" << "\n";
    file_content << "add_E_Glutamic_Acid = " << stat_mods['E'] << "; added to E - avg. 129.1155, mono. 129.04259" << "\n";
    file_content << "add_M_Methionine = " << stat_mods['M'] << "; added to M - avg. 131.1926, mono. 131.04049" << "\n";
    file_content << "add_H_Histidine = " << stat_mods['H'] << "; added to H - avg. 137.1411, mono. 137.05891" << "\n";
    file_content << "add_F_Phenylalanine = " << stat_mods['F'] << "; added to F - avg. 147.1766, mono. 147.06841" << "\n";
    file_content << "add_R_Arginine = " << stat_mods['R'] << "; added to R - avg. 156.1875, mono. 156.10111" << "\n";
    file_content << "add_Y_Tyrosine = " << stat_mods['Y'] << "; added to Y - avg. 163.1760, mono. 163.06333" << "\n";
    file_content << "add_W_Tryptophan = " << stat_mods['W'] << "; added to W - avg. 186.2132, mono. 186.07931" << "\n" << "\n";

    file_content << getEnzymeInfoAsString();

    ofs << file_content.str();

    ofs.close();
    ofs.clear();
  }

  const map<String, vector<String> > & SequestInfile::getEnzymeInfo_() const
  {
    return enzyme_info_;
  }

  const String
  SequestInfile::getEnzymeInfoAsString() const
  {
    stringstream ss;
    Size i(0);
    String::size_type max_name_length(0);
    String::size_type max_cut_before_length(0);
    String::size_type max_doesnt_cut_after_length(0);
    ss << "[SEQUEST_ENZYME_INFO]" << "\n";
    for (map<String, vector<String> >::const_iterator einfo_i = enzyme_info_.begin(); einfo_i != enzyme_info_.end(); ++einfo_i)
    {
      max_name_length = max(max_name_length, einfo_i->first.length());
      max_cut_before_length = max(max_cut_before_length, einfo_i->second[1].length());
      max_doesnt_cut_after_length = max(max_doesnt_cut_after_length, einfo_i->second[2].length());
    }
    for (map<String, vector<String> >::const_iterator einfo_i = enzyme_info_.begin(); einfo_i != enzyme_info_.end(); ++einfo_i, ++i)
    {
      ss << i << ".  " << einfo_i->first << String(max_name_length + 5 - einfo_i->first.length(), ' ') << einfo_i->second[0] << "     " << einfo_i->second[1] << String(max_cut_before_length + 5 - einfo_i->second[1].length(), ' ') << einfo_i->second[2] << "\n";
    }
    return String(ss.str());
  }

  void
  SequestInfile::addEnzymeInfo(vector<String> & enzyme_info)
  {
    // remove duplicates from the concerned amino acids
    set<char> aas;
    for (String::const_iterator s_i = enzyme_info[2].begin(); s_i != enzyme_info[2].end(); ++s_i)
    {
      aas.insert(*s_i);
    }
    if (enzyme_info[2].length() != aas.size())
    {
      enzyme_info[2].clear();
      enzyme_info[2].reserve(aas.size());
      for (set<char>::const_iterator aa_i = aas.begin(); aa_i != aas.end(); ++aa_i)
      {
        enzyme_info[2].append(1, *aa_i);
      }
    }

    String enzyme_name = enzyme_info[0];
    enzyme_info.erase(enzyme_info.begin());
    enzyme_info_[enzyme_name] = enzyme_info;
    enzyme_number_ = 0;
    for (std::map<String, std::vector<String> >::const_iterator einfo_i = enzyme_info_.begin(); einfo_i != enzyme_info_.end(); ++einfo_i, ++enzyme_number_)
    {
      if (einfo_i->first == enzyme_name)
        break;
    }
  }

  const String & SequestInfile::getDatabase() const
  {
    return database_;
  }

  void SequestInfile::setDatabase(const String & database)
  {
    database_ = database;
  }

  const String & SequestInfile::getNeutralLossesForIons() const
  {
    return neutral_losses_for_ions_;
  }

  void SequestInfile::setNeutralLossesForIons(const String & neutral_losses_for_ions)
  {
    neutral_losses_for_ions_ = neutral_losses_for_ions;
  }

  const String & SequestInfile::getIonSeriesWeights() const
  {
    return ion_series_weights_;
  }

  void SequestInfile::setIonSeriesWeights(const String & ion_series_weights)
  {
    ion_series_weights_ = ion_series_weights;
  }

  const String & SequestInfile::getPartialSequence() const
  {
    return partial_sequence_;
  }

  void SequestInfile::setPartialSequence(const String & partial_sequence)
  {
    partial_sequence_ = partial_sequence;
  }

  const String & SequestInfile::getSequenceHeaderFilter() const
  {
    return sequence_header_filter_;
  }

  void SequestInfile::setSequenceHeaderFilter(const String & sequence_header_filter)
  {
    sequence_header_filter_ = sequence_header_filter;
  }

  const String & SequestInfile::getProteinMassFilter() const
  {
    return protein_mass_filter_;
  }

  void SequestInfile::setProteinMassFilter(const String & protein_mass_filter)
  {
    protein_mass_filter_ = protein_mass_filter;
  }

  float SequestInfile::getPrecursorMassTolerance() const
  {
    return precursor_mass_tolerance_;
  }

  void SequestInfile::setPrecursorMassTolerance(float precursor_mass_tolerance)
  {
    precursor_mass_tolerance_ = precursor_mass_tolerance;
  }

  float SequestInfile::getPeakMassTolerance() const
  {
    return peak_mass_tolerance_;
  }

  void SequestInfile::setPeakMassTolerance(float peak_mass_tolerance)
  {
    peak_mass_tolerance_ = peak_mass_tolerance;
  }

  float SequestInfile::getMatchPeakTolerance() const
  {
    return match_peak_tolerance_;
  }

  void SequestInfile::setMatchPeakTolerance(float match_peak_tolerance)
  {
    match_peak_tolerance_ = match_peak_tolerance;
  }

  float SequestInfile::getIonCutoffPercentage() const
  {
    return ion_cutoff_percentage_;
  }

  void SequestInfile::setIonCutoffPercentage(float ion_cutoff_percentage)
  {
    ion_cutoff_percentage_ = ion_cutoff_percentage;
  }

  Size SequestInfile::getPeptideMassUnit() const
  {
    return peptide_mass_unit_;
  }

  void SequestInfile::setPeptideMassUnit(Size peptide_mass_unit)
  {
    peptide_mass_unit_ = peptide_mass_unit;
  }

  Size SequestInfile::getOutputLines() const
  {
    return output_lines_;
  }

  void SequestInfile::setOutputLines(Size output_lines)
  {
    output_lines_ = output_lines;
  }

  Size SequestInfile::getEnzymeNumber() const
  {
    return enzyme_number_;
  }

  String SequestInfile::getEnzymeName() const
  {
    map<String, vector<String> >::const_iterator einfo_i = enzyme_info_.begin();
    for (Size enzyme_number = 0; enzyme_number < enzyme_number_; ++enzyme_number, ++einfo_i)
    {
    }
    return einfo_i->first;
  }

  Size SequestInfile::setEnzyme(const String& enzyme_name)
  {
    enzyme_number_ = 0;
    map<String, vector<String> >::const_iterator einfo_i;
    for (einfo_i = enzyme_info_.begin(); einfo_i != enzyme_info_.end(); ++einfo_i, ++enzyme_number_)
    {
      if (einfo_i->first == enzyme_name)
      {
        break;
      }
    }
    return (einfo_i == enzyme_info_.end()) ? enzyme_info_.size() : 0;
  }

  Size SequestInfile::getMaxAAPerModPerPeptide() const
  {
    return max_AA_per_mod_per_peptide_;
  }

  void SequestInfile::setMaxAAPerModPerPeptide(Size max_AA_per_mod_per_peptide)
  {
    max_AA_per_mod_per_peptide_ = max_AA_per_mod_per_peptide;
  }

  Size SequestInfile::getMaxModsPerPeptide() const
  {
    return max_mods_per_peptide_;
  }

  void SequestInfile::setMaxModsPerPeptide(Size max_mods_per_peptide)
  {
    max_mods_per_peptide_ = max_mods_per_peptide;
  }

  Size SequestInfile::getNucleotideReadingFrame() const
  {
    return nucleotide_reading_frame_;
  }

  void SequestInfile::setNucleotideReadingFrame(Size nucleotide_reading_frame)
  {
    nucleotide_reading_frame_ = nucleotide_reading_frame;
  }

  Size SequestInfile::getMaxInternalCleavageSites() const
  {
    return max_internal_cleavage_sites_;
  }

  void SequestInfile::setMaxInternalCleavageSites(Size max_internal_cleavage_sites)
  {
    max_internal_cleavage_sites_ = max_internal_cleavage_sites;
  }

  Size SequestInfile::getMatchPeakCount() const
  {
    return match_peak_count_;
  }

  void SequestInfile::setMatchPeakCount(Size match_peak_count)
  {
    match_peak_count_ = match_peak_count;
  }

  Size SequestInfile::getMatchPeakAllowedError() const
  {
    return match_peak_allowed_error_;
  }

  void SequestInfile::setMatchPeakAllowedError(Size match_peak_allowed_error)
  {
    match_peak_allowed_error_ = match_peak_allowed_error;
  }

  bool SequestInfile::getShowFragmentIons() const
  {
    return show_fragment_ions_;
  }

  void SequestInfile::setShowFragmentIons(bool show_fragment_ions)
  {
    show_fragment_ions_ = show_fragment_ions;
  }

  bool SequestInfile::getPrintDuplicateReferences() const
  {
    return print_duplicate_references_;
  }

  void SequestInfile::setPrintDuplicateReferences(bool print_duplicate_references)
  {
    print_duplicate_references_ = print_duplicate_references;
  }

  bool SequestInfile::getRemovePrecursorNearPeaks() const
  {
    return remove_precursor_near_peaks_;
  }

  void SequestInfile::setRemovePrecursorNearPeaks(bool remove_precursor_near_peaks)
  {
    remove_precursor_near_peaks_ = remove_precursor_near_peaks;
  }

  bool SequestInfile::getMassTypeParent() const
  {
    return mass_type_parent_;
  }

  void SequestInfile::setMassTypeParent(bool mass_type_parent)
  {
    mass_type_parent_ = mass_type_parent;
  }

  bool SequestInfile::getMassTypeFragment() const
  {
    return mass_type_fragment_;
  }

  void SequestInfile::setMassTypeFragment(bool mass_type_fragment)
  {
    mass_type_fragment_ = mass_type_fragment;
  }

  bool SequestInfile::getNormalizeXcorr() const
  {
    return normalize_xcorr_;
  }

  void SequestInfile::setNormalizeXcorr(bool normalize_xcorr)
  {
    normalize_xcorr_ = normalize_xcorr;
  }

  bool SequestInfile::getResiduesInUpperCase() const
  {
    return residues_in_upper_case_;
  }

  void SequestInfile::setResiduesInUpperCase(bool residues_in_upper_case)
  {
    residues_in_upper_case_ = residues_in_upper_case;
  }

  const map<String, vector<String> > & SequestInfile::getModifications() const
  {
    return PTMname_residues_mass_type_;
  }

  void SequestInfile::handlePTMs(const String & modification_line, const String & modifications_filename, const bool monoisotopic)
  {
    PTMname_residues_mass_type_.clear();
    // to store the information about modifications from the ptm xml file
    map<String, pair<String, String> > ptm_informations;
    if (!modification_line.empty())       // if modifications are used look whether whether composition and residues (and type and name) is given, the name (type) is used (then the modifications file is needed) or only the mass and residues (and type and name) is given
    {
      vector<String> modifications, mod_parts;
      modification_line.split(':', modifications);       // get the single modifications

      // to get masses from a formula
      EmpiricalFormula add_formula, substract_formula;

      String types = "OPT#FIX#";
      String name, residues, mass, type;

      // 0 - mass; 1 - composition; 2 - ptm name
      Int mass_or_composition_or_name(-1);

      for (const String& mod_i : modifications)
      {
        if (mod_i.empty())
        {
          continue;
        }
        // clear the formulae
        add_formula = substract_formula = EmpiricalFormula();
        name = residues = mass = type = "";

        // get the single parts of the modification string
        mod_i.split(',', mod_parts);
        mass_or_composition_or_name = -1;

        // check whether the first part is a mass, composition or name

        // check whether it is a mass
        try
        {
          mass = mod_parts.front();
          // to check whether the first part is a mass, it is converted into a float and then back into a string and compared to the given string
          // remove + signs because they don't appear in a float
          if (mass.hasPrefix("+"))
          {
            mass.erase(0, 1);
          }
          if (mass.hasSuffix("+"))
          {
            mass.erase(mass.length() - 1, 1);
          }
          if (mass.hasSuffix("-"))             // a - sign at the end will not be converted
          {
            mass.erase(mass.length() - 1, 1);
            mass.insert(0, "-");
          }
          // if it is a mass
          if (!String(mass.toFloat()).empty()) // just check if conversion does not throw, i.e. consumes the whole string
          {
            mass_or_composition_or_name = 0;
          }
        }
        catch (Exception::ConversionError & /*c_e*/)
        {
          mass_or_composition_or_name = -1;
        }

        // check whether it is a name (look it up in the corresponding file)
        if (mass_or_composition_or_name == -1)
        {
          if (ptm_informations.empty())             // if the ptm xml file has not been read yet, read it
          {
            if (!File::exists(modifications_filename))
            {
              throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, modifications_filename);
            }
            if (!File::readable(modifications_filename))
            {
              throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, modifications_filename);
            }

            // getting all available modifications from a file
            PTMXMLFile().load(modifications_filename, ptm_informations);
          }
          // if the modification cannot be found
          if (ptm_informations.find(mod_parts.front()) != ptm_informations.end())
          {
            mass = ptm_informations[mod_parts.front()].first;             // composition
            residues = ptm_informations[mod_parts.front()].second;             // residues
            name = mod_parts.front();             // name

            mass_or_composition_or_name = 2;
          }
        }

        // check whether it's an empirical formula / if a composition was given, get the mass
        if (mass_or_composition_or_name == -1)
        {
          mass = mod_parts.front();
        }
        if (mass_or_composition_or_name == -1 || mass_or_composition_or_name == 2)
        {
          // check whether there is a positive and a negative formula
          String::size_type pos = mass.find("-");
          try
          {
            if (pos != String::npos)
            {
              add_formula = EmpiricalFormula(mass.substr(0, pos));
              substract_formula = EmpiricalFormula(mass.substr(++pos));
            }
            else
            {
              add_formula = EmpiricalFormula(mass);
            }
            // sum up the masses
            if (monoisotopic)
            {
              mass = String(add_formula.getMonoWeight() - substract_formula.getMonoWeight());
            }
            else
            {
              mass = String(add_formula.getAverageWeight() - substract_formula.getAverageWeight());
            }
            if (mass_or_composition_or_name == -1)
            {
              mass_or_composition_or_name = 1;
            }
          }
          catch (Exception::ParseError & /*pe*/)
          {
            PTMname_residues_mass_type_.clear();
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mod_i, "There's something wrong with this modification. Aborting!");
          }
        }

        // now get the residues
        mod_parts.erase(mod_parts.begin());
        if (mass_or_composition_or_name < 2)
        {
          if (mod_parts.empty())
          {
            PTMname_residues_mass_type_.clear();
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mod_i, "No residues for modification given. Aborting!");
          }

          // get the residues
          residues = mod_parts.front();
          residues.substitute('*', 'X');
          residues.toUpper();
          mod_parts.erase(mod_parts.begin());
        }

        // get the type
        if (mod_parts.empty())
        {
          type = "OPT";
        }
        else
        {
          type = mod_parts.front();
          type.toUpper();
          if (types.find(type) != String::npos)
          {
            mod_parts.erase(mod_parts.begin());
          }
          else
          {
            type = "OPT";
          }
        }

        if (mod_parts.size() > 1)
        {
          PTMname_residues_mass_type_.clear();
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mod_i, "There's something wrong with the type of this modification. Aborting!");
        }

        // get the name
        if (mass_or_composition_or_name < 2)
        {
          if (mod_parts.empty())
          {
            name = "PTM_" + String(PTMname_residues_mass_type_.size());
          }
          else
          {
            name = mod_parts.front();
          }
        }

        // insert the modification
        if (PTMname_residues_mass_type_.find(name) == PTMname_residues_mass_type_.end())
        {
          PTMname_residues_mass_type_[name] = vector<String>(3);
          PTMname_residues_mass_type_[name][0] = residues;
          // mass must not have more than 5 digits after the . (otherwise the test may fail)
          PTMname_residues_mass_type_[name][1] = mass.substr(0, mass.find(".") + 6);
          PTMname_residues_mass_type_[name][2] = type;
        }
        else
        {
          PTMname_residues_mass_type_.clear();
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mod_i, "There's already a modification with this name. Aborting!");
        }
      }
    }
  }

  void SequestInfile::setStandardEnzymeInfo_()
  {
    vector<String> info;
    // cuts n to c? cuts before doesn't cut after
    info.emplace_back("0"); info.emplace_back("-"); info.emplace_back("-"); enzyme_info_["No_Enzyme"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("KR"); info.emplace_back("-"); enzyme_info_["Trypsin_Strict"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("KRLNH"); info.emplace_back("-"); enzyme_info_["Trypsin"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("FWYL"); info.emplace_back("-"); enzyme_info_["Chymotrypsin"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("FWY"); info.emplace_back("-"); enzyme_info_["Chymotrypsin_WYF"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("R"); info.emplace_back("-"); enzyme_info_["Clostripain"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("M"); info.emplace_back("-"); enzyme_info_["Cyanogen_Bromide"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("W"); info.emplace_back("-"); enzyme_info_["IodosoBenzoate"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("P"); info.emplace_back("-"); enzyme_info_["Proline_Endopept"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("E"); info.emplace_back("-"); enzyme_info_["GluC"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("ED"); info.emplace_back("-"); enzyme_info_["GluC_ED"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("K"); info.emplace_back("-"); enzyme_info_["LysC"] = info; info.clear();
    info.emplace_back("0"); info.emplace_back("D"); info.emplace_back("-"); enzyme_info_["AspN"] = info; info.clear();
    info.emplace_back("0"); info.emplace_back("DE"); info.emplace_back("-"); enzyme_info_["AspN_DE"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("ALIV"); info.emplace_back("P"); enzyme_info_["Elastase"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("ALIVKRWFY"); info.emplace_back("P"); enzyme_info_["Elastase/Tryp/Chymo"] = info; info.clear();
    info.emplace_back("1"); info.emplace_back("KRLFWYN"); info.emplace_back("-"); enzyme_info_["Trypsin/Chymo"] = info; info.clear();
  }

}
