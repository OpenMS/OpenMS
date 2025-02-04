// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Timo Sachsenberg, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OpenPepXLLFAlgorithm.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/ANALYSIS/XLMS/XQuestScores.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <OpenMS/CHEMISTRY/Tagger.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>

#include <iostream>

using namespace std;
using namespace OpenMS;

// turn on additional debug output
// #define DEBUG_OPENPEPXLLFALGO

#ifdef _OPENMP
#include <omp.h>
#endif

  OpenPepXLLFAlgorithm::OpenPepXLLFAlgorithm()
    : DefaultParamHandler("OpenPepXLLFAlgorithm")
  {
    defaults_.setValue("decoy_string", "DECOY_", "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.");
    std::vector<std::string> bool_strings = {"true", "false"};
    defaults_.setValue("decoy_prefix", "true", "Set to true, if the decoy_string is a prefix of accessions in the protein database. Otherwise it is a suffix.");
    defaults_.setValidStrings("decoy_prefix", {"true", "false"});

    defaults_.setValue("precursor:mass_tolerance", 10.0, "Width of precursor mass tolerance window");
    std::vector<std::string> mass_tolerance_unit_valid_strings = {"ppm", "Da"};
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", mass_tolerance_unit_valid_strings);
    defaults_.setValue("precursor:min_charge", 2, "Minimum precursor charge to be considered.");
    defaults_.setValue("precursor:max_charge", 8, "Maximum precursor charge to be considered.");
    defaults_.setValue("precursor:corrections", IntList({4, 3, 2, 1, 0}), "Monoisotopic peak correction. Matches candidates for possible monoisotopic precursor peaks for experimental mass m and given numbers n at masses (m - n * (C13-C12)). These should be ordered from more extreme to less extreme corrections. Numbers later in the list will be preferred in case of ambiguities.");
    defaults_.setSectionDescription("precursor", "Precursor filtering settings");

    defaults_.setValue("fragment:mass_tolerance", 20.0, "Fragment mass tolerance");
    defaults_.setValue("fragment:mass_tolerance_xlinks", 20.0, "Fragment mass tolerance for cross-link ions");
    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", mass_tolerance_unit_valid_strings);
    defaults_.setSectionDescription("fragment", "Fragment peak matching settings");

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("modifications:fixed", std::vector<std::string>({"Carbamidomethyl (C)"}), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>({"Oxidation (M)"}), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
    defaults_.setValidStrings("modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable_max_per_peptide", 3, "Maximum number of residues carrying a variable modification per candidate peptide");
    defaults_.setSectionDescription("modifications", "Peptide modification settings");

    defaults_.setValue("peptide:min_size", 5, "Minimum size a peptide must have after digestion to be considered in the search.");
    defaults_.setValue("peptide:missed_cleavages", 3, "Number of missed cleavages.");
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);

    defaults_.setValue("peptide:enzyme", "Trypsin", "The enzyme used for peptide digestion.");
    defaults_.setValidStrings("peptide:enzyme", ListUtils::create<std::string>(all_enzymes));
    defaults_.setSectionDescription("peptide", "Settings for digesting proteins into peptides");

    defaults_.setValue("cross_linker:residue1", std::vector<std::string>({"K", "N-term"}), "Comma separated residues, that the first side of a bifunctional cross-linker can attach to");
    defaults_.setValue("cross_linker:residue2", std::vector<std::string>({"K", "N-term"}), "Comma separated residues, that the second side of a bifunctional cross-linker can attach to");
    defaults_.setValue("cross_linker:mass", 138.0680796, "Mass of the light cross-linker, linking two residues on one or two peptides");
    defaults_.setValue("cross_linker:mass_mono_link", DoubleList({156.07864431, 155.094628715}), "Possible masses of the linker, when attached to only one peptide");
    defaults_.setValue("cross_linker:name", "DSS", "Name of the searched cross-link, used to resolve ambiguity of equal masses (e.g. DSS or BS3)");
    defaults_.setSectionDescription("cross_linker", "Description of the cross-linker reagent");

    defaults_.setValue("algorithm:number_top_hits", 1, "Number of top hits reported for each spectrum pair");
    std::vector<std::string> deisotope_strings = std::vector<std::string>({"true", "false", "auto"});
    defaults_.setValue("algorithm:deisotope", "auto", "Set to true, if the input spectra should be deisotoped before any other processing steps. If set to auto the spectra will be deisotoped, if the fragment mass tolerance is < 0.1 Da or < 100 ppm (0.1 Da at a mass of 1000)", std::vector<std::string>({"advanced"}));
    defaults_.setValidStrings("algorithm:deisotope", deisotope_strings);
    defaults_.setValue("algorithm:use_sequence_tags", "false", "Use sequence tags (de novo sequencing of short fragments) to filter out candidates before scoring. This will make the search faster, but can impact the sensitivity positively or negatively, depending on the dataset.");
    defaults_.setValidStrings("algorithm:use_sequence_tags", bool_strings);
    defaults_.setValue("algorithm:sequence_tag_min_length", 2, "Minimal length of sequence tags to use for filtering candidates. Longer tags will make the search faster but much less sensitive. Ignored if 'algorithm:use_sequence_tags' is false.", std::vector<std::string>({"advanced"}));
    defaults_.setSectionDescription("algorithm", "Additional algorithm settings");

    defaults_.setValue("ions:b_ions", "true", "Search for peaks of b-ions.", std::vector<std::string>({"advanced"}));
    defaults_.setValue("ions:y_ions", "true", "Search for peaks of y-ions.", std::vector<std::string>({"advanced"}));
    defaults_.setValue("ions:a_ions", "false", "Search for peaks of a-ions.", std::vector<std::string>({"advanced"}));
    defaults_.setValue("ions:x_ions", "false", "Search for peaks of x-ions.", std::vector<std::string>({"advanced"}));
    defaults_.setValue("ions:c_ions", "false", "Search for peaks of c-ions.", std::vector<std::string>({"advanced"}));
    defaults_.setValue("ions:z_ions", "false", "Search for peaks of z-ions.", std::vector<std::string>({"advanced"}));
    defaults_.setValue("ions:neutral_losses", "true", "Search for neutral losses of H2O and H3N.", std::vector<std::string>({"advanced"}));
    defaults_.setValidStrings("ions:b_ions", bool_strings);
    defaults_.setValidStrings("ions:y_ions", bool_strings);
    defaults_.setValidStrings("ions:a_ions", bool_strings);
    defaults_.setValidStrings("ions:x_ions", bool_strings);
    defaults_.setValidStrings("ions:c_ions", bool_strings);
    defaults_.setValidStrings("ions:z_ions", bool_strings);
    defaults_.setValidStrings("ions:neutral_losses", bool_strings);
    defaults_.setSectionDescription("ions", "Ion types to search for in MS/MS spectra");

    defaultsToParam_();
  }

  OpenPepXLLFAlgorithm::~OpenPepXLLFAlgorithm() = default;

  void OpenPepXLLFAlgorithm::updateMembers_()
  {
    decoy_string_ = static_cast<String>(param_.getValue("decoy_string").toString());
    decoy_prefix_ = param_.getValue("decoy_prefix") == "true";

    min_precursor_charge_ = static_cast<Int>(param_.getValue("precursor:min_charge"));
    max_precursor_charge_ = static_cast<Int>(param_.getValue("precursor:max_charge"));
    precursor_mass_tolerance_ = static_cast<double>(param_.getValue("precursor:mass_tolerance"));
    precursor_mass_tolerance_unit_ppm_ = (static_cast<String>(param_.getValue("precursor:mass_tolerance_unit").toString()) == "ppm");
    precursor_correction_steps_ = param_.getValue("precursor:corrections");

    fragment_mass_tolerance_ = static_cast<double>(param_.getValue("fragment:mass_tolerance"));
    fragment_mass_tolerance_xlinks_ = static_cast<double>(param_.getValue("fragment:mass_tolerance_xlinks"));
    fragment_mass_tolerance_unit_ppm_  = (static_cast<String>(param_.getValue("fragment:mass_tolerance_unit").toString()) == "ppm");

    cross_link_residue1_ = ListUtils::toStringList<std::string>(param_.getValue("cross_linker:residue1"));
    cross_link_residue2_ = ListUtils::toStringList<std::string>(param_.getValue("cross_linker:residue2"));
    cross_link_mass_ = static_cast<double>(param_.getValue("cross_linker:mass"));
    cross_link_mass_mono_link_ = param_.getValue("cross_linker:mass_mono_link");
    cross_link_name_ = static_cast<String>(param_.getValue("cross_linker:name").toString());

    fixedModNames_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:fixed"));
    varModNames_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    max_variable_mods_per_peptide_ = static_cast<Size>(param_.getValue("modifications:variable_max_per_peptide"));
    peptide_min_size_ = static_cast<Size>(param_.getValue("peptide:min_size"));
    missed_cleavages_ = static_cast<Size>(param_.getValue("peptide:missed_cleavages"));
    enzyme_name_ = static_cast<String>(param_.getValue("peptide:enzyme").toString());

    number_top_hits_ = static_cast<Int>(param_.getValue("algorithm:number_top_hits"));
    deisotope_mode_ = static_cast<String>(param_.getValue("algorithm:deisotope").toString());
    use_sequence_tags_ = param_.getValue("algorithm:use_sequence_tags") == "true";
    sequence_tag_min_length_ = static_cast<Size>(param_.getValue("algorithm:sequence_tag_min_length"));

    add_y_ions_ = param_.getValue("ions:y_ions").toString();
    add_b_ions_ = param_.getValue("ions:b_ions").toString();
    add_x_ions_ = param_.getValue("ions:x_ions").toString();
    add_a_ions_ = param_.getValue("ions:a_ions").toString();
    add_c_ions_ = param_.getValue("ions:c_ions").toString();
    add_z_ions_ = param_.getValue("ions:z_ions").toString();
    add_losses_ = param_.getValue("ions:neutral_losses").toString();
  }

  OpenPepXLLFAlgorithm::ExitCodes OpenPepXLLFAlgorithm::run(PeakMap& unprocessed_spectra, std::vector<FASTAFile::FASTAEntry>& fasta_db, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, PeakMap& spectra)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(this->getLogType());

    // preprocess parameters for convenience
    if (fragment_mass_tolerance_xlinks_ < fragment_mass_tolerance_)
    {
      fragment_mass_tolerance_xlinks_ = fragment_mass_tolerance_;
    }
    std::sort(cross_link_mass_mono_link_.begin(), cross_link_mass_mono_link_.end(), std::greater< double >());
    set<String> fixed_unique(fixedModNames_.begin(), fixedModNames_.end());

    // deisotope if "true" or if "auto" and the tolerance is below the threshold (0.1 Da or 100 ppm)
    bool deisotope = (deisotope_mode_ == "true") ||
                      (deisotope_mode_ == "auto" &&
                      ((!fragment_mass_tolerance_unit_ppm_ && fragment_mass_tolerance_ < 0.1) ||
                      (fragment_mass_tolerance_unit_ppm_ && fragment_mass_tolerance_ < 100)));

    if (fixed_unique.size() != fixedModNames_.size())
    {
      OPENMS_LOG_WARN << "duplicate fixed modification provided." << endl;
      return ExitCodes::ILLEGAL_PARAMETERS;
    }

    set<String> var_unique(varModNames_.begin(), varModNames_.end());
    if (var_unique.size() != varModNames_.size())
    {
      OPENMS_LOG_WARN << "duplicate variable modification provided." << endl;
      return ExitCodes::ILLEGAL_PARAMETERS;
    }
    ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(fixedModNames_);
    ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(varModNames_);

    protein_ids[0].setPrimaryMSRunPath({}, unprocessed_spectra);

    if (unprocessed_spectra.empty() && unprocessed_spectra.getChromatograms().empty())
    {
      OPENMS_LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry." << endl;
      return INCOMPATIBLE_INPUT_DATA;
    }

    //check if spectra are sorted
    for (Size i = 0; i < unprocessed_spectra.size(); ++i)
    {
      if (!unprocessed_spectra[i].isSorted())
      {
        OPENMS_LOG_WARN << "Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!" << endl;
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    // Peak Picking, check if all levels are picked and pick uncentroided MS levels
    PeakPickerHiRes pp;
    PeakMap picked_spectra;
    progresslogger.startProgress(0, 1, "Centroiding data (if necessary)...");
    pp.pickExperiment(unprocessed_spectra, picked_spectra, true);
    progresslogger.endProgress();
    unprocessed_spectra.clear(true);

    // Precursor Purity precalculation
    map<String, PrecursorPurity::PurityScores> precursor_purities = PrecursorPurity::computePrecursorPurities(picked_spectra, precursor_mass_tolerance_, precursor_mass_tolerance_unit_ppm_);

    // preprocess spectra (filter out 0 values, sort by position)
    progresslogger.startProgress(0, 1, "Filtering spectra...");
    spectra = OPXLSpectrumProcessingAlgorithms::preprocessSpectra(picked_spectra, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, peptide_min_size_, min_precursor_charge_, max_precursor_charge_, deisotope, false);
    progresslogger.endProgress();
    picked_spectra.clear(true);

    // sort the spectra by RT, the order might have been changed by parallel preprocessing
    spectra.sortSpectra(false);

    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme_name_);
    digestor.setMissedCleavages(missed_cleavages_);

    ProteinIdentification::SearchParameters search_params = protein_ids[0].getSearchParameters();
    String searched_charges((String(min_precursor_charge_)));
    for (int ch = min_precursor_charge_+1; ch <= max_precursor_charge_; ++ch)
    {
      searched_charges += "," + String(ch);
    }
    search_params.charges = searched_charges;
    search_params.digestion_enzyme = (*ProteaseDB::getInstance()->getEnzyme(enzyme_name_));
    search_params.fixed_modifications = fixedModNames_;
    search_params.variable_modifications = varModNames_;
    search_params.mass_type = ProteinIdentification::MONOISOTOPIC;
    search_params.missed_cleavages = missed_cleavages_;
    search_params.fragment_mass_tolerance = fragment_mass_tolerance_;
    search_params.fragment_mass_tolerance_ppm =  fragment_mass_tolerance_unit_ppm_;
    search_params.precursor_mass_tolerance = precursor_mass_tolerance_;
    search_params.precursor_mass_tolerance_ppm = precursor_mass_tolerance_unit_ppm_;

    // As MetaValues
    search_params.setMetaValue("decoy_prefix", decoy_prefix_);
    search_params.setMetaValue("decoy_string", decoy_string_);

    search_params.setMetaValue("precursor:min_charge", min_precursor_charge_);
    search_params.setMetaValue("precursor:max_charge", max_precursor_charge_);

    search_params.setMetaValue("fragment:mass_tolerance_xlinks", fragment_mass_tolerance_xlinks_);
    search_params.setMetaValue("peptide:min_size", peptide_min_size_);

    search_params.setMetaValue("cross_link:residue1", cross_link_residue1_);
    search_params.setMetaValue("cross_link:residue2", cross_link_residue2_);
    search_params.setMetaValue("cross_link:mass", cross_link_mass_);
    search_params.setMetaValue("cross_link:mass_monolink", cross_link_mass_mono_link_);
    search_params.setMetaValue("cross_link:name", cross_link_name_);
    search_params.setMetaValue("precursor:corrections", precursor_correction_steps_);

    search_params.setMetaValue("modifications:variable_max_per_peptide", max_variable_mods_per_peptide_);
    protein_ids[0].setSearchParameters(search_params);
    protein_ids[0].setScoreType("OpenPepXL_Protein_Score");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    vector<OPXLDataStructs::AASeqWithMass> peptide_masses;
    peptide_masses = OPXLHelper::digestDatabase(fasta_db, digestor, peptide_min_size_, cross_link_residue1_, cross_link_residue2_, fixed_modifications,  variable_modifications, max_variable_mods_per_peptide_);

    // declare and set up spectrum generators
    TheoreticalSpectrumGeneratorXLMS specGen_full;
    SimpleTSGXLMS specGen_mainscore;

    // settings for full-scoring, annotations, 2nd isotopic peaks, losses and precursors
    Param specGenParams_full = specGen_full.getParameters();
    specGenParams_full.setValue("add_b_ions", add_b_ions_, "Add peaks of y-ions to the spectrum");
    specGenParams_full.setValue("add_y_ions", add_y_ions_, "Add peaks of b-ions to the spectrum");
    specGenParams_full.setValue("add_a_ions", add_a_ions_, "Add peaks of a-ions to the spectrum");
    specGenParams_full.setValue("add_x_ions", add_x_ions_, "Add peaks of c-ions to the spectrum");
    specGenParams_full.setValue("add_c_ions", add_c_ions_, "Add peaks of x-ions to the spectrum");
    specGenParams_full.setValue("add_z_ions", add_z_ions_, "Add peaks of z-ions to the spectrum");
    specGenParams_full.setValue("add_losses", add_losses_, "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");

    specGenParams_full.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    specGenParams_full.setValue("add_metainfo", "true");
    specGenParams_full.setValue("add_charges", "true");
    specGenParams_full.setValue("add_isotopes", "true", "If set to 1 isotope peaks of the product ion peaks are added");
    specGenParams_full.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    specGenParams_full.setValue("add_precursor_peaks", "true");
    specGenParams_full.setValue("add_k_linked_ions", "false");
    specGen_full.setParameters(specGenParams_full);

    Param specGenParams_mainscore = specGen_mainscore.getParameters();
    specGenParams_mainscore.setValue("add_b_ions", add_b_ions_, "Add peaks of y-ions to the spectrum");
    specGenParams_mainscore.setValue("add_y_ions", add_y_ions_, "Add peaks of b-ions to the spectrum");
    specGenParams_mainscore.setValue("add_a_ions", add_a_ions_, "Add peaks of a-ions to the spectrum");
    specGenParams_mainscore.setValue("add_x_ions", add_x_ions_, "Add peaks of c-ions to the spectrum");
    specGenParams_mainscore.setValue("add_c_ions", add_c_ions_, "Add peaks of x-ions to the spectrum");
    specGenParams_mainscore.setValue("add_z_ions", add_z_ions_, "Add peaks of z-ions to the spectrum");
    specGenParams_mainscore.setValue("add_losses", add_losses_, "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    specGenParams_mainscore.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    specGenParams_mainscore.setValue("add_isotopes", "true", "If set to 1 isotope peaks of the product ion peaks are added");
    specGenParams_mainscore.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");
    specGenParams_mainscore.setValue("add_precursor_peaks", "true");
    specGenParams_mainscore.setValue("add_k_linked_ions", "false");
    specGen_mainscore.setParameters(specGenParams_mainscore);

    // use a less stringent tolerance for the tagger
    double tagger_tol;
    if (fragment_mass_tolerance_unit_ppm_) // ppm: increase tolerance by 50%
    {
      tagger_tol = fragment_mass_tolerance_xlinks_ * 1.5;
    }
    else // dalton: turn into ppm at 500 mz
    {
      tagger_tol = Math::getPPM(500.0, 500+fragment_mass_tolerance_xlinks_);
    }

    Tagger tagger(sequence_tag_min_length_, tagger_tol, sequence_tag_min_length_, 1, max_precursor_charge_, fixedModNames_, varModNames_);
    Size all_candidates_count(0);

#ifdef DEBUG_OPENPEPXLLFALGO
    OPENMS_LOG_DEBUG << "Peptide candidates: " << peptide_masses.size() << endl;
#endif

    search_params = protein_ids[0].getSearchParameters();
    search_params.setMetaValue("MS:1001029", peptide_masses.size()); // number of sequences searched = MS:1001029
    protein_ids[0].setSearchParameters(search_params);

    // Collect precursor MZs for filtering enumerated peptide pairs
    vector< double > spectrum_precursors;
    for (Size i = 0; i < spectra.size(); i++)
    {
      double current_precursor_mz = spectra[i].getPrecursors()[0].getMZ();
      double current_precursor_charge = spectra[i].getPrecursors()[0].getCharge();
      double current_precursor_mass = (current_precursor_mz * current_precursor_charge) - (current_precursor_charge * Constants::PROTON_MASS_U);
      spectrum_precursors.push_back(current_precursor_mass);
    }
    sort(spectrum_precursors.begin(), spectrum_precursors.end());
    sort(peptide_masses.begin(), peptide_masses.end(), OPXLDataStructs::AASeqWithMassComparator());
    // The largest peptides given a fixed maximal precursor mass are possible with loop links
    // Filter peptides using maximal loop link mass first
    double max_precursor_mass = spectrum_precursors[spectrum_precursors.size()-1];

    // compute absolute tolerance from relative, if necessary
    double max_peptide_allowed_error = 0;
    if (precursor_mass_tolerance_unit_ppm_) // ppm
    {
      max_peptide_allowed_error = max_precursor_mass * precursor_mass_tolerance_ * 1e-6;
    }
    else // Dalton
    {
      max_peptide_allowed_error= precursor_mass_tolerance_;
    }

    double max_peptide_mass = max_precursor_mass - cross_link_mass_ + max_peptide_allowed_error;

    // search for the first mass greater than the maximum, cut off everything larger
    vector<OPXLDataStructs::AASeqWithMass>::iterator last = upper_bound(peptide_masses.begin(), peptide_masses.end(), max_peptide_mass, OPXLDataStructs::AASeqWithMassComparator());
    vector<OPXLDataStructs::AASeqWithMass> filtered_peptide_masses;
    filtered_peptide_masses.assign(peptide_masses.begin(), last);

    // iterate over all spectra
    progresslogger.startProgress(0, 1, "Matching to theoretical spectra and scoring...");

    Size spectrum_counter(0);

    for (SignedSize scan_index = 0; scan_index < static_cast<SignedSize>(spectra.size()); ++scan_index)
    {
      const PeakSpectrum& spectrum = spectra[scan_index];

      const double precursor_charge = spectrum.getPrecursors()[0].getCharge();
      const double precursor_mz = spectrum.getPrecursors()[0].getMZ();
      const double precursor_mass = (precursor_mz * static_cast<double>(precursor_charge)) - (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U);

      std::vector<std::string> tags;
      if (use_sequence_tags_)
      {
        tagger.setMaxCharge(precursor_charge-1);
        tagger.getTag(spectrum, tags);
      }

      vector< OPXLDataStructs::CrossLinkSpectrumMatch > top_csms_spectrum;
      vector< OPXLDataStructs::ProteinProteinCrossLink > cross_link_candidates = OPXLHelper::collectPrecursorCandidates(precursor_correction_steps_, precursor_mass, precursor_mass_tolerance_, precursor_mass_tolerance_unit_ppm_, filtered_peptide_masses, cross_link_mass_, cross_link_mass_mono_link_, cross_link_residue1_, cross_link_residue2_, cross_link_name_, use_sequence_tags_, tags);
      all_candidates_count += cross_link_candidates.size();

#ifdef DEBUG_OPENPEPXLLFALGO
      OPENMS_LOG_DEBUG << "Size of enumerated candidates: " << double(cross_link_candidates.size()) * sizeof(OPXLDataStructs::ProteinProteinCrossLink) / 1024.0 / 1024.0 << " mb" << endl;
#endif

      spectrum_counter++;
      cout << "Processing spectrum " << spectrum_counter << " / " << spectra.size() << " |\tSpectrum ID: " << spectrum.getNativeID() << "\t| at: " << DateTime::now().getTime() << endl;
      cout << "Number of peaks: " << spectrum.size() << " |\tNumber of candidates: " << cross_link_candidates.size() << endl;

      if (cross_link_candidates.empty())
      {
        continue;
      }

      // lists for one spectrum, to determine best match to the spectrum
      vector< OPXLDataStructs::CrossLinkSpectrumMatch > all_csms_spectrum;
      vector< OPXLDataStructs::CrossLinkSpectrumMatch > mainscore_csms_spectrum;

#pragma omp parallel for schedule(guided)
      for (SignedSize i = 0; i < static_cast<SignedSize>(cross_link_candidates.size()); ++i)
      {
        OPXLDataStructs::ProteinProteinCrossLink cross_link_candidate = cross_link_candidates[i];

        std::vector< SimpleTSGXLMS::SimplePeak > theoretical_spec_linear_alpha;
        theoretical_spec_linear_alpha.reserve(1500);
        std::vector< SimpleTSGXLMS::SimplePeak > theoretical_spec_linear_beta;
        std::vector< SimpleTSGXLMS::SimplePeak > theoretical_spec_xlinks_alpha;
        std::vector< SimpleTSGXLMS::SimplePeak > theoretical_spec_xlinks_beta;

        bool type_is_cross_link = cross_link_candidate.getType() == OPXLDataStructs::CROSS;
        bool type_is_loop = cross_link_candidate.getType() == OPXLDataStructs::LOOP;
        Size link_pos_B = 0;
        if (type_is_loop)
        {
          link_pos_B = cross_link_candidate.cross_link_position.second;
        }
        AASequence alpha;
        AASequence beta;
        if (cross_link_candidate.alpha) { alpha = *cross_link_candidate.alpha; }
        if (cross_link_candidate.beta) { beta = *cross_link_candidate.beta; }

        specGen_mainscore.getLinearIonSpectrum(theoretical_spec_linear_alpha, alpha, cross_link_candidate.cross_link_position.first, 2, link_pos_B);
        if (type_is_cross_link)
        {
          theoretical_spec_linear_beta.reserve(1500);
          specGen_mainscore.getLinearIonSpectrum(theoretical_spec_linear_beta, beta, cross_link_candidate.cross_link_position.second, 2);
        }

        // Something like this can happen, e.g. with a loop link connecting the first and last residue of a peptide
        if ( theoretical_spec_linear_alpha.empty() )
        {
          continue;
        }

        vector< pair< Size, Size > > matched_spec_linear_alpha;
        vector< pair< Size, Size > > matched_spec_linear_beta;
        vector< pair< Size, Size > > matched_spec_xlinks_alpha;
        vector< pair< Size, Size > > matched_spec_xlinks_beta;

        PeakSpectrum::IntegerDataArray exp_charges;
        if (!spectrum.getIntegerDataArrays().empty())
        {
          exp_charges = spectrum.getIntegerDataArrays()[0];
        }
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(matched_spec_linear_alpha, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_linear_alpha, spectrum, exp_charges);
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(matched_spec_linear_beta, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_linear_beta, spectrum, exp_charges);

        // drop candidates with almost no linear fragment peak matches before making the more complex theoretical spectra and aligning them
        // this removes hits that no one would trust after manual validation anyway and reduces time wasted on really bad spectra or candidates without any matching peaks
        if (matched_spec_linear_alpha.size() < 2 || (type_is_cross_link && matched_spec_linear_beta.size() < 2) )
        {
          continue;
        }
        theoretical_spec_xlinks_alpha.reserve(1500);

        if (type_is_cross_link)
        {
          theoretical_spec_xlinks_beta.reserve(1500);
          specGen_mainscore.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate, true, 2, precursor_charge);
          specGen_mainscore.getXLinkIonSpectrum(theoretical_spec_xlinks_beta, cross_link_candidate, false, 2, precursor_charge);
        }
        else
        {
          // Function for mono-links or loop-links
          specGen_mainscore.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, alpha, cross_link_candidate.cross_link_position.first, precursor_mass, 1, precursor_charge, link_pos_B);
        }
        if (theoretical_spec_xlinks_alpha.empty())
        {
          continue;
        }

        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(matched_spec_xlinks_alpha, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_xlinks_alpha, spectrum, exp_charges);
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(matched_spec_xlinks_beta, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_xlinks_beta, spectrum, exp_charges);

        // the maximal xlink ion charge is (precursor charge - 1) and the minimal xlink ion charge is 2.
        // we need the difference between min and max here, which is (precursor_charge - 3) in most cases
        // but we also need a number > 0, we set 1 as the minimum, in case the precursor charge is only 3 or smaller
        Size n_xlink_charges = 1;
        if (precursor_charge > 3)
        {
          n_xlink_charges = precursor_charge - 3;
        }

        // compute match odds (unweighted), the 3 is the number of charge states in the theoretical spectra
        double match_odds_c_alpha = XQuestScores::matchOddsScoreSimpleSpec(theoretical_spec_linear_alpha, matched_spec_linear_alpha.size(), fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_);
        double match_odds_x_alpha = XQuestScores::matchOddsScoreSimpleSpec(theoretical_spec_xlinks_alpha, matched_spec_xlinks_alpha.size(), fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, true, n_xlink_charges);
        double match_odds = 0;
        double match_odds_alpha = 0;
        double match_odds_beta = 0;

        if (type_is_cross_link)
        {
          double match_odds_c_beta = XQuestScores::matchOddsScoreSimpleSpec(theoretical_spec_linear_beta, matched_spec_linear_beta.size(), fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_);
          double match_odds_x_beta = XQuestScores::matchOddsScoreSimpleSpec(theoretical_spec_xlinks_beta, matched_spec_xlinks_beta.size(), fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, true, n_xlink_charges);
          match_odds = (match_odds_c_alpha + match_odds_x_alpha + match_odds_c_beta + match_odds_x_beta) / 4;
          match_odds_alpha = (match_odds_c_alpha + match_odds_x_alpha) / 2;
          match_odds_beta = (match_odds_c_beta + match_odds_x_beta) / 2;
        }
        else
        {
          match_odds = (match_odds_c_alpha + match_odds_x_alpha) / 2;
          match_odds_alpha = match_odds;
        }

        OPXLDataStructs::CrossLinkSpectrumMatch csm;
        csm.cross_link = cross_link_candidate;
        csm.precursor_correction = cross_link_candidate.precursor_correction;
        double rel_error = OPXLHelper::computePrecursorError(csm, precursor_mz, precursor_charge);

        double new_match_odds_weight = 0.2;
        double new_rel_error_weight = -0.03;
        double new_score = new_match_odds_weight * std::log(1e-7 + match_odds) + new_rel_error_weight * abs(rel_error);

        csm.score = new_score;
        csm.match_odds = match_odds;
        csm.match_odds_alpha = match_odds_alpha;
        csm.match_odds_beta = match_odds_beta;
        csm.precursor_error_ppm = rel_error;

#pragma omp critical (mainscore_csms_spectrum_access)
        mainscore_csms_spectrum.push_back(csm);

      }
      std::sort(mainscore_csms_spectrum.rbegin(), mainscore_csms_spectrum.rend(), OPXLDataStructs::CLSMScoreComparator());

      int last_candidate_index = static_cast<int>(mainscore_csms_spectrum.size());
      last_candidate_index = std::min(last_candidate_index, number_top_hits_);

#pragma omp parallel for schedule(guided)
      for (int i = 0; i < last_candidate_index ; ++i)
      {
        OPXLDataStructs::ProteinProteinCrossLink cross_link_candidate = mainscore_csms_spectrum[i].cross_link;
        AASequence alpha;
        AASequence beta;
        if (cross_link_candidate.alpha) { alpha = *cross_link_candidate.alpha; }
        if (cross_link_candidate.beta) { beta = *cross_link_candidate.beta; }

#ifdef DEBUG_OPENPEPXLLFALGO
        double candidate_mz = (alpha.getMonoWeight() + beta.getMonoWeight() +  cross_link_candidate.cross_linker_mass + (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / precursor_charge;
#pragma omp critical (LOG_DEBUG_access)
        {
          OPENMS_LOG_DEBUG << "Pair: " << alpha.toString() << "-" << beta.toString() << " matched to light spectrum " << scan_index << "\t and heavy spectrum " << scan_index
            << " with m/z: " << precursor_mz << "\t" << "and candidate m/z: " << candidate_mz << "\tK Positions: " << cross_link_candidate.cross_link_position.first << "\t" << cross_link_candidate.cross_link_position.second << endl;
        }
#endif
        OPXLDataStructs::CrossLinkSpectrumMatch csm = mainscore_csms_spectrum[i];

        PeakSpectrum theoretical_spec_linear_alpha;
        theoretical_spec_linear_alpha.reserve(1500);
        PeakSpectrum theoretical_spec_linear_beta;
        PeakSpectrum theoretical_spec_xlinks_alpha;
        theoretical_spec_xlinks_alpha.reserve(1500);
        PeakSpectrum theoretical_spec_xlinks_beta;

        bool type_is_cross_link = cross_link_candidate.getType() == OPXLDataStructs::CROSS;
        bool type_is_loop = cross_link_candidate.getType() == OPXLDataStructs::LOOP;
        Size link_pos_B = 0;
        if (type_is_loop)
        {
          link_pos_B = cross_link_candidate.cross_link_position.second;
        }
        specGen_full.getLinearIonSpectrum(theoretical_spec_linear_alpha, alpha, cross_link_candidate.cross_link_position.first, true, 2, link_pos_B);
        if (type_is_cross_link)
        {
          theoretical_spec_linear_beta.reserve(1500);
          theoretical_spec_xlinks_beta.reserve(1500);
          specGen_full.getLinearIonSpectrum(theoretical_spec_linear_beta, beta, cross_link_candidate.cross_link_position.second, false, 2);
          specGen_full.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate, true, 1, precursor_charge);
          specGen_full.getXLinkIonSpectrum(theoretical_spec_xlinks_beta, cross_link_candidate, false, 1, precursor_charge);
        }
        else
        {
          // Function for mono-links or loop-links
          specGen_full.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, alpha, cross_link_candidate.cross_link_position.first, precursor_mass, true, 2, precursor_charge, link_pos_B);
        }

        // Something like this can happen, e.g. with a loop link connecting the first and last residue of a peptide
        if ( (theoretical_spec_linear_alpha.empty()) || (theoretical_spec_xlinks_alpha.empty()) )
        {
          continue;
        }

        vector< pair< Size, Size > > matched_spec_linear_alpha;
        vector< pair< Size, Size > > matched_spec_linear_beta;
        vector< pair< Size, Size > > matched_spec_xlinks_alpha;
        vector< pair< Size, Size > > matched_spec_xlinks_beta;

        DataArrays::FloatDataArray ppm_error_array_linear_alpha;
        DataArrays::FloatDataArray ppm_error_array_xlinks_alpha;
        DataArrays::FloatDataArray ppm_error_array_linear_beta;
        DataArrays::FloatDataArray ppm_error_array_xlinks_beta;

        PeakSpectrum::IntegerDataArray& theo_charges_la = theoretical_spec_linear_alpha.getIntegerDataArrays()[0];
        PeakSpectrum::IntegerDataArray theo_charges_xa;
        if (!theoretical_spec_xlinks_alpha.getIntegerDataArrays().empty())
        {
          theo_charges_xa = theoretical_spec_xlinks_alpha.getIntegerDataArrays()[0];
        }
        PeakSpectrum::IntegerDataArray theo_charges_lb;
        PeakSpectrum::IntegerDataArray theo_charges_xb;
        if (!theoretical_spec_linear_beta.getIntegerDataArrays().empty())
        {
          theo_charges_lb = theoretical_spec_linear_beta.getIntegerDataArrays()[0];
        }
        if (!theoretical_spec_xlinks_beta.getIntegerDataArrays().empty())
        {
          theo_charges_xb = theoretical_spec_xlinks_beta.getIntegerDataArrays()[0];
        }
        PeakSpectrum::IntegerDataArray exp_charges;
        if (!spectrum.getIntegerDataArrays().empty())
        {
          exp_charges = spectrum.getIntegerDataArrays()[0];
        }
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_linear_alpha, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_linear_alpha, spectrum, theo_charges_la, exp_charges, ppm_error_array_linear_alpha);
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_linear_beta, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_linear_beta, spectrum, theo_charges_lb, exp_charges, ppm_error_array_linear_beta);
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_xlinks_alpha, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_xlinks_alpha, spectrum, theo_charges_xa, exp_charges, ppm_error_array_xlinks_alpha);
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_xlinks_beta, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_xlinks_beta, spectrum, theo_charges_xb, exp_charges, ppm_error_array_xlinks_beta);

#ifdef DEBUG_OPENPEPXLLFALGO
#pragma omp critical (LOG_DEBUG_access)
        {
          OPENMS_LOG_DEBUG << "Spectrum sizes: " << spectrum.size() << " || " << theoretical_spec_linear_alpha.size() <<  " | " << theoretical_spec_linear_beta.size()
                                <<  " | " << theoretical_spec_xlinks_alpha.size() <<  " | " << theoretical_spec_xlinks_beta.size() << endl;
          OPENMS_LOG_DEBUG << "Matched peaks: " << matched_spec_linear_alpha.size() << " | " << matched_spec_linear_beta.size()
                                <<  " | " << matched_spec_xlinks_alpha.size() <<  " | " << matched_spec_xlinks_beta.size() << endl;
        }
#endif

        // TODO define good exclusion criteria for total crap
        Size matched_peaks = matched_spec_linear_alpha.size() + matched_spec_linear_beta.size() + matched_spec_xlinks_alpha.size() + matched_spec_xlinks_beta.size();
        if (matched_peaks < 1)
        {
          continue;
        }

        // compute intsum score
        double intsum = XQuestScores::totalMatchedCurrent(matched_spec_linear_alpha, matched_spec_linear_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta, spectrum, spectrum);

        // Total ion intensity of light spectrum
        // sum over linear and xlink ion spectra instead of unfiltered
        double total_current = 0;
        for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum.size()); ++j)
        {
          total_current += spectrum[j].getIntensity();
        }
        double TIC = intsum / total_current;

        // TIC_alpha and _beta (total ion current)
        double intsum_alpha = XQuestScores::matchedCurrentChain(matched_spec_linear_alpha, matched_spec_xlinks_alpha, spectrum, spectrum);
        double intsum_beta = 0;
        if (type_is_cross_link)
        {
          intsum_beta = XQuestScores::matchedCurrentChain(matched_spec_linear_beta, matched_spec_xlinks_beta, spectrum, spectrum);
        }

        // normalize TIC_alpha and  _beta
        if ((intsum_alpha + intsum_beta) > 0.0)
        {
          intsum_alpha = intsum_alpha * intsum / (intsum_alpha + intsum_beta);
          intsum_beta = intsum_beta *  intsum / (intsum_alpha + intsum_beta);
        }

        // compute weighted TIC
        double wTIC = XQuestScores::weightedTICScore(alpha.size(), beta.size(), intsum_alpha, intsum_beta, total_current, type_is_cross_link);
        double wTICold = XQuestScores::weightedTICScoreXQuest(alpha.size(), beta.size(), intsum_alpha, intsum_beta, total_current, type_is_cross_link);

        // compute match odds (unweighted), the 3 is the number of charge states in the theoretical spectra
        double log_occu_c_alpha = XQuestScores::logOccupancyProb(theoretical_spec_linear_alpha, matched_spec_linear_alpha.size(), fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_);
        double log_occu_x_alpha = XQuestScores::logOccupancyProb(theoretical_spec_xlinks_alpha, matched_spec_xlinks_alpha.size(), fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_);
        double log_occu = 0;
        double log_occu_alpha = 0;
        double log_occu_beta = 0;

        if (type_is_cross_link)
        {
          double log_occu_c_beta = XQuestScores::logOccupancyProb(theoretical_spec_linear_beta, matched_spec_linear_beta.size(), fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_);
          double log_occu_x_beta = XQuestScores::logOccupancyProb(theoretical_spec_xlinks_beta, matched_spec_xlinks_beta.size(), fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_);
          log_occu = (log_occu_c_alpha + log_occu_x_alpha + log_occu_c_beta + log_occu_x_beta) / 4;
          log_occu_alpha = (log_occu_c_alpha + log_occu_x_alpha) / 2;
          log_occu_beta = (log_occu_c_beta + log_occu_x_beta) / 2;
        }
        else
        {
          log_occu = (log_occu_c_alpha + log_occu_x_alpha) / 2;
          log_occu_alpha = log_occu;
        }

        //Cross-correlation
        PeakSpectrum theoretical_spec_linear;
        PeakSpectrum theoretical_spec_xlinks;
        if (type_is_cross_link)
        {
          theoretical_spec_linear = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_linear_alpha, theoretical_spec_linear_beta);
          theoretical_spec_xlinks = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta);
        }
        else
        {
          theoretical_spec_linear = theoretical_spec_linear_alpha;
          theoretical_spec_xlinks = theoretical_spec_xlinks_alpha;
        }
        //
        PeakSpectrum theoretical_spec = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_linear, theoretical_spec_xlinks);
        double xcorrx_max = XQuestScores::xCorrelationPrescore(spectrum, theoretical_spec_xlinks, 0.1);
        double xcorrc_max = XQuestScores::xCorrelationPrescore(spectrum, theoretical_spec_linear, 0.1);

        // Compute score from the 4 scores and 4 weights
        // The weights are adapted from the xQuest algorithm (O. Rinner et al., 2008, "Identification of cross-linked peptides from large sequence databases"),
        // they were determined by an Linear Discriminant Analysis on CID fragmentation data.
        // The match-odds score does not work very well on HCD data and label-free cross-linkers (has the maximal possible value very often), so its weight was drastically reduced here.
        double xcorrx_weight = 2.488;
        double xcorrc_weight = 21.279;
        double match_odds_weight = 1.973;
        double wTIC_weight = 12.829;
        double intsum_weight = 1.8;

        double xquest_score = xcorrx_weight * xcorrx_max + xcorrc_weight * xcorrc_max + match_odds_weight * csm.match_odds + wTIC_weight * wTICold + intsum_weight * intsum;
        csm.xquest_score = xquest_score;

        // csm.precursor_correction = cross_link_candidate.precursor_correction;
        // double rel_error = OPXLHelper::computePrecursorError(csm, precursor_mz, precursor_charge);

        csm.percTIC = TIC;
        csm.wTIC = wTIC;
        csm.wTICold = wTICold;
        csm.int_sum = intsum;
        csm.intsum_alpha = intsum_alpha;
        csm.intsum_beta = intsum_beta;
        csm.total_current = total_current;
        // csm.precursor_error_ppm = rel_error;

        csm.log_occupancy = log_occu;
        csm.log_occupancy_alpha = log_occu_alpha;
        csm.log_occupancy_beta = log_occu_beta;

        csm.xcorrx_max = xcorrx_max;
        csm.xcorrc_max = xcorrc_max;

        csm.matched_linear_alpha = matched_spec_linear_alpha.size();
        csm.matched_linear_beta = matched_spec_linear_beta.size();
        csm.matched_xlink_alpha = matched_spec_xlinks_alpha.size();
        csm.matched_xlink_beta = matched_spec_xlinks_beta.size();
        csm.scan_index_light = scan_index;
        csm.scan_index_heavy = -1;

        csm.precursor_correction = cross_link_candidate.precursor_correction;

        if (!precursor_purities.empty())
        {
          csm.precursor_total_intensity = precursor_purities[spectrum.getNativeID()].total_intensity;
          csm.precursor_target_intensity = precursor_purities[spectrum.getNativeID()].target_intensity;
          csm.precursor_signal_proportion = precursor_purities[spectrum.getNativeID()].signal_proportion;
          csm.precursor_target_peak_count = precursor_purities[spectrum.getNativeID()].target_peak_count;
          csm.precursor_residual_peak_count = precursor_purities[spectrum.getNativeID()].interfering_peak_count;
        }

        // num_iso_peaks array from deisotoping
        if (deisotope)
        {
          DataArrays::IntegerDataArray num_iso_peaks_array;
          auto num_iso_peaks_array_it = getDataArrayByName(spectrum.getIntegerDataArrays(), "iso_peak_count");
          num_iso_peaks_array = *num_iso_peaks_array_it;

          OPXLHelper::isoPeakMeans(csm, num_iso_peaks_array, matched_spec_linear_alpha, matched_spec_linear_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta);
        }
        csm.ppm_error_abs_sum_linear_alpha = 0;
        csm.ppm_error_abs_sum_linear_beta = 0;
        csm.ppm_error_abs_sum_xlinks_alpha = 0;
        csm.ppm_error_abs_sum_xlinks_beta = 0;
        csm.ppm_error_abs_sum_linear = 0;
        csm.ppm_error_abs_sum_xlinks = 0;
        csm.ppm_error_abs_sum_alpha = 0;
        csm.ppm_error_abs_sum_beta = 0;
        csm.ppm_error_abs_sum = 0;

        // TODO find a better way to compute the absolute sum
        if (!ppm_error_array_linear_alpha.empty())
        {
          for (Size k = 0; k < ppm_error_array_linear_alpha.size(); ++k)
          {
            csm.ppm_error_abs_sum_linear_alpha += abs(ppm_error_array_linear_alpha[k]);
          }
          csm.ppm_error_abs_sum_linear_alpha = csm.ppm_error_abs_sum_linear_alpha / ppm_error_array_linear_alpha.size();
        }

        if (!ppm_error_array_linear_beta.empty())
        {
          for (Size k = 0; k < ppm_error_array_linear_beta.size(); ++k)
          {
            csm.ppm_error_abs_sum_linear_beta += abs(ppm_error_array_linear_beta[k]);
          }
          csm.ppm_error_abs_sum_linear_beta = csm.ppm_error_abs_sum_linear_beta / ppm_error_array_linear_beta.size();
        }

        if (!ppm_error_array_xlinks_alpha.empty())
        {
          for (Size k = 0; k < ppm_error_array_xlinks_alpha.size(); ++k)
          {
            csm.ppm_error_abs_sum_xlinks_alpha += abs(ppm_error_array_xlinks_alpha[k]);
          }
          csm.ppm_error_abs_sum_xlinks_alpha = csm.ppm_error_abs_sum_xlinks_alpha / ppm_error_array_xlinks_alpha.size();
        }

        if (!ppm_error_array_xlinks_beta.empty())
        {
          for (Size k = 0; k < ppm_error_array_xlinks_beta.size(); ++k)
          {
            csm.ppm_error_abs_sum_xlinks_beta += abs(ppm_error_array_xlinks_beta[k]);
          }
          csm.ppm_error_abs_sum_xlinks_beta = csm.ppm_error_abs_sum_xlinks_beta / ppm_error_array_xlinks_beta.size();
        }

        DataArrays::FloatDataArray ppm_error_array_linear;
        DataArrays::FloatDataArray ppm_error_array_xlinks;
        DataArrays::FloatDataArray ppm_error_array_alpha;
        DataArrays::FloatDataArray ppm_error_array_beta;
        DataArrays::FloatDataArray ppm_error_array;
        ppm_error_array_linear.insert(ppm_error_array_linear.end(), ppm_error_array_linear_alpha.begin(), ppm_error_array_linear_alpha.end());
        ppm_error_array_linear.insert(ppm_error_array_linear.end(), ppm_error_array_linear_beta.begin(), ppm_error_array_linear_beta.end());
        ppm_error_array_xlinks.insert(ppm_error_array_xlinks.end(), ppm_error_array_xlinks_alpha.begin(), ppm_error_array_xlinks_alpha.end());
        ppm_error_array_xlinks.insert(ppm_error_array_xlinks.end(), ppm_error_array_xlinks_beta.begin(), ppm_error_array_xlinks_beta.end());
        ppm_error_array_alpha.insert(ppm_error_array_alpha.end(), ppm_error_array_linear_alpha.begin(), ppm_error_array_linear_alpha.end());
        ppm_error_array_alpha.insert(ppm_error_array_alpha.end(), ppm_error_array_xlinks_alpha.begin(), ppm_error_array_xlinks_alpha.end());
        ppm_error_array_beta.insert(ppm_error_array_beta.end(), ppm_error_array_linear_beta.begin(), ppm_error_array_linear_beta.end());
        ppm_error_array_beta.insert(ppm_error_array_beta.end(), ppm_error_array_xlinks_beta.begin(), ppm_error_array_xlinks_beta.end());
        ppm_error_array.insert(ppm_error_array.end(), ppm_error_array_linear.begin(), ppm_error_array_linear.end());
        ppm_error_array.insert(ppm_error_array.end(), ppm_error_array_xlinks.begin(), ppm_error_array_xlinks.end());

        if (!ppm_error_array_linear.empty())
        {
          for (double ppm_error : ppm_error_array_linear)
          {
            csm.ppm_error_abs_sum_linear += abs(ppm_error);
          }
          csm.ppm_error_abs_sum_linear = csm.ppm_error_abs_sum_linear / ppm_error_array_linear.size();
        }

        if (!ppm_error_array_xlinks.empty())
        {
          for (double ppm_error : ppm_error_array_xlinks)
          {
            csm.ppm_error_abs_sum_xlinks += abs(ppm_error);
          }
          csm.ppm_error_abs_sum_xlinks = csm.ppm_error_abs_sum_xlinks / ppm_error_array_xlinks.size();
        }

        if (!ppm_error_array_alpha.empty())
        {
          for (double ppm_error : ppm_error_array_alpha)
          {
            csm.ppm_error_abs_sum_alpha += abs(ppm_error);
          }
          csm.ppm_error_abs_sum_alpha = csm.ppm_error_abs_sum_alpha / ppm_error_array_alpha.size();
        }

        if (!ppm_error_array_beta.empty())
        {
          for (double ppm_error : ppm_error_array_beta)
          {
            csm.ppm_error_abs_sum_beta += abs(ppm_error);
          }
          csm.ppm_error_abs_sum_beta = csm.ppm_error_abs_sum_beta / ppm_error_array_beta.size();
        }

        if (!ppm_error_array.empty())
        {
          for (double ppm_error : ppm_error_array)
          {
            csm.ppm_error_abs_sum += abs(ppm_error);
          }
          csm.ppm_error_abs_sum = csm.ppm_error_abs_sum / ppm_error_array.size();
        }

        // write fragment annotations
        vector<PeptideHit::PeakAnnotation> frag_annotations;

        OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_linear_alpha, theoretical_spec_linear_alpha, spectrum);
        OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_linear_beta, theoretical_spec_linear_beta, spectrum);
        OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, spectrum);
        OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, spectrum);

        // make annotations unique
        sort(frag_annotations.begin(), frag_annotations.end());
        vector<PeptideHit::PeakAnnotation>::iterator last_unique_anno = unique(frag_annotations.begin(), frag_annotations.end());
        if (last_unique_anno != frag_annotations.end())
        {
          frag_annotations.erase(last_unique_anno, frag_annotations.end());
        }

        csm.frag_annotations = frag_annotations;

#pragma omp critical (all_csms_spectrum_access)
        {
          all_csms_spectrum.push_back(csm);
        }
      } // end of parallel loop over top X candidates

      // collect top n matches to spectrum
      sort(all_csms_spectrum.rbegin(), all_csms_spectrum.rend(), OPXLDataStructs::CLSMScoreComparator());
      Size max_hit = min(all_csms_spectrum.size(), static_cast<Size>(number_top_hits_));

      for (Size top = 0; top < max_hit; top++)
      {
        all_csms_spectrum[top].rank = top+1;
        top_csms_spectrum.push_back(all_csms_spectrum[top]);
      }

      Size all_top_csms_current_index = 0;

#pragma omp critical (all_top_csms_access)
      {
        if (!top_csms_spectrum.empty())
        {
          all_top_csms.push_back(top_csms_spectrum);
          all_top_csms_current_index = all_top_csms.size()-1;
        }
      }

      // Write PeptideIdentifications and PeptideHits for n top hits
      if (!top_csms_spectrum.empty())
      {
        OPXLHelper::buildPeptideIDs(peptide_ids, top_csms_spectrum, all_top_csms, all_top_csms_current_index, spectra, scan_index, scan_index);
      }
#ifdef DEBUG_OPENPEPXLLFALGO
#pragma omp critical (LOG_DEBUG_access)
      OPENMS_LOG_DEBUG << "Next Spectrum ##################################" << endl;
#endif
    }

    std::cout << std::endl << "Total number of potential candidates searched: " << all_candidates_count << std::endl;

    // end of matching / scoring
    progresslogger.endProgress();

    // Add protein identifications
    PeptideIndexing pep_indexing;
    Param indexing_param = pep_indexing.getParameters();

    String d_prefix = decoy_prefix_ ? "prefix" : "suffix";
    indexing_param.setValue("decoy_string_position", d_prefix, "If set, protein accessions in the database contain 'decoy_string' as prefix.");
    indexing_param.setValue("decoy_string", decoy_string_, "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.");
    indexing_param.setValue("missing_decoy_action", "warn");
    indexing_param.setValue("enzyme:name", enzyme_name_);
    pep_indexing.setParameters(indexing_param);

    pep_indexing.run(fasta_db, protein_ids, peptide_ids);

    OPXLHelper::addProteinPositionMetaValues(peptide_ids);
    OPXLHelper::addBetaAccessions(peptide_ids);
    OPXLHelper::addXLTargetDecoyMV(peptide_ids);
    OPXLHelper::removeBetaPeptideHits(peptide_ids);
    OPXLHelper::computeDeltaScores(peptide_ids);
    OPXLHelper::addPercolatorFeatureList(protein_ids[0]);
    return OpenPepXLLFAlgorithm::ExitCodes::EXECUTION_OK;
  }
