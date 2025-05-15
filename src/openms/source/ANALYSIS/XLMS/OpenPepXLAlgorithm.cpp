// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Timo Sachsenberg, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OpenPepXLAlgorithm.h>

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
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
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>

#include <iostream>

using namespace std;
using namespace OpenMS;

// turn on additional debug output
// #define DEBUG_OPENPEPXLALGO

#ifdef _OPENMP
#include <omp.h>
#endif

  OpenPepXLAlgorithm::OpenPepXLAlgorithm()
    : DefaultParamHandler("OpenPepXLAlgorithm")
  {
    defaults_.setValue("decoy_string", "DECOY_", "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.");
    std::vector<std::string> bool_strings = {"true","false"};
    defaults_.setValue("decoy_prefix", "true", "Set to true, if the decoy_string is a prefix of accessions in the protein database. Otherwise it is a suffix.");
    defaults_.setValidStrings("decoy_prefix", bool_strings);

    defaults_.setValue("precursor:mass_tolerance", 10.0, "Width of precursor mass tolerance window");
    std::vector<std::string> mass_tolerance_unit_valid_strings = {"ppm","Da"};
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", mass_tolerance_unit_valid_strings);
    defaults_.setValue("precursor:min_charge", 2, "Minimum precursor charge to be considered.");
    defaults_.setValue("precursor:max_charge", 8, "Maximum precursor charge to be considered.");
    defaults_.setValue("precursor:corrections", ListUtils::create<int>("4, 3, 2, 1, 0"), "Monoisotopic peak correction. Matches candidates for possible monoisotopic precursor peaks for experimental mass m and given numbers n at masses (m - n * (C13-C12)). These should be ordered from more extreme to less extreme corrections. Numbers later in the list will be preferred in case of ambiguities.");
    defaults_.setSectionDescription("precursor", "Precursor filtering settings");

    defaults_.setValue("fragment:mass_tolerance", 20.0, "Fragment mass tolerance");
    defaults_.setValue("fragment:mass_tolerance_xlinks", 20.0, "Fragment mass tolerance for cross-link ions");
    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", mass_tolerance_unit_valid_strings);
    defaults_.setSectionDescription("fragment", "Fragment peak matching settings");

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
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

    defaults_.setValue("cross_linker:residue1", std::vector<std::string>{"K","N-term"}, "Comma separated residues, that the first side of a bifunctional cross-linker can attach to");
    defaults_.setValue("cross_linker:residue2", std::vector<std::string>{"K","N-term"}, "Comma separated residues, that the second side of a bifunctional cross-linker can attach to");
    defaults_.setValue("cross_linker:mass_light", 138.0680796, "Mass of the light cross-linker, linking two residues on one or two peptides");
    defaults_.setValue("cross_linker:mass_iso_shift", 12.075321, "Mass of the isotopic shift between the light and heavy linkers");
    defaults_.setValue("cross_linker:mass_mono_link", ListUtils::create<double>("156.07864431, 155.094628715"), "Possible masses of the linker, when attached to only one peptide");
    defaults_.setValue("cross_linker:name", "DSS", "Name of the searched cross-link, used to resolve ambiguity of equal masses (e.g. DSS or BS3)");
    defaults_.setSectionDescription("cross_linker", "Description of the cross-linker reagent");

    defaults_.setValue("algorithm:number_top_hits", 1, "Number of top hits reported for each spectrum pair");
    std::vector<std::string> deisotope_strings = {"true","false","auto"};
    defaults_.setValue("algorithm:deisotope", "auto", "Set to true, if the input spectra should be deisotoped before any other processing steps. If set to auto the spectra will be deisotoped, if the fragment mass tolerance is < 0.1 Da or < 100 ppm (0.1 Da at a mass of 1000)", {"advanced"});
    defaults_.setValidStrings("algorithm:deisotope", deisotope_strings);
    defaults_.setSectionDescription("algorithm", "Additional algorithm settings");

    defaults_.setValue("ions:b_ions", "true", "Search for peaks of b-ions.", {"advanced"});
    defaults_.setValue("ions:y_ions", "true", "Search for peaks of y-ions.", {"advanced"});
    defaults_.setValue("ions:a_ions", "false", "Search for peaks of a-ions.", {"advanced"});
    defaults_.setValue("ions:x_ions", "false", "Search for peaks of x-ions.", {"advanced"});
    defaults_.setValue("ions:c_ions", "false", "Search for peaks of c-ions.", {"advanced"});
    defaults_.setValue("ions:z_ions", "false", "Search for peaks of z-ions.", {"advanced"});
    defaults_.setValue("ions:neutral_losses", "true", "Search for neutral losses of H2O and H3N.", {"advanced"});
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

  OpenPepXLAlgorithm::~OpenPepXLAlgorithm() = default;

  void OpenPepXLAlgorithm::updateMembers_()
  {
    decoy_string_ = param_.getValue("decoy_string").toString();
    decoy_prefix_ = param_.getValue("decoy_prefix").toBool();

    min_precursor_charge_ = param_.getValue("precursor:min_charge");
    max_precursor_charge_ = param_.getValue("precursor:max_charge");
    precursor_mass_tolerance_ = param_.getValue("precursor:mass_tolerance");
    precursor_mass_tolerance_unit_ppm_ = (param_.getValue("precursor:mass_tolerance_unit") == "ppm");
    precursor_correction_steps_ = param_.getValue("precursor:corrections");

    fragment_mass_tolerance_ = param_.getValue("fragment:mass_tolerance");
    fragment_mass_tolerance_xlinks_ = param_.getValue("fragment:mass_tolerance_xlinks");
    fragment_mass_tolerance_unit_ppm_  = (param_.getValue("fragment:mass_tolerance_unit") == "ppm");

    cross_link_residue1_ = ListUtils::toStringList<std::string>(param_.getValue("cross_linker:residue1"));
    cross_link_residue2_ = ListUtils::toStringList<std::string>(param_.getValue("cross_linker:residue2"));
    cross_link_mass_light_ = param_.getValue("cross_linker:mass_light");
    cross_link_mass_iso_shift_ = param_.getValue("cross_linker:mass_iso_shift");
    cross_link_mass_mono_link_ = param_.getValue("cross_linker:mass_mono_link");
    cross_link_name_ = param_.getValue("cross_linker:name").toString();

    fixedModNames_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:fixed"));
    varModNames_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    max_variable_mods_per_peptide_ = static_cast<Size>(param_.getValue("modifications:variable_max_per_peptide"));
    peptide_min_size_ = static_cast<Size>(param_.getValue("peptide:min_size"));
    missed_cleavages_ = static_cast<Size>(param_.getValue("peptide:missed_cleavages"));
    enzyme_name_ = param_.getValue("peptide:enzyme").toString();

    number_top_hits_ = param_.getValue("algorithm:number_top_hits");
    deisotope_mode_ = param_.getValue("algorithm:deisotope").toString();

    add_y_ions_ = param_.getValue("ions:y_ions").toString();
    add_b_ions_ = param_.getValue("ions:b_ions").toString();
    add_x_ions_ = param_.getValue("ions:x_ions").toString();
    add_a_ions_ = param_.getValue("ions:a_ions").toString();
    add_c_ions_ = param_.getValue("ions:c_ions").toString();
    add_z_ions_ = param_.getValue("ions:z_ions").toString();
    add_losses_ = param_.getValue("ions:neutral_losses").toString();
  }

  OpenPepXLAlgorithm::ExitCodes OpenPepXLAlgorithm::run(PeakMap& unprocessed_spectra, ConsensusMap& cfeatures, std::vector<FASTAFile::FASTAEntry>& fasta_db, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, OPXLDataStructs::PreprocessedPairSpectra& preprocessed_pair_spectra, std::vector< std::pair<Size, Size> >& spectrum_pairs, std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, PeakMap& spectra)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(this->getLogType());

    // preprocess parameters for convenience
    if (fragment_mass_tolerance_xlinks_ < fragment_mass_tolerance_)
    {
      fragment_mass_tolerance_xlinks_ = fragment_mass_tolerance_;
    }
#ifdef DEBUG_OPENPEPXLALGO
    OPENMS_LOG_DEBUG << "XLinks Tolerance: " << fragment_mass_tolerance_xlinks_ << endl;
#endif

    std::sort(cross_link_mass_mono_link_.begin(), cross_link_mass_mono_link_.end(), std::greater< double >());

    // deisotope if "true" or if "auto" and the tolerance is below the threshold (0.1 Da or 100 ppm)
    bool deisotope = (deisotope_mode_ == "true") ||
                      (deisotope_mode_ == "auto" &&
                      ((!fragment_mass_tolerance_unit_ppm_ && fragment_mass_tolerance_ < 0.1) ||
                      (fragment_mass_tolerance_unit_ppm_ && fragment_mass_tolerance_ < 100)));

    set<String> fixed_unique(fixedModNames_.begin(), fixedModNames_.end());
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
    spectra = OPXLSpectrumProcessingAlgorithms::preprocessSpectra(picked_spectra, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, peptide_min_size_, min_precursor_charge_, max_precursor_charge_, deisotope, true);
    progresslogger.endProgress();
    picked_spectra.clear(true);

    // sort the spectra by RT, the order might have been changed by parallel preprocessing
    spectra.sortSpectra(false);

    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme_name_);
    digestor.setMissedCleavages(missed_cleavages_);

    IDMapper idmapper;
    Param p = idmapper.getParameters();
    p.setValue("rt_tolerance", 30.0);
    p.setValue("mz_tolerance", precursor_mass_tolerance_);
    String mz_measure = precursor_mass_tolerance_unit_ppm_ ? "ppm" : "Da";
    p.setValue("mz_measure", mz_measure);
    p.setValue("mz_reference", "precursor");
    p.setValue("ignore_charge", "false");
    idmapper.setParameters(p);

    progresslogger.startProgress(0, 1, "Map spectrum precursors to linked features...");
    idmapper.annotate(cfeatures, vector<PeptideIdentification>(), vector<ProteinIdentification>(), true, true, spectra);
    progresslogger.endProgress();

    vector< double > spectrum_precursors;

    // find pairs of MS2 spectra, that correspond to MS1 features linked by the consensus map / FeatureFinderMultiplex
    for (ConsensusMap::const_iterator cit = cfeatures.begin(); cit != cfeatures.end(); ++cit)
    {
      if (cit->getFeatures().size() == 2 && cit->getPeptideIdentifications().size() >= 2)
      {
        for (Size x = 0; x < cit->getPeptideIdentifications().size(); ++x)
        {
          if (static_cast<Size>(cit->getPeptideIdentifications()[x].getMetaValue("map_index")) == 0)
          {
            for (Size y = 0; y < cit->getPeptideIdentifications().size(); ++y)
            {
              if (static_cast<Size>(cit->getPeptideIdentifications()[y].getMetaValue("map_index")) == 1)
              {
                const PeptideIdentification& pi_0 = cit->getPeptideIdentifications()[x];
                const PeptideIdentification& pi_1 = cit->getPeptideIdentifications()[y];
                spectrum_pairs.emplace_back(pi_0.getMetaValue("spectrum_index"), pi_1.getMetaValue("spectrum_index"));
                double current_precursor_mz0 = spectra[pi_0.getMetaValue("spectrum_index")].getPrecursors()[0].getMZ();
                double current_precursor_mz1 = spectra[pi_1.getMetaValue("spectrum_index")].getPrecursors()[0].getMZ();
                double current_precursor_charge0 = spectra[pi_0.getMetaValue("spectrum_index")].getPrecursors()[0].getCharge();
                double current_precursor_charge1 = spectra[pi_1.getMetaValue("spectrum_index")].getPrecursors()[0].getCharge();

                double current_precursor_mass0 = (current_precursor_mz0 * current_precursor_charge0) - (current_precursor_charge0 * Constants::PROTON_MASS_U);
                double current_precursor_mass1 = (current_precursor_mz1 * current_precursor_charge1) - (current_precursor_charge1 * Constants::PROTON_MASS_U);
                spectrum_precursors.push_back(current_precursor_mass0);
                spectrum_precursors.push_back(current_precursor_mass1);
              }
            }
          }
        }
      }
    }
    sort(spectrum_precursors.begin(), spectrum_precursors.end());

    // create linear peak / shifted peak spectra for all pairs
    progresslogger.startProgress(0, 1, "Preprocessing Spectra Pairs...");
    preprocessed_pair_spectra = OpenPepXLAlgorithm::preprocessPairs_(spectra, spectrum_pairs, cross_link_mass_iso_shift_, fragment_mass_tolerance_, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, deisotope);
    progresslogger.endProgress();

    ProteinIdentification::SearchParameters search_params = protein_ids[0].getSearchParameters();
    String searched_charges((String(min_precursor_charge_)));
    for (int ch = min_precursor_charge_+1; ch <= max_precursor_charge_; ++ch)
    {
      searched_charges += "," + String(ch);
    }
    search_params.charges = searched_charges;
    search_params.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_name_));
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
    search_params.setMetaValue("cross_link:mass", cross_link_mass_light_);
    search_params.setMetaValue("cross_link:mass_isoshift", cross_link_mass_iso_shift_);
    search_params.setMetaValue("cross_link:mass_monolink", cross_link_mass_mono_link_);
    search_params.setMetaValue("cross_link:name", cross_link_name_);
    search_params.setMetaValue("precursor:corrections", precursor_correction_steps_);

    search_params.setMetaValue("modifications:variable_max_per_peptide", max_variable_mods_per_peptide_);
    protein_ids[0].setSearchParameters(search_params);
    protein_ids[0].setScoreType("OpenPepXL_Protein_Score");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    vector<OPXLDataStructs::AASeqWithMass> peptide_masses;
    peptide_masses = OPXLHelper::digestDatabase(fasta_db, digestor, peptide_min_size_, cross_link_residue1_, cross_link_residue2_, fixed_modifications,  variable_modifications, max_variable_mods_per_peptide_);

    // create spectrum generator
    TheoreticalSpectrumGeneratorXLMS specGen;
    SimpleTSGXLMS specGen_mainscore;

    // Set parameters for cross-link fragmentation
    Param specGenParams = specGen.getParameters();
    specGenParams.setValue("add_y_ions", add_y_ions_, "Add peaks of b-ions to the spectrum");
    specGenParams.setValue("add_b_ions", add_b_ions_, "Add peaks of y-ions to the spectrum");
    specGenParams.setValue("add_a_ions", add_a_ions_, "Add peaks of a-ions to the spectrum");
    specGenParams.setValue("add_x_ions", add_x_ions_, "Add peaks of c-ions to the spectrum");
    specGenParams.setValue("add_c_ions", add_c_ions_, "Add peaks of x-ions to the spectrum");
    specGenParams.setValue("add_z_ions", add_z_ions_, "Add peaks of z-ions to the spectrum");
    specGenParams.setValue("add_losses", add_losses_, "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    specGenParams.setValue("add_metainfo", "true");
    specGenParams.setValue("add_isotopes", "true", "If set to 1 isotope peaks of the product ion peaks are added");
    specGenParams.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");
    specGenParams.setValue("add_precursor_peaks", "true", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    specGenParams.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    specGenParams.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    specGenParams.setValue("add_k_linked_ions", "false");
    specGen.setParameters(specGenParams);

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

#ifdef DEBUG_OPENPEPXLALGO
    OPENMS_LOG_DEBUG << "Peptide candidates: " << peptide_masses.size() << endl;
#endif

    search_params = protein_ids[0].getSearchParameters();
    search_params.setMetaValue("MS:1001029", peptide_masses.size()); // number of sequences searched = MS:1001029
    protein_ids[0].setSearchParameters(search_params);

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
      max_peptide_allowed_error = precursor_mass_tolerance_;
    }

    // maximal possible peptide mass given the largest precursor
    double max_peptide_mass = max_precursor_mass - cross_link_mass_light_ + max_peptide_allowed_error;

    // search for the first mass greater than the maximum, use everything before that peptide
    vector<OPXLDataStructs::AASeqWithMass>::iterator last = upper_bound(peptide_masses.begin(), peptide_masses.end(), max_peptide_mass, OPXLDataStructs::AASeqWithMassComparator());
    vector<OPXLDataStructs::AASeqWithMass> filtered_peptide_masses;
    filtered_peptide_masses.assign(peptide_masses.begin(), last);
    peptide_masses.clear();

    // iterate over all spectra
    progresslogger.startProgress(0, 1, "Matching to theoretical spectra and scoring...");
    Size spectrum_counter = 0;

    for (SignedSize pair_index = 0; pair_index < static_cast<SignedSize>(spectrum_pairs.size()); ++pair_index)
    {
      Size scan_index = spectrum_pairs[pair_index].first;
      Size scan_index_heavy = spectrum_pairs[pair_index].second;

#ifdef DEBUG_OPENPEPXLALGO
      OPENMS_LOG_DEBUG << "New scan indices: " << scan_index << "\t" << scan_index_heavy << endl;
#endif
      const PeakSpectrum& spectrum_light = spectra[scan_index];
      const double precursor_charge = spectrum_light.getPrecursors()[0].getCharge();
      const double precursor_mz = spectrum_light.getPrecursors()[0].getMZ();
      const double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

      const PeakSpectrum& linear_peaks = preprocessed_pair_spectra.spectra_linear_peaks[pair_index];
      const PeakSpectrum& xlink_peaks = preprocessed_pair_spectra.spectra_xlink_peaks[pair_index];
      const PeakSpectrum& all_peaks = preprocessed_pair_spectra.spectra_all_peaks[pair_index];

      vector< OPXLDataStructs::CrossLinkSpectrumMatch > top_csms_spectrum;

      // ignore this spectrum pair, if they have less paired peaks than the minimal peptide size
      if (all_peaks.size() < peptide_min_size_)
      {
        continue;
      }

      vector <OPXLDataStructs::ProteinProteinCrossLink> cross_link_candidates = OPXLHelper::collectPrecursorCandidates(precursor_correction_steps_, precursor_mass, precursor_mass_tolerance_, precursor_mass_tolerance_unit_ppm_, filtered_peptide_masses, cross_link_mass_light_, cross_link_mass_mono_link_, cross_link_residue1_, cross_link_residue2_, cross_link_name_);

      spectrum_counter++;
      cout << "Processing spectrum pair " << spectrum_counter << " / " << spectrum_pairs.size() << endl;
      cout << "Light Spectrum ID: " << spectrum_light.getNativeID() << " |\tHeavy Spectrum ID: " << spectra[scan_index_heavy].getNativeID() << "\t| at: " << DateTime::now().getTime() << endl;
      cout << "Number of peaks in light spectrum: " << spectrum_light.size() << " |\tNumber of candidates: " << cross_link_candidates.size() << endl;

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
        if (theoretical_spec_linear_alpha.empty())
        {
          continue;
        }

        vector< pair< Size, Size > > matched_spec_linear_alpha;
        vector< pair< Size, Size > > matched_spec_linear_beta;
        vector< pair< Size, Size > > matched_spec_xlinks_alpha;
        vector< pair< Size, Size > > matched_spec_xlinks_beta;

        if (!linear_peaks.empty())
        {
          DataArrays::IntegerDataArray exp_charges;
          if (!linear_peaks.getIntegerDataArrays().empty())
          {
            exp_charges = linear_peaks.getIntegerDataArrays()[0];
          }
          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(matched_spec_linear_alpha, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_linear_alpha, linear_peaks, exp_charges);
          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(matched_spec_linear_beta, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_linear_beta, linear_peaks, exp_charges);
        }
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

        if (!xlink_peaks.empty())
        {
          DataArrays::IntegerDataArray exp_charges;
          if (!xlink_peaks.getIntegerDataArrays().empty())
          {
            exp_charges = xlink_peaks.getIntegerDataArrays()[0];
          }
          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(matched_spec_xlinks_alpha, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_xlinks_alpha, xlink_peaks, exp_charges);
          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentSimple(matched_spec_xlinks_beta, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_xlinks_beta, xlink_peaks, exp_charges);
        }

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
      // progresslogger.endProgress();
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

#ifdef DEBUG_OPENPEPXLALGO
        double candidate_mz = (alpha.getMonoWeight() + beta.getMonoWeight() +  cross_link_candidate.cross_linker_mass+ (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / precursor_charge;
#pragma omp critical (LOG_DEBUG_access)
        {
          OPENMS_LOG_DEBUG << "Pair: " << alpha.toString() << "-" << beta.toString() << " matched to light spectrum " << scan_index << "\t and heavy spectrum " << scan_index_heavy
              << " with m/z: " << precursor_mz << "\t" << "and candidate m/z: " << candidate_mz << "\tK Positions: " << cross_link_candidate.cross_link_position.first << "\t" << cross_link_candidate.cross_link_position.second << endl;
        }
#endif

        OPXLDataStructs::CrossLinkSpectrumMatch csm = mainscore_csms_spectrum[i];
        csm.cross_link = cross_link_candidate;

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

        specGen.getLinearIonSpectrum(theoretical_spec_linear_alpha, alpha, cross_link_candidate.cross_link_position.first, true, 2, link_pos_B);
        if (type_is_cross_link)
        {
          theoretical_spec_linear_beta.reserve(1500);
          theoretical_spec_xlinks_beta.reserve(1500);
          specGen.getLinearIonSpectrum(theoretical_spec_linear_beta, beta, cross_link_candidate.cross_link_position.second, false, 2);
          specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate, true, 1, precursor_charge);
          specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_beta, cross_link_candidate, false, 1, precursor_charge);
        }
        else
        {
          // Function for mono-links or loop-links
          specGen.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, alpha, cross_link_candidate.cross_link_position.first, precursor_mass, true, 2, precursor_charge, link_pos_B);
        }

        vector< pair< Size, Size > > matched_spec_linear_alpha;
        vector< pair< Size, Size > > matched_spec_linear_beta;
        vector< pair< Size, Size > > matched_spec_xlinks_alpha;
        vector< pair< Size, Size > > matched_spec_xlinks_beta;

        DataArrays::FloatDataArray ppm_error_array_linear_alpha;
        DataArrays::FloatDataArray ppm_error_array_xlinks_alpha;
        DataArrays::FloatDataArray ppm_error_array_linear_beta;
        DataArrays::FloatDataArray ppm_error_array_xlinks_beta;

        if (!linear_peaks.empty())
        {
          DataArrays::IntegerDataArray theo_charges_alpha;
          DataArrays::IntegerDataArray theo_charges_beta;
          DataArrays::IntegerDataArray exp_charges;

          auto theo_alpha_it = getDataArrayByName(theoretical_spec_linear_alpha.getIntegerDataArrays(), "charge");
          theo_charges_alpha = *theo_alpha_it;
          if (!theoretical_spec_linear_beta.empty())
          {
            auto theo_beta_it = getDataArrayByName(theoretical_spec_linear_beta.getIntegerDataArrays(), "charge");
            theo_charges_beta = *theo_beta_it;
          }

          auto exp_it = getDataArrayByName(linear_peaks.getIntegerDataArrays(), "charge");
          if (exp_it != linear_peaks.getIntegerDataArrays().end())
          {
            if (!exp_it->empty())
            {
              exp_charges = *exp_it;
            }
          }

          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_linear_alpha, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_linear_alpha, linear_peaks, theo_charges_alpha, exp_charges, ppm_error_array_linear_alpha);
          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_linear_beta, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_linear_beta, linear_peaks, theo_charges_beta, exp_charges, ppm_error_array_linear_beta);
        }
        if (!xlink_peaks.empty())
        {
          DataArrays::IntegerDataArray theo_charges_alpha;
          DataArrays::IntegerDataArray theo_charges_beta;
          DataArrays::IntegerDataArray exp_charges;

          auto theo_alpha_it = getDataArrayByName(theoretical_spec_xlinks_alpha.getIntegerDataArrays(), "charge");
          theo_charges_alpha = *theo_alpha_it;
          if (!theoretical_spec_xlinks_beta.empty())
          {
            auto theo_beta_it = getDataArrayByName(theoretical_spec_xlinks_beta.getIntegerDataArrays(), "charge");
            theo_charges_beta = *theo_beta_it;
          }

          auto exp_it = getDataArrayByName(xlink_peaks.getIntegerDataArrays(), "charge");
          exp_charges = *exp_it;

          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_xlinks_alpha, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_xlinks_alpha, xlink_peaks, theo_charges_alpha, exp_charges, ppm_error_array_xlinks_alpha);
          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_xlinks_beta, fragment_mass_tolerance_xlinks_, fragment_mass_tolerance_unit_ppm_, theoretical_spec_xlinks_beta, xlink_peaks, theo_charges_beta, exp_charges, ppm_error_array_xlinks_beta);
        }

        // Pre-Score calculations
        Size matched_alpha_count = matched_spec_linear_alpha.size() + matched_spec_xlinks_alpha.size();
        Size theor_alpha_count = theoretical_spec_linear_alpha.size() + theoretical_spec_xlinks_alpha.size();
        Size matched_beta_count = matched_spec_linear_beta.size() + matched_spec_xlinks_beta.size();
        Size theor_beta_count = theoretical_spec_linear_beta.size() + theoretical_spec_xlinks_beta.size();

#ifdef DEBUG_OPENPEPXLALGO
#pragma omp critical (LOG_DEBUG_access)
        {
          OPENMS_LOG_DEBUG << "matched peaks: " << matched_alpha_count + matched_beta_count << endl;
          OPENMS_LOG_DEBUG << "theoretical peaks: " << theor_alpha_count + theor_beta_count << endl;
          OPENMS_LOG_DEBUG << "exp peaks: " << all_peaks.size() << endl;
        }
#endif

        if (matched_alpha_count + matched_beta_count > 0)
        {
          // Simplified pre-Score
          double pre_score = 0;
          if (type_is_cross_link)
          {
            pre_score = XQuestScores::preScore(matched_alpha_count, theor_alpha_count, matched_beta_count, theor_beta_count);
          }
          else
          {
            pre_score = XQuestScores::preScore(matched_alpha_count, theor_alpha_count);
          }

          // compute intsum score
          double intsum = XQuestScores::totalMatchedCurrent(matched_spec_linear_alpha, matched_spec_linear_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta, linear_peaks, xlink_peaks);


          // Total ion intensity of light spectrum
          // sum over linear and xlink ion spectra instead of unfiltered
          double total_current = 0;
          for (SignedSize j = 0; j < static_cast<SignedSize>(linear_peaks.size()); ++j)
          {
            total_current += linear_peaks[j].getIntensity();
          }
          for (SignedSize j = 0; j < static_cast<SignedSize>(xlink_peaks.size()); ++j)
          {
            total_current += xlink_peaks[j].getIntensity();
          }
          double TIC = intsum / total_current;

          // TIC_alpha and _beta
          double intsum_alpha = XQuestScores::matchedCurrentChain(matched_spec_linear_alpha, matched_spec_xlinks_alpha, linear_peaks, xlink_peaks);
          double intsum_beta = 0;
          if (type_is_cross_link)
          {
            intsum_beta = XQuestScores::matchedCurrentChain(matched_spec_linear_beta, matched_spec_xlinks_beta, linear_peaks, xlink_peaks);
          }

          // normalize TIC_alpha and  _beta
          if ((intsum_alpha + intsum_beta) > 0.0)
          {
            intsum_alpha = intsum_alpha * intsum / (intsum_alpha + intsum_beta);
            intsum_beta = intsum_beta *  intsum / (intsum_alpha + intsum_beta);
          }

          // compute wTIC
          double wTIC = XQuestScores::weightedTICScore(alpha.size(), beta.size(), intsum_alpha, intsum_beta, total_current, type_is_cross_link);
          double wTICold = XQuestScores::weightedTICScoreXQuest(alpha.size(), beta.size(), intsum_alpha, intsum_beta, total_current, type_is_cross_link);

          // compute match odds (unweighted), the 3 is the number of charge states in the theoretical spectra
          double log_occu_c_alpha = XQuestScores::logOccupancyProb(theoretical_spec_linear_alpha, matched_spec_linear_alpha.size(), fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm_);
          double log_occu_x_alpha = XQuestScores::logOccupancyProb(theoretical_spec_xlinks_alpha, matched_spec_xlinks_alpha.size(), fragment_mass_tolerance_xlinks_ , fragment_mass_tolerance_unit_ppm_);
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

          PeakSpectrum theoretical_spec = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_linear, theoretical_spec_xlinks);
          PeakSpectrum theoretical_spec_alpha = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_linear_alpha, theoretical_spec_xlinks_alpha);

          PeakSpectrum theoretical_spec_beta;
          if (type_is_cross_link)
          {
            theoretical_spec_beta = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_linear_beta, theoretical_spec_xlinks_beta);
          }

          double xcorrx_max = XQuestScores::xCorrelationPrescore(xlink_peaks, theoretical_spec_xlinks, 0.1);
          double xcorrc_max = XQuestScores::xCorrelationPrescore(linear_peaks, theoretical_spec_linear, 0.1);

          // Compute score from the 4 scores and 4 weights
          // The weights are adapted from the xQuest algorithm (O. Rinner et al., 2008, "Identification of cross-linked peptides from large sequence databases"),
          // they were determined by an Linear Discriminant Analysis on CID fragmentation data.
          double xcorrx_weight = 2.488;
          double xcorrc_weight = 21.279;
          double match_odds_weight = 1.973;
          double wTIC_weight = 12.829;
          double intsum_weight = 1.8;

          double xquest_score = xcorrx_weight * xcorrx_max + xcorrc_weight * xcorrc_max + match_odds_weight * csm.match_odds + wTIC_weight * wTICold + intsum_weight * intsum;
          csm.xquest_score = xquest_score;

          csm.pre_score = pre_score;
          csm.percTIC = TIC;
          csm.wTIC = wTIC;
          csm.wTICold = wTICold;
          csm.int_sum = intsum;
          csm.intsum_alpha = intsum_alpha;
          csm.intsum_beta = intsum_beta;
          csm.total_current = total_current;
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
          csm.scan_index_heavy = scan_index_heavy;

          if (precursor_purities.size() > scan_index)
          {
            csm.precursor_total_intensity = precursor_purities[spectrum_light.getNativeID()].total_intensity;
            csm.precursor_target_intensity = precursor_purities[spectrum_light.getNativeID()].target_intensity;
            csm.precursor_signal_proportion = precursor_purities[spectrum_light.getNativeID()].signal_proportion;
            csm.precursor_target_peak_count = precursor_purities[spectrum_light.getNativeID()].target_peak_count;
            csm.precursor_residual_peak_count = precursor_purities[spectrum_light.getNativeID()].interfering_peak_count;
          }

          // num_iso_peaks array from deisotoping
          if (deisotope)
          {
            auto num_iso_peaks_array_it = getDataArrayByName(all_peaks.getIntegerDataArrays(), "iso_peak_count");
            DataArrays::IntegerDataArray num_iso_peaks_array = *num_iso_peaks_array_it;
            auto num_iso_peaks_array_linear_it = getDataArrayByName(linear_peaks.getIntegerDataArrays(), "iso_peak_count");
            DataArrays::IntegerDataArray num_iso_peaks_array_linear = *num_iso_peaks_array_linear_it;
            auto num_iso_peaks_array_xlinks_it = getDataArrayByName(xlink_peaks.getIntegerDataArrays(), "iso_peak_count");
            DataArrays::IntegerDataArray num_iso_peaks_array_xlinks = *num_iso_peaks_array_xlinks_it;

            csm.num_iso_peaks_mean = Math::mean(num_iso_peaks_array.begin(), num_iso_peaks_array.end());

            vector< double > iso_peaks_linear_alpha;
            vector< double > iso_peaks_linear_beta;
            vector< double > iso_peaks_xlinks_alpha;
            vector< double > iso_peaks_xlinks_beta;

            if (!matched_spec_linear_alpha.empty())
            {
              for (const auto& match : matched_spec_linear_alpha)
              {
                iso_peaks_linear_alpha.push_back(num_iso_peaks_array_linear[match.second]);
              }
              csm.num_iso_peaks_mean_linear_alpha = Math::mean(iso_peaks_linear_alpha.begin(), iso_peaks_linear_alpha.end());
            }

            if (!matched_spec_linear_beta.empty())
            {
              for (const auto& match : matched_spec_linear_beta)
              {
                iso_peaks_linear_beta.push_back(num_iso_peaks_array_linear[match.second]);
              }
              csm.num_iso_peaks_mean_linear_beta = Math::mean(iso_peaks_linear_beta.begin(), iso_peaks_linear_beta.end());
            }

            if (!matched_spec_xlinks_alpha.empty())
            {
              for (const auto& match : matched_spec_xlinks_alpha)
              {
                iso_peaks_xlinks_alpha.push_back(num_iso_peaks_array_xlinks[match.second]);
              }
              if (!iso_peaks_xlinks_alpha.empty())
              {
                csm.num_iso_peaks_mean_xlinks_alpha = Math::mean(iso_peaks_xlinks_alpha.begin(), iso_peaks_xlinks_alpha.end());
              }
            }

            if (!matched_spec_xlinks_beta.empty())
            {
              for (const auto& match : matched_spec_xlinks_beta)
              {
                iso_peaks_xlinks_beta.push_back(num_iso_peaks_array_xlinks[match.second]);
              }
              if (!iso_peaks_xlinks_beta.empty())
              {
                csm.num_iso_peaks_mean_xlinks_beta = Math::mean(iso_peaks_xlinks_beta.begin(), iso_peaks_xlinks_beta.end());
              }
            }
          }

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
          OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_linear_alpha, theoretical_spec_linear_alpha, linear_peaks);
          OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_linear_beta, theoretical_spec_linear_beta, linear_peaks);
          OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, xlink_peaks);
          OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, xlink_peaks);

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

      // Write PeptideIdentifications and PeptideHits for n top hits of this spectrum
      if (!top_csms_spectrum.empty())
      {
        OPXLHelper::buildPeptideIDs(peptide_ids, top_csms_spectrum, all_top_csms, all_top_csms_current_index, spectra, scan_index, scan_index_heavy);
      }

#ifdef DEBUG_OPENPEPXLALGO
#pragma omp critical (LOG_DEBUG_access)
      OPENMS_LOG_DEBUG << "Next Spectrum #############################################" << endl;
#endif
    } // end of matching / scoring, end of parallel for-loop

    progresslogger.endProgress();

    peptide_ids = OPXLHelper::combineTopRanksFromPairs(peptide_ids, number_top_hits_);

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
    return OpenPepXLAlgorithm::ExitCodes::EXECUTION_OK;
  }

  // create linear / shifted peak spectra for all pairs
  OPXLDataStructs::PreprocessedPairSpectra OpenPepXLAlgorithm::preprocessPairs_(const PeakMap& spectra, const vector< pair<Size, Size> >& spectrum_pairs, const double cross_link_mass_iso_shift, double fragment_mass_tolerance, double fragment_mass_tolerance_xlinks, bool fragment_mass_tolerance_unit_ppm, bool deisotope)
  {
    OPXLDataStructs::PreprocessedPairSpectra preprocessed_pair_spectra(spectrum_pairs.size());

#pragma omp parallel for
    for (SignedSize pair_index = 0; pair_index < static_cast<SignedSize>(spectrum_pairs.size()); ++pair_index)
    {
      Size scan_index = spectrum_pairs[pair_index].first;
      const PeakSpectrum& spectrum_light = spectra[scan_index];
      const Size scan_index_heavy = spectrum_pairs[pair_index].second;
      Size max_charge_xlink = spectrum_light.getPrecursors()[0].getCharge();

      const PeakSpectrum& spectrum_heavy = spectra[scan_index_heavy];
      vector< pair< Size, Size > > matched_fragments_without_shift;
      DataArrays::FloatDataArray dummy_array;
      DataArrays::IntegerDataArray dummy_charges;
      OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_fragments_without_shift, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, spectrum_light, spectrum_heavy, dummy_charges, dummy_charges, dummy_array, 0.3);

      // transform by m/z difference between unlabeled and labeled cross-link to make heavy and light comparable.
      PeakSpectrum xlink_peaks;
      PeakSpectrum::IntegerDataArray spectrum_heavy_charges;
      PeakSpectrum::IntegerDataArray spectrum_light_iso_peaks;

      auto spectrum_heavy_charges_it = getDataArrayByName(spectrum_heavy.getIntegerDataArrays(), "charge");
      if (spectrum_heavy_charges_it != spectrum_heavy.getIntegerDataArrays().end())
      {
        if (!spectrum_heavy_charges_it->empty())
        {
          spectrum_heavy_charges = *spectrum_heavy_charges_it;
        }
      }
      auto spectrum_light_iso_peaks_it = getDataArrayByName(spectrum_light.getIntegerDataArrays(), "iso_peak_count");
      if (spectrum_light_iso_peaks_it != spectrum_light.getIntegerDataArrays().end())
      {
        if (!spectrum_light_iso_peaks_it->empty())
        {
          spectrum_light_iso_peaks = *spectrum_light_iso_peaks_it;
        }
      }

      if (deisotope)
      {
        xlink_peaks.getIntegerDataArrays().resize(2);
        xlink_peaks.getIntegerDataArrays()[0].setName("charge");
        xlink_peaks.getIntegerDataArrays()[1].setName("iso_peak_count");
      }
      else
      {
        xlink_peaks.getIntegerDataArrays().resize(1);
        xlink_peaks.getIntegerDataArrays()[0].setName("charge");
      }

      // keep track of matched peaks
      vector<Size> used_peaks;

      // transform all peaks in the heavy spectrum by shifting them, considering all expected charge states
      for (Size charge = 1; charge <= max_charge_xlink; ++charge)
      {
        PeakSpectrum spectrum_heavy_to_light;
        PeakSpectrum::IntegerDataArray spectrum_heavy_to_light_charges;
        spectrum_heavy_to_light_charges.setName("charge");
        double mass_shift = cross_link_mass_iso_shift / charge;

        // transform heavy spectrum
        for (Size i = 0; i != spectrum_heavy.size(); ++i)
        {
          bool charge_fits = true;
          // check if the charge for the heavy peak determined by deisotoping matches the currently considered charge
          if (deisotope && spectrum_heavy_charges[i] != 0 && static_cast<unsigned int>(spectrum_heavy_charges[i]) != charge)
          {
            charge_fits = false;
          }

          if (charge_fits)
          {
            Peak1D p = spectrum_heavy[i];
            p.setMZ(p.getMZ() - mass_shift);
            spectrum_heavy_to_light.push_back(p);
            spectrum_heavy_to_light_charges.push_back(charge);
          }
        }
        spectrum_heavy_to_light.getIntegerDataArrays().push_back(spectrum_heavy_to_light_charges);

        // align peaks from light spectrum with shifted peaks from heavy spectrum
        // matching fragments are potentially carrying the cross-linker
        vector< pair< Size, Size > > matched_fragments_with_shift;

        spectrum_heavy_to_light.sortByPosition();
        if (!spectrum_heavy_to_light.empty())
        {
          dummy_array.clear();
          OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_fragments_with_shift, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, spectrum_light, spectrum_heavy_to_light, dummy_charges, dummy_charges, dummy_array, 0.3);

          // fill xlink_peaks spectrum with matched peaks from the light spectrum and add the currently considered charge
          for (Size i = 0; i < matched_fragments_with_shift.size(); ++i)
          {
            // test whether this peak was matched with a lower charge before (biased towards lower charge matches, if one light peak matches to multiple heavy peaks with different charges)
            vector<Size>::iterator it = find(used_peaks.begin(), used_peaks.end(), matched_fragments_with_shift[i].first);
            if (it == used_peaks.end())
            {
              xlink_peaks.push_back(spectrum_light[matched_fragments_with_shift[i].first]);
              xlink_peaks.getIntegerDataArrays()[0].push_back(charge);
              used_peaks.push_back(matched_fragments_with_shift[i].first);
              if (deisotope)
              {
                xlink_peaks.getIntegerDataArrays()[1].push_back(spectrum_light_iso_peaks[matched_fragments_with_shift[i].first]);
              }
            }
          }
        }
      }

      // generate linear peaks spectrum, include charges determined through deisotoping in preprocessing
      PeakSpectrum linear_peaks;
      PeakSpectrum::IntegerDataArray spectrum_light_charges;

      auto spectrum_light_charges_it = getDataArrayByName(spectrum_light.getIntegerDataArrays(), "charge");
      if (spectrum_light_charges_it != spectrum_light.getIntegerDataArrays().end())
      {
        if (!spectrum_light_charges_it->empty())
        {
          spectrum_light_charges = *spectrum_light_charges_it;
          linear_peaks.getIntegerDataArrays().resize(2);
          linear_peaks.getIntegerDataArrays()[0].setName("charge");
          linear_peaks.getIntegerDataArrays()[1].setName("iso_peak_count");
        }
      }

      for (Size i = 0; i != matched_fragments_without_shift.size(); ++i)
      {
        linear_peaks.push_back(spectrum_light[matched_fragments_without_shift[i].first]);
        if (!spectrum_light_charges.empty())
        {
          linear_peaks.getIntegerDataArrays()[0].push_back(spectrum_light_charges[matched_fragments_without_shift[i].first]);
          linear_peaks.getIntegerDataArrays()[1].push_back(spectrum_light_iso_peaks[matched_fragments_without_shift[i].first]);
        }
      }

      // TODO replace with window mower
      Size max_peak_number = 250;
      NLargest nfilter(max_peak_number);
      nfilter.filterSpectrum(linear_peaks);
      nfilter.filterSpectrum(xlink_peaks);

      PeakSpectrum all_peaks = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(linear_peaks, xlink_peaks);

      linear_peaks.setPrecursors(spectrum_light.getPrecursors());
      xlink_peaks.setPrecursors(spectrum_light.getPrecursors());
      all_peaks.setPrecursors(spectrum_light.getPrecursors());

      linear_peaks.sortByPosition();
      xlink_peaks.sortByPosition();
      all_peaks.sortByPosition();

  #pragma omp critical (preprocessed_pair_spectra_access)
      {
        swap(preprocessed_pair_spectra.spectra_linear_peaks[pair_index], linear_peaks);
        swap(preprocessed_pair_spectra.spectra_xlink_peaks[pair_index], xlink_peaks);
        swap(preprocessed_pair_spectra.spectra_all_peaks[pair_index], all_peaks);
      }
    }
    return preprocessed_pair_spectra;
  }
