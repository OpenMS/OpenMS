// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann, Leon Bichmann $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/APPLICATIONS/SearchEngineBase.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/SYSTEM/JavaInfo.h>

#include <QStringList>

#include <iostream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MSFraggerAdapter MSFraggerAdapter

@brief Peptide Identification with MSFragger

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; MSFraggerAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

@em MSFragger must be installed before this adapter can be used. This adapter is fully compatible with version 3.2 of MSFragger
and later versions of MSFragger were tested up to version 3.5.

All MSFragger parameters (as specified in the fragger.params file) have been transcribed to parameters of this OpenMS util.
It is not possible to provide an explicit fragger.params file to avoid redundancy with the ini file.
This adapter creates an fragger.params file prior to calling MSFragger. If the fragger.params file should be inspected, set the
-debug option to 2. MSFraggerAdapter will print the path to the working directory to standard out.

MSFragger can process multiple input files (mzML, mzXML) one after another. The number of output files specified must match
the number of input spectra files. The output file is then matched to the input file by index. The default parameters of the
adapter are the same as given by the official MSFragger manual.

Please cite:
Andy T Kong, Felipe V Leprevost, Dmitry M Avtonomov, Dattatreya Mellacheruvu & Alexey I Nesvizhskii
MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry–based proteomics
Nature Methods volume 14, pages 513–520 (2017) doi:10.1038/nmeth.4256


<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MSFraggerAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MSFraggerAdapter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPMSFraggerAdapter final :
  public SearchEngineBase
{
public:

  static const String license;

  static const String java_executable;
  static const String java_heapmemory;
  static const String executable;
  static const String in;
  static const String out;
  static const String opt_out;
  static const String database;

  // tolerance
  static const String precursor_mass_tolerance_lower;
  static const String precursor_mass_tolerance_upper;
  static const String precursor_mass_unit;
  static const String precursor_true_tolerance;
  static const String precursor_true_unit;
  static const String fragment_mass_tolerance;
  static const String fragment_mass_unit;
  static const String isotope_error;

  // digest
  static const String search_enzyme_name;
  static const String search_enzyme_cutafter;
  static const String search_enzyme_nocutbefore;
  static const String num_enzyme_termini;
  static const String allowed_missed_cleavage;
  static const String digest_min_length;
  static const String digest_max_length;
  static const String digest_mass_range_min;
  static const String digest_mass_range_max;

  // varmod
  static const String clip_nterm_m;
  static const String varmod_masses;
  static const String varmod_syntax;
  static const String varmod_enable_common;
  static const String variable_modifications_unimod;
  static const String not_allow_multiple_variable_mods_on_residue;
  static const String max_variable_mods_per_peptide;
  static const String max_variable_mods_combinations;

  // spectrum
  static const String minimum_peaks;
  static const String use_topn_peaks;
  static const String minimum_ratio;
  static const String clear_mz_range_min;
  static const String clear_mz_range_max;
  static const String max_fragment_charge;
  static const String override_charge;
  static const String precursor_charge_min;
  static const String precursor_charge_max;

  // search
  static const String track_zero_topn;
  static const String zero_bin_accept_expect;
  static const String zero_bin_mult_expect;
  static const String add_topn_complementary;
  static const String min_fragments_modeling;
  static const String min_matched_fragments;
  static const String output_report_topn;
  static const String output_max_expect;
  static const String localize_delta_mass;

  // statmod
  static const String add_cterm_peptide;
  static const String add_nterm_peptide;
  static const String add_cterm_protein;
  static const String add_nterm_protein;
  static const String add_G_glycine;
  static const String add_A_alanine;
  static const String add_S_serine;
  static const String add_P_proline;
  static const String add_V_valine;
  static const String add_T_threonine;
  static const String add_C_cysteine;
  static const String add_L_leucine;
  static const String add_I_isoleucine;
  static const String add_N_asparagine;
  static const String add_D_aspartic_acid;
  static const String add_Q_glutamine;
  static const String add_K_lysine;
  static const String add_E_glutamic_acid;
  static const String add_M_methionine;
  static const String add_H_histidine;
  static const String add_F_phenylalanine;
  static const String add_R_arginine;
  static const String add_Y_tyrosine;
  static const String add_W_tryptophan;
  static const String fixed_modifications_unimod;

  // Log level for verbose output
  static const int LOG_LEVEL_VERBOSE;


  TOPPMSFraggerAdapter() :
    SearchEngineBase("MSFraggerAdapter",  "Peptide Identification with MSFragger.\n"
                                  "Important note:\n"
                                  "The Regents of the University of Michigan (\"Michigan\") grants us permission to redistribute    \n"
                                  "the MS Fragger application developed by Michigan within the OpenMS Pipeline and make available \n"
                                  "for use on related service offerings supported by the University of Tubingen and the Center for\n"
                                  "Integrative Bioinformatics.                                                                    \n"
                                  "Per the license agreement the use of the pipeline and associated materials is for academic     \n"
                                  "research, non-commercial or educational purposes. Any commercial use inquiries                 \n"
                                  "must be directed to the University of Michigan Technology Transfer Office at                   \n"
                                  "techtransfer@umich.edu. All right title and interest in MS Fragger shall remain with the       \n"
                                  "University of Michigan.\n"
                                  "\n"
                                  "For details, please see the supplied license file or                                           \n"
				  "https://raw.githubusercontent.com/OpenMS/THIRDPARTY/master/All/MSFragger/License.txt           \n"    
    , false,
             {
                 {"Kong AT, Leprevost FV, Avtonomov DM, Mellacheruvu D, Nesvizhskii AI",
                  "MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry–based proteomics",
                  "Nature Methods volume 14, pages 513–520 (2017)",
                  "doi:10.1038/nmeth.4256"}
             }),
    java_exe(""),
    exe(""),
    parameter_file_path(""),
    input_file(),
    output_file()
  {
  }


protected:

  void registerOptionsAndFlags_() override
  {
    const StringList emptyStrings;
    const std::vector< double > emptyDoubles;

    const StringList validUnits = ListUtils::create<String>("Da,ppm");
    const StringList zero_to_five = ListUtils::create<String>("0,1,2,3,4,5");

    // License agreement
    registerStringOption_(TOPPMSFraggerAdapter::license, "<license>", "", "Set to yes, if you have read and agreed to the MSFragger license terms.", true, false);
    setValidStrings_(TOPPMSFraggerAdapter::license, {"yes","no"});

    // Java executable
    registerInputFile_(TOPPMSFraggerAdapter::java_executable, "<file>", "java", "The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java", false, false, {"skipexists"});
    registerIntOption_(TOPPMSFraggerAdapter::java_heapmemory, "<num>", 3500, "Maximum Java heap size (in MB)", false);

    // Handle executable
    registerInputFile_(TOPPMSFraggerAdapter::executable, "<path_to_executable>", "MSFragger.jar", "Path to the MSFragger executable to use; may be empty if the executable is globally available.", true, false, {"is_executable"});

    // Input file
    registerInputFile_(TOPPMSFraggerAdapter::in, "<file>", "", "Input File with specta for MSFragger");
    setValidFormats_(TOPPMSFraggerAdapter::in, ListUtils::create<String>("mzML,mzXML"));

    // Output file
    registerOutputFile_(TOPPMSFraggerAdapter::out, "<file>", "", "MSFragger output file");
    setValidFormats_(TOPPMSFraggerAdapter::out, ListUtils::create<String>("idXML"), true);

    // Optional output file
    registerOutputFile_(TOPPMSFraggerAdapter::opt_out, "<file>", "", "MSFragger optional output file", false);
    setValidFormats_(TOPPMSFraggerAdapter::opt_out, ListUtils::create<String>("pepXML"), true);
    
    // Path to database to search
    registerInputFile_(TOPPMSFraggerAdapter::database, "<path_to_fasta>", "", "Protein FASTA database file path", true, false);
    setValidFormats_(TOPPMSFraggerAdapter::database, ListUtils::create<String>("FASTA,fasta,fa,fas"), false);

    // TOPP tolerance
    registerTOPPSubsection_("tolerance", "Search Tolerances");

    // Precursor mass tolerance and unit
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::precursor_mass_tolerance_lower, "<precursor_mass_tolerance>", 20.0, "Lower precursor mass tolerance", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::precursor_mass_tolerance_upper, "<precursor_mass_tolerance>", 20.0, "Upper precursor mass tolerance", false, false);
    registerStringOption_(TOPPMSFraggerAdapter::precursor_mass_unit, "<precursor_mass_unit>", "ppm", "Unit of precursor mass tolerance", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::precursor_mass_unit, validUnits);

    // Precursor true tolerance
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::precursor_true_tolerance, "<precursor_true_tolerance>", 0.0, "True precursor mass tolerance (window is +/- this value). Used for tie breaker of results (in spectrally ambiguous cases) and zero bin boosting in open searches (0 disables these features). This option is STRONGLY recommended for open searches.", false, false);
    registerStringOption_(TOPPMSFraggerAdapter::precursor_true_unit, "<precursor_true_unit>", "ppm", "Unit of precursor true tolerance", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::precursor_true_unit, validUnits);

    // Fragment mass tolerance
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::fragment_mass_tolerance, "<fragment_mass_tolerance>", 20.0, "Fragment mass tolerance (window is +/- this value)", false, false);
    registerStringOption_(TOPPMSFraggerAdapter::fragment_mass_unit, "<fragment_mass_unit>", "ppm", "Unit of fragment mass tolerance", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::fragment_mass_unit, validUnits);

    // Isotope error
    registerStringOption_(TOPPMSFraggerAdapter::isotope_error, "<isotope_error>", "0", "Isotope correction for MS/MS events triggered on isotopic peaks. Should be set to 0 (disabled) for open search or 0/1/2 for correction of narrow window searches. Shifts the precursor mass window to multiples of this value multiplied by the mass of C13-C12.", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::isotope_error, ListUtils::create<String>("0,1,2,0/1/2"));

    // TOPP digest
    registerTOPPSubsection_("digest", "In-Silico Digestion Parameters");

    // Enzyme
    StringList enzyme_names;
    ProteaseDB::getInstance()->getAllNames(enzyme_names);
    registerStringOption_(TOPPMSFraggerAdapter::search_enzyme_name, "<search_enzyme_name>", "Trypsin", "Name of the enzyme to be written to the pepXML file", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::search_enzyme_name, enzyme_names);

    // Cut after
    registerStringOption_(TOPPMSFraggerAdapter::search_enzyme_cutafter, "<search_enzyme_cutafter>", "KR", "Residues after which the enzyme cuts (specified as a string of amino acids)", false , false);

    // No cut before
    registerStringOption_(TOPPMSFraggerAdapter::search_enzyme_nocutbefore, "<search_enzyme_nocutbefore>", "P", "Residues that the enzyme will not cut before", false, false);

    // Number of enzyme termini
    registerStringOption_(TOPPMSFraggerAdapter::num_enzyme_termini, "<num_enzyme_termini>", "fully", "Number of enzyme termini (non-enzymatic (0), semi (1), fully (2)", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::num_enzyme_termini, ListUtils::create<String>("non-enzymatic,semi,fully"));

    // Allowed missed cleavages
    registerStringOption_(TOPPMSFraggerAdapter::allowed_missed_cleavage, "<allowed_missed_cleavage>", "2", "Allowed number of missed cleavages", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::allowed_missed_cleavage, zero_to_five); // 5 is the max. allowed value according to MSFragger

    // Digest min length
    _registerNonNegativeInt(TOPPMSFraggerAdapter::digest_min_length, "<digest_min_length>", 7, "Minimum length of peptides to be generated during in-silico digestion", false, false);

    // Digest max length
    _registerNonNegativeInt(TOPPMSFraggerAdapter::digest_max_length, "<digest_max_length>", 64, "Maximum length of peptides to be generated during in-silico digestion", false, false);

    // Digest min mass range
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::digest_mass_range_min, "<digest_mass_range_min>", 500.0, "Min mass of peptides to be generated (Da)", false, false);

    // Digest max mass range
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::digest_mass_range_max, "<digest_mass_range_max>", 5000.0, "Max mass of peptides to be generated (Da)", false, false);

    // TOPP varmod
    registerTOPPSubsection_("varmod", "Variable Modification Parameters");

    // Clip nterm M
    registerFlag_(TOPPMSFraggerAdapter::clip_nterm_m, "Specifies the trimming of a protein N-terminal methionine as a variable modification", false);

    // Modifications
    registerDoubleList_(TOPPMSFraggerAdapter::varmod_masses, "<varmod1_mass .. varmod7_mass>", emptyDoubles , "Masses for variable modifications", false, false);
    registerStringList_(TOPPMSFraggerAdapter::varmod_syntax, "<varmod1_syntax .. varmod7_syntax>", emptyStrings, "Syntax Strings for variable modifications", false, false);
    registerStringList_(TOPPMSFraggerAdapter::variable_modifications_unimod, "<varmod1_unimod .. varmod7_unimod>", emptyStrings, "Variable modifications in unimod syntax, is added to mass+syntax varmod list", false, false);
    registerFlag_(TOPPMSFraggerAdapter::varmod_enable_common, "Enable common variable modifications (15.9949 M and 42.0106 [^)", false);

    // allow_multiple_variable_mods_on_residue
    registerFlag_(TOPPMSFraggerAdapter::not_allow_multiple_variable_mods_on_residue, "Do not allow any one amino acid to be modified by multiple variable modifications", false);

    // Max variable mods per mod
    registerStringOption_(TOPPMSFraggerAdapter::max_variable_mods_per_peptide, "<max_variable_mods_per_peptide>", "2", "Maximum total number of variable modifications per peptide", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::max_variable_mods_per_peptide, zero_to_five);

    // Max variable mods combinations
    _registerNonNegativeInt(TOPPMSFraggerAdapter::max_variable_mods_combinations, "<max_variable_mods_combinations>", 5000, "Maximum allowed number of modified variably modified peptides from each peptide sequence, (maximum of 65534). If a greater number than the maximum is generated, only the unmodified peptide is considered", false, false);
    setMaxInt_(TOPPMSFraggerAdapter::max_variable_mods_combinations, 65534);

    // TOPP spectrum
    registerTOPPSubsection_("spectrum", "Spectrum Processing Parameters");

    _registerNonNegativeInt(TOPPMSFraggerAdapter::minimum_peaks, "<minimum_peaks>", 10, "Minimum number of peaks in experimental spectrum for matching", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::use_topn_peaks, "<use_topN_peaks>", 50, "Pre-process experimental spectrum to only use top N peaks", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::minimum_ratio, "<minimum_ratio>", 0.0, "Filters out all peaks in experimental spectrum less intense than this multiple of the base peak intensity", false, false);
    setMaxFloat_(TOPPMSFraggerAdapter::minimum_ratio, 1.0);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::clear_mz_range_min, "<clear_mz_range_min>", 0.0, "Removes peaks in this m/z range prior to matching (minimum value). Useful for iTRAQ/TMT experiments (i.e. 0.0 150.0)", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::clear_mz_range_max, "<clear_mz_range_max>", 0.0, "Removes peaks in this m/z range prior to matching (maximum value). Useful for iTRAQ/TMT experiments (i.e. 0.0 150.0)", false, false);

    registerStringOption_(TOPPMSFraggerAdapter::max_fragment_charge, "<max_fragment_charge>", "2", "Maximum charge state for theoretical fragments to match", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::max_fragment_charge, ListUtils::create<String>("1,2,3,4"));

    registerFlag_(TOPPMSFraggerAdapter::override_charge, "Ignores precursor charge and uses charge state specified in precursor_charge range (parameters: spectrum:precursor_charge_min and spectrum:precursor_charge_max)" , false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::precursor_charge_min, "<precursor_charge_min>", 1, "Min charge of precursor charge range to consider. If specified, also spectrum:override_charge must be set)" , false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::precursor_charge_max, "<precursor_charge_max>", 4, "Max charge of precursor charge range to consider. If specified, also spectrum:override_charge must be set)" , false, false);

    registerTOPPSubsection_("search", "Open Search Features");

    _registerNonNegativeInt(TOPPMSFraggerAdapter::track_zero_topn, "<track_zero_topn>", 0, "Track top N unmodified peptide results separately from main results internally for boosting features. Should be set to a number greater than search:output_report_topN if zero bin boosting is desired", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::zero_bin_accept_expect, "<zero_bin_accept_expect>", 0.0, "Ranks a zero-bin hit above all non-zero-bin hit if it has expectation less than this value", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::zero_bin_mult_expect, "<zero_bin_mult_expect>", 1.0, "Multiplies expect value of PSMs in the zero-bin during results ordering (set to less than 1 for boosting)", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::add_topn_complementary, "<add_topn_complementary>", 0, "Inserts complementary ions corresponding to the top N most intense fragments in each experimental spectrum. Useful for recovery of modified peptides near C-terminus in open search. 0 disables this option", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::min_fragments_modeling, "<min_fragments_modeling>", 3, "Minimum number of matched peaks in PSM for inclusion in statistical modeling", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::min_matched_fragments, "<min_matched_fragments>", 4, "Minimum number of matched peaks for PSM to be reported. MSFragger recommends a minimum of 4 for narrow window searching and 6 for open searches", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::output_report_topn, "<output_report_topn>", 1, "Reports top N PSMs per input spectrum", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::output_max_expect, "<output_max_expect>", 50.0, "Suppresses reporting of PSM if top hit has expectation greater than this threshold", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::localize_delta_mass, "<localize_delta_mass>", 0, "Include fragment ions mass-shifted by unknown modifications (recommended for open and mass offset searches) (0 for OFF, 1 for ON)", false, false);

    registerTOPPSubsection_("statmod", "Static Modification Parameters");

    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_cterm_peptide, "<add_cterm_peptide>", 0.0, "Statically add mass in Da to C-terminal of peptide", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_nterm_peptide, "<add_nterm_peptide>", 0.0, "Statically add mass in Da to N-terminal of peptide", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_cterm_protein, "<add_cterm_protein>", 0.0, "Statically add mass in Da to C-terminal of protein", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_nterm_protein, "<add_nterm_protein>", 0.0, "Statically add mass in Da to N-terminal of protein", false, false);

    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_G_glycine,       "<add_G_glycine>",       0.0, "Statically add mass to glycine",       false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_A_alanine,       "<add_A_alanine>",       0.0, "Statically add mass to alanine",       false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_S_serine,        "<add_S_serine>",        0.0, "Statically add mass to serine",        false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_P_proline,       "<add_P_proline>",       0.0, "Statically add mass to proline",       false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_V_valine,        "<add_V_valine>",        0.0, "Statically add mass to valine",        false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_T_threonine,     "<add_T_threonine>",     0.0, "Statically add mass to threonine",     false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_C_cysteine,      "<add_C_cysteine>",      57.021464, "Statically add mass to cysteine",      false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_L_leucine,       "<add_L_leucine>",       0.0, "Statically add mass to leucine",       false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_I_isoleucine,    "<add_I_isoleucine>",    0.0, "Statically add mass to isoleucine",    false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_N_asparagine,    "<add_N_asparagine>",    0.0, "Statically add mass to asparagine",    false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_D_aspartic_acid, "<add_D_aspartic_acid>", 0.0, "Statically add mass to aspartic_acid", false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_Q_glutamine,     "<add_Q_glutamine>",     0.0, "Statically add mass to glutamine",     false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_K_lysine,        "<add_K_lysine>",        0.0, "Statically add mass to lysine",        false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_E_glutamic_acid, "<add_E_glutamic_acid>", 0.0, "Statically add mass to glutamic_acid", false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_M_methionine,    "<add_M_methionine>",    0.0, "Statically add mass to methionine",    false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_H_histidine,     "<add_H_histidine>",     0.0, "Statically add mass to histidine",     false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_F_phenylalanine, "<add_F_phenylalanine>", 0.0, "Statically add mass to phenylalanine", false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_R_arginine,      "<add_R_arginine>",      0.0, "Statically add mass to arginine",      false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_Y_tyrosine,      "<add_Y_tyrosine>",      0.0, "Statically add mass to tyrosine",      false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::add_W_tryptophan,    "<add_W_tryptophan>",    0.0, "Statically add mass to tryptophan",    false, true);
    registerStringList_(TOPPMSFraggerAdapter::fixed_modifications_unimod, "<fixedmod1_unimod .. fixedmod7_unimod>", emptyStrings, "Fixed modifications in unimod syntax if specific mass is unknown, e.g. Carbamidomethylation (C). When multiple different masses are given for one aminoacid this parameter (unimod) will have priority.", false, false);

    // register peptide indexing parameter (with defaults for this search engine) TODO: check if search engine defaults are needed
    registerPeptideIndexingParameter_(PeptideIndexing().getParameters());  
  }

  ExitCodes main_(int, const char**) override
  {
    if (getStringOption_(TOPPMSFraggerAdapter::license) != "yes" && !getFlag_("test"))
    {
      _fatalError("MSFragger may only be used upon acceptance of license terms.");
    }

    File::TempDir working_directory(debug_level_ >= 2);
    try
    {
      // java executable
      this->java_exe = this->getStringOption_(TOPPMSFraggerAdapter::java_executable);

      if (!JavaInfo::canRun(this->java_exe, true))
      {
        return ExitCodes::EXTERNAL_PROGRAM_NOTFOUND;
      }

      // executable
      this->exe = this->getStringOption_(TOPPMSFraggerAdapter::executable);

      if (this->exe.empty())
      {
        // looks for MSFRAGGER_PATH in the environment
        QString qmsfragger_path = getenv("MSFRAGGER_PATH");
        if (qmsfragger_path.isEmpty())
        {
          std::cerr << "No executable for MSFragger could be found (also not in MSFRAGGER_PATH)!";
          return ExitCodes::EXTERNAL_PROGRAM_NOTFOUND;
        }
        this->exe = qmsfragger_path;
      }

      // input, output, database name
      const String database = File::absolutePath(this->getStringOption_(TOPPMSFraggerAdapter::database)); // the working dir will be a TMP-dir, so we need absolute paths
      input_file = (this->getStringOption_(TOPPMSFraggerAdapter::in)).toQString();
      output_file = this->getStringOption_(TOPPMSFraggerAdapter::out);
      optional_output_file = this->getStringOption_(TOPPMSFraggerAdapter::opt_out);

      // tolerance
      const double arg_precursor_mass_tolerance_lower(this->getDoubleOption_(TOPPMSFraggerAdapter::precursor_mass_tolerance_lower));
      const double arg_precursor_mass_tolerance_upper(this->getDoubleOption_(TOPPMSFraggerAdapter::precursor_mass_tolerance_upper));
      const String & arg_precursor_mass_unit = this->getStringOption_(TOPPMSFraggerAdapter::precursor_mass_unit);
      const double arg_precursor_true_tolerance(this->getDoubleOption_(TOPPMSFraggerAdapter::precursor_true_tolerance));
      const String & arg_precursor_true_unit = this->getStringOption_(TOPPMSFraggerAdapter::precursor_true_unit);
      const double arg_fragment_mass_tolerance(this->getDoubleOption_(TOPPMSFraggerAdapter::fragment_mass_tolerance));
      const String & arg_fragment_mass_unit = this->getStringOption_(TOPPMSFraggerAdapter::fragment_mass_unit);
      const String & arg_isotope_error = this->getStringOption_(TOPPMSFraggerAdapter::isotope_error);

      // digest
      const String & arg_search_enzyme_name = this->getStringOption_(TOPPMSFraggerAdapter::search_enzyme_name);
      const String & arg_search_enzyme_cutafter = this->getStringOption_(TOPPMSFraggerAdapter::search_enzyme_cutafter);
      const String & arg_search_enzyme_nocutbefore = this->getStringOption_(TOPPMSFraggerAdapter::search_enzyme_nocutbefore);

      std::map< String,int > num_enzyme_termini;
      num_enzyme_termini["non-enzymatic"] = 0;
      num_enzyme_termini["semi"] = 1;
      num_enzyme_termini["fully"] = 2;
      const int arg_num_enzyme_termini = num_enzyme_termini[this->getStringOption_(TOPPMSFraggerAdapter::num_enzyme_termini)];

      const String & arg_allowed_missed_cleavage = this->getStringOption_(TOPPMSFraggerAdapter::allowed_missed_cleavage);
      const int arg_digest_min_length = this->getIntOption_(TOPPMSFraggerAdapter::digest_min_length);
      const int arg_digest_max_length = this->getIntOption_(TOPPMSFraggerAdapter::digest_max_length);
      ensureRange(arg_digest_min_length, arg_digest_max_length, "Maximum length of digest is not allowed to be smaller than minimum length of digest");

      const double arg_digest_mass_range_min = this->getDoubleOption_(TOPPMSFraggerAdapter::digest_mass_range_min);
      const double arg_digest_mass_range_max = this->getDoubleOption_(TOPPMSFraggerAdapter::digest_mass_range_max);
      ensureRange(arg_digest_mass_range_min, arg_digest_mass_range_max, "Maximum digest mass is not allowed to be smaller than minimum digest mass!");

      // varmod
      const bool arg_clip_nterm_m = this->getFlag_(clip_nterm_m);
      std::vector< double > arg_varmod_masses = this->getDoubleList_(TOPPMSFraggerAdapter::varmod_masses);
      std::vector< String > arg_varmod_syntax = this->getStringList_(TOPPMSFraggerAdapter::varmod_syntax);
      std::vector< String > arg_varmod_unimod = this->getStringList_(TOPPMSFraggerAdapter::variable_modifications_unimod);

      // assignment of mass to syntax is by index, so the vectors have to be the same length
      if (arg_varmod_masses.size() != arg_varmod_syntax.size())
      {
        _fatalError("List of arguments for the parameters 'varmod_masses' and 'varmod_syntax' must have the same length!");
      }
      // only up to 7 variable modifications are allowed
      if (arg_varmod_masses.size() > 7)
      {
        _fatalError("MSFragger is restricted to at most 7 variable modifications.");
      }

      // add common variable modifications if requested
      if (this->getFlag_(varmod_enable_common))
      {
        // oxidation on methionine
        this->_addVarMod(arg_varmod_masses, arg_varmod_syntax, 15.9949, "M");

        // N-terminal acetylation
        this->_addVarMod(arg_varmod_masses, arg_varmod_syntax, 42.0106, "[^");
      }

      const bool arg_not_allow_multiple_variable_mods_on_residue = this->getFlag_(TOPPMSFraggerAdapter::not_allow_multiple_variable_mods_on_residue);
      const String & arg_max_variable_mods_per_peptide  = this->getStringOption_(TOPPMSFraggerAdapter::max_variable_mods_per_peptide);
      const int arg_max_variable_mods_combinations = this->getIntOption_(TOPPMSFraggerAdapter::max_variable_mods_combinations);

      // spectrum
      const int arg_minimum_peaks = this->getIntOption_(TOPPMSFraggerAdapter::minimum_peaks);
      const int arg_use_topn_peaks  = this->getIntOption_(TOPPMSFraggerAdapter::use_topn_peaks);
      const double arg_minimum_ratio = this->getDoubleOption_(TOPPMSFraggerAdapter::minimum_ratio);
      const double arg_clear_mz_range_min = this->getDoubleOption_(TOPPMSFraggerAdapter::clear_mz_range_min);
      const double arg_clear_mz_range_max = this->getDoubleOption_(TOPPMSFraggerAdapter::clear_mz_range_max);
      ensureRange(arg_clear_mz_range_min, arg_clear_mz_range_max, "Maximum clear mz value is not allowed to be smaller than minimum clear mz value!");
      const String & arg_max_fragment_charge = this->getStringOption_(TOPPMSFraggerAdapter::max_fragment_charge);
      const bool arg_override_charge = this->getFlag_(TOPPMSFraggerAdapter::override_charge);
      const int arg_precursor_charge_min = this->getIntOption_(TOPPMSFraggerAdapter::precursor_charge_min);
      const int arg_precursor_charge_max = this->getIntOption_(TOPPMSFraggerAdapter::precursor_charge_max);
      ensureRange(arg_precursor_charge_min, arg_precursor_charge_max, "Maximum precursor charge is not allowed to be smaller than minimum precursor charge!");

      // ensures that the user is aware of overriding the precursoe charges
      if ((arg_precursor_charge_min != 1 || arg_precursor_charge_max != 4) && !arg_override_charge)
      {
        _fatalError("If you want to ignore the precursor charge, please also set the -" + override_charge + " flag!");
      }

      // search
      const int arg_track_zero_topn = this->getIntOption_(TOPPMSFraggerAdapter::track_zero_topn);
      const double arg_zero_bin_accept_expect = this->getDoubleOption_(TOPPMSFraggerAdapter::zero_bin_accept_expect);
      const double arg_zero_bin_mult_expect = this->getDoubleOption_(TOPPMSFraggerAdapter::zero_bin_mult_expect);
      const int arg_add_topn_complementary = this->getIntOption_(TOPPMSFraggerAdapter::add_topn_complementary);
      const int arg_min_fragments_modeling = this->getIntOption_(TOPPMSFraggerAdapter::min_fragments_modeling);
      const int arg_min_matched_fragments = this->getIntOption_(TOPPMSFraggerAdapter::min_matched_fragments);
      const int arg_output_report_topn = this->getIntOption_(TOPPMSFraggerAdapter::output_report_topn);
      const double arg_output_max_expect = this->getDoubleOption_(TOPPMSFraggerAdapter::output_max_expect);
      const int arg_localize_delta_mass = this->getIntOption_(TOPPMSFraggerAdapter::localize_delta_mass);
      
      // statmod
      double arg_add_cterm_peptide = this->getDoubleOption_(TOPPMSFraggerAdapter::add_cterm_peptide);
      double arg_add_nterm_peptide = this->getDoubleOption_(TOPPMSFraggerAdapter::add_nterm_peptide);
      double arg_add_cterm_protein = this->getDoubleOption_(TOPPMSFraggerAdapter::add_cterm_protein);
      double arg_add_nterm_protein = this->getDoubleOption_(TOPPMSFraggerAdapter::add_nterm_protein);
      double arg_add_G_glycine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_G_glycine);
      double arg_add_A_alanine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_A_alanine);
      double arg_add_S_serine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_S_serine);
      double arg_add_P_proline = this->getDoubleOption_(TOPPMSFraggerAdapter::add_P_proline);
      double arg_add_V_valine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_V_valine);
      double arg_add_T_threonine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_T_threonine);
      double arg_add_C_cysteine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_C_cysteine);
      double arg_add_L_leucine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_L_leucine);
      double arg_add_I_isoleucine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_I_isoleucine);
      double arg_add_N_asparagine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_N_asparagine);
      double arg_add_D_aspartic_acid = this->getDoubleOption_(TOPPMSFraggerAdapter::add_D_aspartic_acid);
      double arg_add_Q_glutamine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_Q_glutamine);
      double arg_add_K_lysine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_K_lysine);
      double arg_add_E_glutamic_acid = this->getDoubleOption_(TOPPMSFraggerAdapter::add_E_glutamic_acid);
      double arg_add_M_methionine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_M_methionine);
      double arg_add_H_histidine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_H_histidine);
      double arg_add_F_phenylalanine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_F_phenylalanine);
      double arg_add_R_arginine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_R_arginine);
      double arg_add_Y_tyrosine = this->getDoubleOption_(TOPPMSFraggerAdapter::add_Y_tyrosine);
      double arg_add_W_tryptophan = this->getDoubleOption_(TOPPMSFraggerAdapter::add_W_tryptophan);
      std::vector< String > arg_fixmod_unimod = this->getStringList_(TOPPMSFraggerAdapter::fixed_modifications_unimod);

      // parameters have been read in and verified, they are now going to be written into the fragger.params file in a temporary directory
      this->parameter_file_path = File::getTemporaryFile();

      writeDebug_("Parameter file for MSFragger: '" + this->parameter_file_path + "'", TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE);
      writeDebug_("Working Directory: '" + working_directory.getPath() + "'", TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE);
      writeDebug_("If you want to keep the working directory and the parameter file, set the -debug to 2", 1);
      ofstream os(this->parameter_file_path.c_str());


      // Write all the parameters into the file
      os << "database_name = " << String(database)
                               << "\nnum_threads = " << this->getIntOption_("threads")
                               << "\n\nprecursor_mass_lower = " << (-arg_precursor_mass_tolerance_lower)
                               << "\nprecursor_mass_upper = " << arg_precursor_mass_tolerance_upper
                               << "\nprecursor_mass_units = " << (arg_precursor_mass_unit == "Da" ? 0 : 1)
                               << "\nprecursor_true_tolerance = " << arg_precursor_true_tolerance
                               << "\nprecursor_true_units = " << (arg_precursor_true_unit == "Da" ? 0 : 1)
                               << "\nfragment_mass_tolerance = " << arg_fragment_mass_tolerance
                               << "\nfragment_mass_units = " << (arg_fragment_mass_unit == "Da" ? 0 : 1)
                               << "\n\nisotope_error = " << arg_isotope_error
                               << "\n\nsearch_enzyme_name = " << arg_search_enzyme_name
                               << "\nsearch_enzyme_cutafter = " << arg_search_enzyme_cutafter
                               << "\nsearch_enzyme_butnotafter = " << arg_search_enzyme_nocutbefore
                               << "\n\nnum_enzyme_termini = " << arg_num_enzyme_termini
                               << "\nallowed_missed_cleavage = " << arg_allowed_missed_cleavage
                               << "\n\nclip_nTerm_M = " << arg_clip_nterm_m << '\n';

      // Write variable modifications from masses/syntax and unimod to unique set (and also write to log)
      writeLogInfo_("Variable Modifications set to:");
      std::set< std::pair< double, String > > varmods_combined;
      Size i;
      for (i = 0; i < arg_varmod_masses.size(); ++i)
      { 
        std::pair <double, String> tmp_mod = std::make_pair (arg_varmod_masses[i], arg_varmod_syntax[i]);
        varmods_combined.insert(tmp_mod);
      }


      // TODO Move to Modified Peptide Generator
      if (!arg_varmod_unimod.empty())
      {
        // String filter for terminal aminoacid modification, delete mod from String list, continue with other unimods
        std::vector< String > n_terminal_aa_mods;
        std::vector< String > c_terminal_aa_mods;
        std::vector< int > n_terminal_aa_mods_toDel;
        std::vector< int > c_terminal_aa_mods_toDel;
        for (Size i=0; i<arg_varmod_unimod.size(); i++)
        {
          int nt = arg_varmod_unimod[i].find(" (N-term");
          int ct = arg_varmod_unimod[i].find(" (C-term");

          if (!(nt == -1 && ct == -1)) // has -term modification
          {
            int closed_arg = arg_varmod_unimod[i].find("term)"); // Check if the terminal argument is closed or continued with aminoacid
            if (closed_arg == -1)
            {
              int j = arg_varmod_unimod[i].find("-term");
              if (arg_varmod_unimod[i].substr(j+7)!=")")
              {
                _fatalError("Multiple aminoacids in terminal modification are not allowed");
              }
              String res = arg_varmod_unimod[i].substr(j+6, 1);
              String mod = arg_varmod_unimod[i].substr(0, j-3);
              String modificationString = mod.append(" (").append(res).append(")");
              if (nt != -1)
              {
                n_terminal_aa_mods.push_back(modificationString);
                n_terminal_aa_mods_toDel.push_back(i);
              }
              if (ct != -1)
              {
                c_terminal_aa_mods.push_back(modificationString);
                c_terminal_aa_mods_toDel.push_back(i);
              }
            }
          }         
        }

        // Write the variable modification in correct syntax to a combined list and delete element from parameter list
        const ModifiedPeptideGenerator::MapToResidueType n_var_mod_temp = ModifiedPeptideGenerator::getModifications(n_terminal_aa_mods);
        for (auto const & r : n_var_mod_temp.val)
        {
          const double deltamass = r.first->getDiffMonoMass();
          const String res = r.second->getOneLetterCode();
          std::pair <double, String> tmp_mod = std::make_pair (deltamass, "n" + res);
          varmods_combined.insert(tmp_mod);
        }
        
        for (auto const & i : n_terminal_aa_mods_toDel)
        {
          arg_varmod_unimod.erase(arg_varmod_unimod.begin()+i);
        }

        const ModifiedPeptideGenerator::MapToResidueType c_var_mod_temp = ModifiedPeptideGenerator::getModifications(c_terminal_aa_mods);
        for (auto const & r : c_var_mod_temp.val)
        {
          const double deltamass = r.first->getDiffMonoMass();
          const String res = r.second->getOneLetterCode();
          std::pair <double, String> tmp_mod = std::make_pair (deltamass, "c" + res);
          varmods_combined.insert(tmp_mod);
        }
        
        for (auto const & i : c_terminal_aa_mods_toDel)
        {
          arg_varmod_unimod.erase(arg_varmod_unimod.begin()+i);
        }

        // Collect all other modifications and filter true terminal modifications for correct syntax in MSFragger
        const ModifiedPeptideGenerator::MapToResidueType variable_mod = ModifiedPeptideGenerator::getModifications(arg_varmod_unimod);
        for (auto const & r : variable_mod.val)
        { 
          String res;
          const double deltamass = r.first->getDiffMonoMass();
          if (r.first->getTermSpecificity() == ResidueModification::N_TERM)
          {
            res = "n^";
          }
          else if (r.first->getTermSpecificity() == ResidueModification::C_TERM)
          {
            res = "c^";
          }
          else if (r.first->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
          {
            res = "[^";
          }
          else if (r.first->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
          {
            res = "]^";
          }
          else
          {
            res = r.second->getOneLetterCode();
          }
          std::pair <double, String> tmp_mod = std::make_pair (deltamass, res);
          varmods_combined.insert(tmp_mod);
        }
      }
      i = 0;
      for (auto const & m : varmods_combined)
      {
        const String varmod = "variable_mod_0" + String(i+1) + " = " + String(m.first) + " " + String(m.second);
        os << "\n" << varmod;
        writeLogInfo_(varmod);
        i++;
      }

      // collect all unimod fixed modifications and specify deltamass for each aminoacid
      if (!arg_fixmod_unimod.empty())
      {
        const ModifiedPeptideGenerator::MapToResidueType fixed_mod = ModifiedPeptideGenerator::getModifications(arg_fixmod_unimod);
        for (auto const & r : fixed_mod.val)
        {
          const double deltamass = r.first->getDiffMonoMass();
          if (r.first->getTermSpecificity() == ResidueModification::N_TERM)
          {
            arg_add_nterm_peptide = deltamass;
          }
          else if (r.first->getTermSpecificity() == ResidueModification::C_TERM)
          {
            arg_add_cterm_peptide = deltamass;
          }
          else if (r.first->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
          {
            arg_add_nterm_protein = deltamass;
          }
          else if (r.first->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
          {
            arg_add_cterm_protein = deltamass;
          } 
          else
          {
            const String res = r.second->getOneLetterCode();
            switch(res[0])
            {
              case 'G':
                arg_add_G_glycine = deltamass;
                break;
              case 'A':
                arg_add_A_alanine = deltamass;
                break;
              case 'S':
                arg_add_S_serine = deltamass;
                break;
              case 'P':
                arg_add_P_proline = deltamass;
                break;
              case 'V':
                arg_add_V_valine = deltamass;
                break;
              case 'T':
                arg_add_T_threonine = deltamass;
                break;
              case 'C':
                arg_add_C_cysteine = deltamass;
                break;
              case 'L':
                arg_add_L_leucine = deltamass;
                break;
              case 'I':
                arg_add_I_isoleucine = deltamass;
                break;
              case 'N':
                arg_add_N_asparagine = deltamass;
                break;
              case 'D':
                arg_add_D_aspartic_acid = deltamass;
                break;
              case 'Q':
                arg_add_Q_glutamine = deltamass;
                break;
              case 'K':
                arg_add_K_lysine = deltamass;
                break;
              case 'E':
                arg_add_E_glutamic_acid = deltamass;
                break;
              case 'M':
                arg_add_M_methionine = deltamass;
                break;
              case 'H':
                arg_add_H_histidine = deltamass;
                break;
              case 'F':
                arg_add_F_phenylalanine = deltamass;
                break;
              case 'R':
                arg_add_R_arginine = deltamass;
                break;
              case 'Y':
                arg_add_Y_tyrosine = deltamass;
                break;
              case 'W':
                arg_add_W_tryptophan = deltamass;
                break;
            }
          }
        }
      }

      os << std::endl
          << "\nallow_multiple_variable_mods_on_residue = " << (arg_not_allow_multiple_variable_mods_on_residue ? 0 : 1)
          << "\nmax_variable_mods_per_peptide = " << arg_max_variable_mods_per_peptide
          << "\nmax_variable_mods_combinations = " << arg_max_variable_mods_combinations
          << "\n\noutput_file_extension = " << "pepXML"
          << "\noutput_format = " << "pepXML"
          << "\noutput_report_topN = " << arg_output_report_topn
          << "\noutput_max_expect = " << arg_output_max_expect
          << "\n\nprecursor_charge = " << arg_precursor_charge_min << " " << arg_precursor_charge_max
          << "\noverride_charge = " << (arg_override_charge ? 1  : 0)
          << "\n\ndigest_min_length = " << arg_digest_min_length
          << "\ndigest_max_length = " << arg_digest_max_length
          << "\ndigest_mass_range = " << arg_digest_mass_range_min << " " << arg_digest_mass_range_max
          << "\nmax_fragment_charge = " << arg_max_fragment_charge
          << "\n\ntrack_zero_topN = " << arg_track_zero_topn
          << "\nzero_bin_accept_expect = " << arg_zero_bin_accept_expect
          << "\nzero_bin_mult_expect = " << arg_zero_bin_mult_expect
          << "\nadd_topN_complementary = " << arg_add_topn_complementary
          << "\n\nminimum_peaks = " << arg_minimum_peaks
          << "\nuse_topN_peaks = " << arg_use_topn_peaks
          << "\nlocalize_delta_mass = " << arg_localize_delta_mass
          << "\nmin_fragments_modelling = " << arg_min_fragments_modeling
          << "\nmin_matched_fragments = " << arg_min_matched_fragments
          << "\nminimum_ratio = " << arg_minimum_ratio
          << "\nclear_mz_range = " << arg_clear_mz_range_min << " " << arg_clear_mz_range_max
          << "\nadd_Cterm_peptide = " << arg_add_cterm_peptide
          << "\nadd_Nterm_peptide = " << arg_add_nterm_peptide
          << "\nadd_Cterm_protein = " << arg_add_cterm_protein
          << "\nadd_Nterm_protein = " << arg_add_nterm_protein
          << "\n\nadd_G_glycine = " << arg_add_G_glycine
          << "\nadd_A_alanine = " << arg_add_A_alanine
          << "\nadd_S_serine = " << arg_add_S_serine
          << "\nadd_P_proline = " << arg_add_P_proline
          << "\nadd_V_valine = " << arg_add_V_valine
          << "\nadd_T_threonine = " << arg_add_T_threonine
          << "\nadd_C_cysteine = " << arg_add_C_cysteine
          << "\nadd_L_leucine = " << arg_add_L_leucine
          << "\nadd_I_isoleucine = " << arg_add_I_isoleucine
          << "\nadd_N_asparagine = " << arg_add_N_asparagine
          << "\nadd_D_aspartic_acid = " << arg_add_D_aspartic_acid
          << "\nadd_Q_glutamine = " << arg_add_Q_glutamine
          << "\nadd_K_lysine = " << arg_add_K_lysine
          << "\nadd_E_glutamic_acid = " << arg_add_E_glutamic_acid
          << "\nadd_M_methionine = " << arg_add_M_methionine
          << "\nadd_H_histidine = " << arg_add_H_histidine
          << "\nadd_F_phenylalanine = " << arg_add_F_phenylalanine
          << "\nadd_R_arginine = " << arg_add_R_arginine
          << "\nadd_Y_tyrosine = " << arg_add_Y_tyrosine
          << "\nadd_W_tryptophan = " << arg_add_W_tryptophan;
      os.close();
    }
    catch (int)
    {
      return ILLEGAL_PARAMETERS;
    }

    QStringList process_params; // the actual process is Java, not MSFragger
    process_params << "-Xmx" + QString::number(this->getIntOption_(java_heapmemory)) + "m"
        << "-jar" << this->exe.toQString()
        << this->parameter_file_path.toQString()
        << input_file;

    if (this->debug_level_ >= TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE)
    {
      writeDebug_("COMMAND LINE CALL IS:", 1);
      String command_line = this->java_exe;
      for (const auto& process_param : process_params)
      {
        command_line += (" " + process_param);
      }
      writeDebug_(command_line, TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE);
    }

    TOPPBase::ExitCodes exit_code = runExternalProcess_(java_exe.toQString(), process_params, working_directory.getPath().toQString());
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    // convert from pepXML to idXML
    String pepxmlfile = FileHandler::swapExtension(input_file, FileTypes::PEPXML);
    std::vector<PeptideIdentification> peptide_identifications;
    std::vector<ProteinIdentification> protein_identifications;
    PepXMLFile().load(pepxmlfile, protein_identifications, peptide_identifications);
    for (auto it = protein_identifications.begin(); it != protein_identifications.end(); it++)
    { 
        it->setSearchEngine("MSFragger");
        //Whatever the pepXML says, overwrite origin as the input mzML
        it->setPrimaryMSRunPath({this->getStringOption_(TOPPMSFraggerAdapter::in)}, false);
    }

    // write all (!) parameters as metavalues to the search parameters
    if (!protein_identifications.empty())
    {
      DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_identifications[0].getSearchParameters(), this->getToolPrefix());
    }

    // if "reindex" parameter is set to true will perform reindexing
    if (auto ret = reindex_(protein_identifications, peptide_identifications); ret != EXECUTION_OK) return ret;

    // add percolator features
    StringList feature_set;
    PercolatorFeatureSetHelper::addMSFRAGGERFeatures(feature_set);
    protein_identifications.front().getSearchParameters().setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    FileHandler().storeIdentifications(output_file, protein_identifications, peptide_identifications, {FileTypes::IDXML});

    // remove the msfragger pepXML output from the user location
    if (optional_output_file.empty())
    {
      File::remove(pepxmlfile);
    }
    else
    { // rename the pepXML file to the opt_out
      File::rename(pepxmlfile.toQString(), optional_output_file.toQString()); 
    }

    // remove ".pepindex" database file
    if (this->debug_level_ < 2)
    {
      String db_index = this->getStringOption_(TOPPMSFraggerAdapter::database) + ".1.pepindex"; 
      File::remove(db_index);
    }
    return EXECUTION_OK;
  }

private:
  String java_exe;
  String exe;

  String parameter_file_path;
  QString input_file;
  String output_file;
  String optional_output_file;

  // adds variable modification if not already present
  void _addVarMod(std::vector< double > & masses, std::vector< String > & syntaxes, const double mass, const String & syntax) const
  {
    const std::vector< double >::iterator it1 = std::find(masses.begin(), masses.end(), mass);
    const std::vector< String >::iterator it2 = std::find(syntaxes.begin(), syntaxes.end(), syntax);

    // add the provided variable modification if not already present
    if (   it1 == masses.end()
        || it2 == syntaxes.end()
        || std::distance(masses.begin(), it1) != std::distance(syntaxes.begin(), it2))
    {
      masses.push_back(mass);
      syntaxes.push_back(syntax);
    }
  }

  inline void _registerNonNegativeInt(const String & param_name, const String & argument, const int default_value, const String & description, const bool required, const bool advanced)
  {
    this->registerIntOption_(param_name, argument, default_value, description, required, advanced);
    this->setMinInt_(param_name, 0);
  }

  inline void _registerNonNegativeDouble(const String & param_name, const String & argument, const double default_value, const String & description, const bool required, const bool advanced)
  {
    this->registerDoubleOption_(param_name, argument, default_value, description, required, advanced);
    this->setMinFloat_(param_name, 0.0);
  }


  inline void _fatalError(const String & message)
  {
    OPENMS_LOG_FATAL_ERROR << "FATAL: " << message << std::endl;
    throw 1;
  }


  void checkUnique(const StringList & elements, const String & message)
  {
    for (Size i = 0; i < elements.size(); ++i)
    {
      for (Size j = 0; j < i; ++j)
      {
        if (elements[i] == elements[j])
       {
          _fatalError(message);
        }
      }
    }
  }

  inline void ensureRange(const double left, const double right, const String & message) const
  {
    if (right < left)
    {
      OPENMS_LOG_ERROR << "FATAL: " << message << std::endl;
      throw 1;
    }
  }
};

const String TOPPMSFraggerAdapter::java_executable = "java_executable";
const String TOPPMSFraggerAdapter::java_heapmemory = "java_heapmemory";
const String TOPPMSFraggerAdapter::executable = "executable";
const String TOPPMSFraggerAdapter::in = "in";
const String TOPPMSFraggerAdapter::out = "out";
const String TOPPMSFraggerAdapter::opt_out = "opt_out";
const String TOPPMSFraggerAdapter::database = "database";

// tolerance
const String TOPPMSFraggerAdapter::precursor_mass_tolerance_lower = "tolerance:precursor_mass_tolerance_lower";
const String TOPPMSFraggerAdapter::precursor_mass_tolerance_upper = "tolerance:precursor_mass_tolerance_upper";
const String TOPPMSFraggerAdapter::precursor_mass_unit = "tolerance:precursor_mass_unit";
const String TOPPMSFraggerAdapter::precursor_true_tolerance = "tolerance:precursor_true_tolerance";
const String TOPPMSFraggerAdapter::precursor_true_unit = "tolerance:precursor_true_unit";
const String TOPPMSFraggerAdapter::fragment_mass_tolerance = "tolerance:fragment_mass_tolerance";
const String TOPPMSFraggerAdapter::fragment_mass_unit = "tolerance:fragment_mass_unit";
const String TOPPMSFraggerAdapter::isotope_error = "tolerance:isotope_error";

// digest
const String TOPPMSFraggerAdapter::search_enzyme_name = "digest:search_enzyme_name";
const String TOPPMSFraggerAdapter::search_enzyme_cutafter = "digest:search_enzyme_cutafter";
const String TOPPMSFraggerAdapter::search_enzyme_nocutbefore = "digest:search_enzyme_nocutbefore";
const String TOPPMSFraggerAdapter::num_enzyme_termini = "digest:num_enzyme_termini";
const String TOPPMSFraggerAdapter::allowed_missed_cleavage = "digest:allowed_missed_cleavage";
const String TOPPMSFraggerAdapter::digest_min_length = "digest:min_length";
const String TOPPMSFraggerAdapter::digest_max_length = "digest:max_length";
const String TOPPMSFraggerAdapter::digest_mass_range_min = "digest:mass_range_min";
const String TOPPMSFraggerAdapter::digest_mass_range_max = "digest:mass_range_max";

// varmod
const String TOPPMSFraggerAdapter::clip_nterm_m = "varmod:clip_nterm_m";
const String TOPPMSFraggerAdapter::varmod_masses = "varmod:masses";
const String TOPPMSFraggerAdapter::varmod_syntax = "varmod:syntaxes";
const String TOPPMSFraggerAdapter::varmod_enable_common = "varmod:enable_common";
const String TOPPMSFraggerAdapter::not_allow_multiple_variable_mods_on_residue = "varmod:not_allow_multiple_variable_mods_on_residue";
const String TOPPMSFraggerAdapter::max_variable_mods_per_peptide = "varmod:max_variable_mods_per_peptide";
const String TOPPMSFraggerAdapter::max_variable_mods_combinations = "varmod:max_variable_mods_combinations";
const String TOPPMSFraggerAdapter::variable_modifications_unimod = "varmod:unimod";

// spectrum
const String TOPPMSFraggerAdapter::minimum_peaks = "spectrum:minimum_peaks";
const String TOPPMSFraggerAdapter::use_topn_peaks = "spectrum:use_topn_peaks";
const String TOPPMSFraggerAdapter::minimum_ratio = "spectrum:minimum_ratio";
const String TOPPMSFraggerAdapter::clear_mz_range_min = "spectrum:clear_mz_range_min";
const String TOPPMSFraggerAdapter::clear_mz_range_max = "spectrum:clear_mz_range_max";
const String TOPPMSFraggerAdapter::max_fragment_charge = "spectrum:max_fragment_charge";
const String TOPPMSFraggerAdapter::override_charge = "spectrum:override_charge";
const String TOPPMSFraggerAdapter::precursor_charge_min = "spectrum:precursor_charge_min";
const String TOPPMSFraggerAdapter::precursor_charge_max = "spectrum:precursor_charge_max";

// search
const String TOPPMSFraggerAdapter::track_zero_topn = "search:track_zero_topn";
const String TOPPMSFraggerAdapter::zero_bin_accept_expect = "search:zero_bin_accept_expect";
const String TOPPMSFraggerAdapter::zero_bin_mult_expect = "search:zero_bin_mult_expect";
const String TOPPMSFraggerAdapter::add_topn_complementary = "search:add_topn_complementary";
const String TOPPMSFraggerAdapter::min_fragments_modeling = "search:min_fragments_modeling";
const String TOPPMSFraggerAdapter::min_matched_fragments = "search:min_matched_fragments";
const String TOPPMSFraggerAdapter::output_report_topn = "search:output_report_topn";
const String TOPPMSFraggerAdapter::output_max_expect = "search:output_max_expect";
const String TOPPMSFraggerAdapter::localize_delta_mass = "search:localize_delta_mass";

// statmod
const String TOPPMSFraggerAdapter::add_cterm_peptide = "statmod:add_cterm_peptide";
const String TOPPMSFraggerAdapter::add_nterm_peptide = "statmod:add_nterm_peptide";
const String TOPPMSFraggerAdapter::add_cterm_protein = "statmod:add_cterm_protein";
const String TOPPMSFraggerAdapter::add_nterm_protein = "statmod:add_nterm_protein";
const String TOPPMSFraggerAdapter::add_G_glycine = "statmod:add_G_glycine";
const String TOPPMSFraggerAdapter::add_A_alanine = "statmod:add_A_alanine";
const String TOPPMSFraggerAdapter::add_S_serine = "statmod:add_S_serine";
const String TOPPMSFraggerAdapter::add_P_proline = "statmod:add_P_proline";
const String TOPPMSFraggerAdapter::add_V_valine = "statmod:add_V_valine";
const String TOPPMSFraggerAdapter::add_T_threonine = "statmod:add_T_threonine";
const String TOPPMSFraggerAdapter::add_C_cysteine = "statmod:add_C_cysteine";
const String TOPPMSFraggerAdapter::add_L_leucine = "statmod:add_L_leucine";
const String TOPPMSFraggerAdapter::add_I_isoleucine = "statmod:add_I_isoleucine";
const String TOPPMSFraggerAdapter::add_N_asparagine = "statmod:add_N_asparagine";
const String TOPPMSFraggerAdapter::add_D_aspartic_acid = "statmod:add_D_aspartic_acid";
const String TOPPMSFraggerAdapter::add_Q_glutamine = "statmod:add_Q_glutamine";
const String TOPPMSFraggerAdapter::add_K_lysine = "statmod:add_K_lysine";
const String TOPPMSFraggerAdapter::add_E_glutamic_acid = "statmod:add_E_glutamic_acid";
const String TOPPMSFraggerAdapter::add_M_methionine = "statmod:add_M_methionine";
const String TOPPMSFraggerAdapter::add_H_histidine = "statmod:add_H_histidine";
const String TOPPMSFraggerAdapter::add_F_phenylalanine = "statmod:add_F_phenylalanine";
const String TOPPMSFraggerAdapter::add_R_arginine = "statmod:add_R_arginine";
const String TOPPMSFraggerAdapter::add_Y_tyrosine = "statmod:add_Y_tyrosine";
const String TOPPMSFraggerAdapter::add_W_tryptophan = "statmod:add_W_tryptophan";
const String TOPPMSFraggerAdapter::license = "license";
const String TOPPMSFraggerAdapter::fixed_modifications_unimod = "statmod:unimod";


const int TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE = 1;

int main(int argc, const char** argv)
{
  TOPPMSFraggerAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
