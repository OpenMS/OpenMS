// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann, Leon Bichmann $
// --------------------------------------------------------------------------
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/SYSTEM/JavaInfo.h>
#include <QtCore/QDir>
#include <QtCore/QProcess>
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
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MSFraggerAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em MSFragger must be installed before this adapter can be used.

	All MSFragger parameters (as specified in the fragger.params file) have been transcribed to parameters of this OpenMS util.
	It is not possible to provide an explicit fragger.params file to avoid redundancy with the ini file.
	This adapter creates an fragger.params file prior to calling MSFragger. If the fragger.params file should be inspected, set the
	-debug option to 2. MSFraggerAdapter will print the path to the working directory to standard out.

	MSFragger can process multiple input files (mzML, mzXML) one after another. The number of output files specified must match
	the number of input spectra files. The output file is then matched to the input file by index. The default parameters of the
	adapter are the same as given by the official MSFragger manual.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MSFraggerAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MSFraggerAdapter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPMSFraggerAdapter final :
public TOPPBase
{
public:

  static const String param_java_executable;
  static const String param_java_heapmemory;
  static const String param_executable;
  static const String param_in;
  static const String param_out;
  static const String param_database;

  // tolerance
  static const String param_precursor_mass_tolerance;
  static const String param_precursor_mass_unit;
  static const String param_precursor_true_tolerance;
  static const String param_precursor_true_unit;
  static const String param_fragment_mass_tolerance;
  static const String param_fragment_mass_unit;
  static const String param_isotope_error;

  // digest
  static const String param_search_enzyme_name;
  static const String param_search_enzyme_cutafter;
  static const String param_search_enzyme_nocutbefore;
  static const String param_num_enzyme_termini;
  static const String param_allowed_missed_cleavage;
  static const String param_digest_min_length;
  static const String param_digest_max_length;
  static const String param_digest_mass_range_min;
  static const String param_digest_mass_range_max;

  // varmod
  static const String param_clip_nterm_m;
  static const String param_varmod_masses;
  static const String param_varmod_syntax;
  static const String param_varmod_enable_common;
  static const String param_not_allow_multiple_variable_mods_on_residue;
  static const String param_max_variable_mods_per_mod;
  static const String param_max_variable_mods_combinations;

  // spectrum
  static const String param_minimum_peaks;
  static const String param_use_topn_peaks;
  static const String param_minimum_ratio;
  static const String param_clear_mz_range_min;
  static const String param_clear_mz_range_max;
  static const String param_max_fragment_charge;
  static const String param_override_charge;
  static const String param_precursor_charge_min;
  static const String param_precursor_charge_max;

  // search
  static const String param_track_zero_topn;
  static const String param_zero_bin_accept_expect;
  static const String param_zero_bin_mult_expect;
  static const String param_add_topn_complementary;
  static const String param_min_fragments_modeling;
  static const String param_min_matched_fragments;
  static const String param_output_report_topn;
  static const String param_output_max_expect;

  // statmod
  static const String param_add_cterm_peptide;
  static const String param_add_nterm_peptide;
  static const String param_add_cterm_protein;
  static const String param_add_nterm_protein;
  static const String param_add_G_glycine;
  static const String param_add_A_alanine;
  static const String param_add_S_serine;
  static const String param_add_P_proline;
  static const String param_add_V_valine;
  static const String param_add_T_threonine;
  static const String param_add_C_cysteine;
  static const String param_add_L_leucine;
  static const String param_add_I_isoleucine;
  static const String param_add_N_asparagine;
  static const String param_add_D_aspartic_acid;
  static const String param_add_Q_glutamine;
  static const String param_add_K_lysine;
  static const String param_add_E_glutamic_acid;
  static const String param_add_M_methionine;
  static const String param_add_H_histidine;
  static const String param_add_F_phenylalanine;
  static const String param_add_R_arginine;
  static const String param_add_Y_tyrosine;
  static const String param_add_W_tryptophan;

  // Log level for verbose output
  static const int LOG_LEVEL_VERBOSE;


  TOPPMSFraggerAdapter() :
    TOPPBase("MSFraggerAdapter", "Peptide Identification with MSFragger", false),
    working_directory(""),
    java_executable(""),
    executable(""),
    parameter_file_path(""),
    input_files(),
    output_files(),
    file_type("")
  {
  }

  ~TOPPMSFraggerAdapter()
  {
    // Remove the temp working directory if the debug level is smaller than 2
    // and the working directory has been created
    if (this->debug_level_ < 2 && this->working_directory != "")
    {
      File::removeDir(this->working_directory);
    }
  }


protected:
  void registerOptionsAndFlags_() final override
  {
    const StringList emptyStrings;
    const std::vector< double > emptyDoubles;

    const StringList validUnits = ListUtils::create<String>("Da,ppm");
    const StringList isotope_error_and_enzyme_termini = ListUtils::create<String>("0,1,2");
    const StringList zero_to_five = ListUtils::create<String>("0,1,2,3,4,5");

    // Java executable
    registerInputFile_(TOPPMSFraggerAdapter::param_java_executable, "<file>", "java", "The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java", false, false, ListUtils::create<String>("skipexists"));
    registerIntOption_(TOPPMSFraggerAdapter::param_java_heapmemory, "<num>", 3500, "Maximum Java heap size (in MB)", false);

    // Handle executable
    registerInputFile_(TOPPMSFraggerAdapter::param_executable, "<path_to_executable>", "", "Path to the MSFragger executable to use; may be empty if the executable is globally available.", false, false, ListUtils::create<String>("skipexists"));

    // input spectra
    registerInputFileList_(TOPPMSFraggerAdapter::param_in, "spectra_file_1 [spectra_file_2, .., spectra_file_N]>", emptyStrings, "Spectra files to search with MSFragger", true, false);
    setValidFormats_(TOPPMSFraggerAdapter::param_in, ListUtils::create<String>("mzML,mzXML"));

    //  output files
    registerOutputFileList_(TOPPMSFraggerAdapter::param_out, "output_file_1 [output_file_2, ..., output_file_N]", emptyStrings, "MSFragger output files", true, false);
    setValidFormats_(TOPPMSFraggerAdapter::param_out, ListUtils::create<String>("pep.xml,pepXML,tsv"), true);

    // Path to database to search
    registerInputFile_(TOPPMSFraggerAdapter::param_database, "<path_to_fasta>", "", "Protein FASTA database file path", true, false);
    setValidFormats_(TOPPMSFraggerAdapter::param_database, ListUtils::create<String>("FASTA,fasta,fa,fas"), false);

    // TOPP tolerance
    registerTOPPSubsection_("tolerance", "Search Tolerances");

    // Precursor mass tolerance and unit
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_precursor_mass_tolerance, "<precursor_mass_tolerance>", 20.0, "Precursor mass tolerance (window is +/- this value)", false, false);
    registerStringOption_(TOPPMSFraggerAdapter::param_precursor_mass_unit, "<precursor_mass_unit>", "ppm", "Unit of precursor mass tolerance", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_precursor_mass_unit, validUnits);

    // Precursor true tolerance
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_precursor_true_tolerance, "<precursor_true_tolerance>", 0.0, "True precursor mass tolerance (window is +/- this value). Used for tie breaker of results (in spectrally ambiguous cases) and zero bin boosting in open searches (0 disables these features). This option is STRONGLY recommended for open searches.", false, false);
    registerStringOption_(TOPPMSFraggerAdapter::param_precursor_true_unit, "<precursor_true_unit>", "ppm", "Unit of precursor true tolerance", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_precursor_true_unit, validUnits);

    // Fragment mass tolerance
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_fragment_mass_tolerance, "<fragment_mass_tolerance>", 20.0, "Fragment mass tolerance (window is +/- this value)", false, false);
    registerStringOption_(TOPPMSFraggerAdapter::param_fragment_mass_unit, "<fragment_mass_unit>", "ppm", "Unit of fragment mass tolerance", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_fragment_mass_unit, validUnits);

    // Isotope error
    registerStringOption_(TOPPMSFraggerAdapter::param_isotope_error, "<isotope_error>", "0", "Isotope correction for MS/MS events triggered on isotopic peaks. Should be set to 0 (disabled) for open search or 0/1/2 for correction of narrow window searches. Shifts the precursor mass window to multiples of this value multiplied by the mass of C13-C12.", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_isotope_error, isotope_error_and_enzyme_termini);

    // TOPP digest
    registerTOPPSubsection_("digest", "In-Silico Digestion Parameters");

    // Enzyme
    StringList enzyme_names;
    ProteaseDB::getInstance()->getAllNames(enzyme_names);
    registerStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_name, "<search_enzyme_name>", "Trypsin", "Name of the enzyme to be written to the pepXML file", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_search_enzyme_name, enzyme_names);

    // Cut after
    registerStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_cutafter, "<search_enzyme_cutafter>", "KR", "Residues after which the enzyme cuts (specified as a string of amino acids)", false , false);

    // No cut before
    registerStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_nocutbefore, "<search_enzyme_nocutbefore>", "P", "Residues that the enzyme will not cut before", false, false);

    // Number of enzyme termini
    registerStringOption_(TOPPMSFraggerAdapter::param_num_enzyme_termini, "<num_enzyme_termini>", "fully", "Number of enzyme termini (non-enzymatic (0), semi (1), fully (2)", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_num_enzyme_termini, ListUtils::create<String>("non-enzymatic,semi,fully"));

    // Allowed missed cleavages
    registerStringOption_(TOPPMSFraggerAdapter::param_allowed_missed_cleavage, "<allowed_missed_cleavage>", "2", "Allowed number of missed cleavages", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_allowed_missed_cleavage, zero_to_five); // 5 is the max. allowed value according to MSFragger

    // Digest min length
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_digest_min_length, "<digest_min_length>", 7, "Minimum length of peptides to be generated during in-silico digestion", false, false);

    // Digest max length
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_digest_max_length, "<digest_max_length>", 64, "Maximum length of peptides to be generated during in-silico digestion", false, false);

    // Digest min mass range
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_digest_mass_range_min, "<digest_mass_range_min>", 500.0, "Min mass of peptides to be generated (Da)", false, false);

    // Digest max mass range
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_digest_mass_range_max, "<digest_mass_range_max>", 5000.0, "Max mass of peptides to be generated (Da)", false, false);

    // TOPP varmod
    registerTOPPSubsection_("varmod", "Variable Modification Parameters");

    // Clip nterm M
    registerFlag_(TOPPMSFraggerAdapter::param_clip_nterm_m, "Specifies the trimming of a protein N-terminal methionine as a variable modification", false);

    // Modifications
    registerDoubleList_(TOPPMSFraggerAdapter::param_varmod_masses, "<varmod1_mass .. varmod7_mass>", emptyDoubles , "Masses for variable modifications", false, false);
    registerStringList_(TOPPMSFraggerAdapter::param_varmod_syntax, "<varmod1_syntax .. varmod7_syntax>", emptyStrings, "Syntax Strings for variable modifications", false, false);
    registerFlag_(TOPPMSFraggerAdapter::param_varmod_enable_common, "Enable common variable modifications (15.9949 M and 42.0106 [*)", false);

    // allow_multiple_variable_mods_on_residue
    registerFlag_(TOPPMSFraggerAdapter::param_not_allow_multiple_variable_mods_on_residue, "Do not allow any one amino acid to be modified by multiple variable modifications", false);

    // Max variable mods per mod
    registerStringOption_(TOPPMSFraggerAdapter::param_max_variable_mods_per_mod, "<max_variable_mods_per_mod>", "2", "Maximum number of residues that can be occupied by each variable modification", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_max_variable_mods_per_mod, zero_to_five);

    // Max variable mods combinations
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_max_variable_mods_combinations, "<max_variable_mods_combinations>", 5000, "Maximum allowed number of modified variably modified peptides from each peptide sequence, (maximum of 65534). If a greater number than the maximum is generated, only the unmodified peptide is considered", false, false);
    setMaxInt_(TOPPMSFraggerAdapter::param_max_variable_mods_combinations, 65534);

    // TOPP spectrum
    registerTOPPSubsection_("spectrum", "Spectrum Processing Parameters");

    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_minimum_peaks, "<minimum_peaks>", 10, "Minimum number of peaks in experimental spectrum for matching", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_use_topn_peaks, "<use_topN_peaks>", 50, "Pre-process experimental spectrum to only use top N peaks", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_minimum_ratio, "<minimum_ratio>", 0.0, "Filters out all peaks in experimental spectrum less intense than this multiple of the base peak intensity", false, false);
    setMaxFloat_(TOPPMSFraggerAdapter::param_minimum_ratio, 1.0);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_clear_mz_range_min, "<clear_mz_range_min>", 0.0, "Removes peaks in this m/z range prior to matching (minimum value). Useful for iTRAQ/TMT experiments (i.e. 0.0 150.0)", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_clear_mz_range_max, "<clear_mz_range_max>", 0.0, "Removes peaks in this m/z range prior to matching (maximum value). Useful for iTRAQ/TMT experiments (i.e. 0.0 150.0)", false, false);

    registerStringOption_(TOPPMSFraggerAdapter::param_max_fragment_charge, "<max_fragment_charge>", "2", "Maximum charge state for theoretical fragments to match", false, false);
    setValidStrings_(TOPPMSFraggerAdapter::param_max_fragment_charge, ListUtils::create<String>("1,2,3,4"));

    registerFlag_(TOPPMSFraggerAdapter::param_override_charge, "Ignores precursor charge and uses charge state specified in precursor_charge range (parameters: spectrum:precursor_charge_min and spectrum:precursor_charge_max)" , false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_precursor_charge_min, "<precursor_charge_min>", 1, "Min charge of precursor charge range to consider. If specified, also spectrum:override_charge must be set)" , false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_precursor_charge_max, "<precursor_charge_max>", 4, "Max charge of precursor charge range to consider. If specified, also spectrum:override_charge must be set)" , false, false);

    registerTOPPSubsection_("search", "Open Search Features");

    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_track_zero_topn, "<track_zero_topn>", 0, "Track top N unmodified peptide results separately from main results internally for boosting features. Should be set to a number greater than search:output_report_topN if zero bin boosting is desired", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_zero_bin_accept_expect, "<zero_bin_accept_expect>", 0.0, "Ranks a zero-bin hit above all non-zero-bin hit if it has expectation less than this value", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_zero_bin_mult_expect, "<zero_bin_mult_expect>", 1.0, "Multiplies expect value of PSMs in the zero-bin during results ordering (set to less than 1 for boosting)", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_add_topn_complementary, "<add_topn_complementary>", 0, "Inserts complementary ions corresponding to the top N most intense fragments in each experimental spectrum. Useful for recovery of modified peptides near C-terminus in open search. 0 disables this option", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_min_fragments_modeling, "<min_fragments_modeling>", 3, "Minimum number of matched peaks in PSM for inclusion in statistical modeling", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_min_matched_fragments, "<min_matched_fragments>", 4, "Minimum number of matched peaks for PSM to be reported. MSFragger recommends a minimum of 4 for narrow window searching and 6 for open searches", false, false);
    _registerNonNegativeInt(TOPPMSFraggerAdapter::param_output_report_topn, "<output_report_topn>", 1, "Reports top N PSMs per input spectrum", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_output_max_expect, "<output_max_expect>", 50.0, "Suppresses reporting of PSM if top hit has expectation greater than this threshold", false, false);

    registerTOPPSubsection_("statmod", "Static Modification Parameters");

    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_cterm_peptide, "<add_cterm_peptide>", 0.0, "Statically add mass in Da to C-terminal of peptide", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_nterm_peptide, "<add_nterm_peptide>", 0.0, "Statically add mass in Da to N-terminal of peptide", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_cterm_protein, "<add_cterm_protein>", 0.0, "Statically add mass in Da to C-terminal of protein", false, false);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_nterm_protein, "<add_nterm_protein>", 0.0, "Statically add mass in Da to N-terminal of protein", false, false);

    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_G_glycine,       "<add_G_glycine>",       0.0, "Statically add mass to glycine",       false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_A_alanine,       "<add_A_alanine>",       0.0, "Statically add mass to alanine",       false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_S_serine,        "<add_S_serine>",        0.0, "Statically add mass to serine",        false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_P_proline,       "<add_P_proline>",       0.0, "Statically add mass to proline",       false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_V_valine,        "<add_V_valine>",        0.0, "Statically add mass to valine",        false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_T_threonine,     "<add_T_threonine>",     0.0, "Statically add mass to threonine",     false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_C_cysteine,      "<add_C_cysteine>",      57.021464, "Statically add mass to cysteine",      false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_L_leucine,       "<add_L_leucine>",       0.0, "Statically add mass to leucine",       false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_I_isoleucine,    "<add_I_isoleucine>",    0.0, "Statically add mass to isoleucine",    false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_N_asparagine,    "<add_N_asparagine>",    0.0, "Statically add mass to asparagine",    false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_D_aspartic_acid, "<add_D_aspartic_acid>", 0.0, "Statically add mass to aspartic_acid", false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_Q_glutamine,     "<add_Q_glutamine>",     0.0, "Statically add mass to glutamine",     false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_K_lysine,        "<add_K_lysine>",        0.0, "Statically add mass to lysine",        false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_E_glutamic_acid, "<add_E_glutamic_acid>", 0.0, "Statically add mass to glutamic_acid", false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_M_methionine,    "<add_M_methionine>",    0.0, "Statically add mass to methionine",    false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_H_histidine,     "<add_H_histidine>",     0.0, "Statically add mass to histidine",     false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_F_phenylalanine, "<add_F_phenylalanine>", 0.0, "Statically add mass to phenylalanine", false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_R_arginine,      "<add_R_arginine>",      0.0, "Statically add mass to arginine",      false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_Y_tyrosine,      "<add_Y_tyrosine>",      0.0, "Statically add mass to tyrosine",      false, true);
    _registerNonNegativeDouble(TOPPMSFraggerAdapter::param_add_W_tryptophan,    "<add_W_tryptophan>",    0.0, "Statically add mass to tryptophan",    false, true);
  }


  ExitCodes main_(int, const char**) final override
  {
    try
    {
      // Java executable
      this->java_executable = this->getStringOption_(TOPPMSFraggerAdapter::param_java_executable);

      if (JavaInfo::canRun(this->java_executable, true) == false)
      {
        _fatalError("Java executable cannot be run!");
      }

      // Executable
      this->executable = this->getStringOption_(TOPPMSFraggerAdapter::param_executable);

      if (this->executable.empty())
      {
        // Looks for MSFRAGGER_PATH in the environment
        QProcessEnvironment env;
        QString qmsfragger_path = env.systemEnvironment().value("MSFRAGGER_PATH");
        if (qmsfragger_path.isEmpty())
        {
          _fatalError("No executable for MSFragger could be found!");
        }
        this->executable = qmsfragger_path;
      }

      // input, output, database name
      const String arg_database = this->getStringOption_(TOPPMSFraggerAdapter::param_database);
      const StringList arg_in = this->getStringList_(TOPPMSFraggerAdapter::param_in);
      checkUnique(arg_in, "Some input files were specified several times!");
      this->output_files = this->getStringList_(TOPPMSFraggerAdapter::param_out);
      checkUnique(this->output_files, "Out files must be unique for each input file!");

      // Check that none of the output files already exists
      for (const String & output_file : this->output_files)
      {
        if (File::exists(output_file))
        {
          _fatalError("Output file: " + output_file +  " already exists!");
        }
      }

      // Number of input files must match the number of output files (matched by index)
      if (arg_in.size() != this->output_files.size())
      {
        _fatalError("Number of output files has to match the number of input files!");
      }
      const FileTypes::Type out_file_type = FileHandler::getTypeByFileName(this->output_files[0]);

      // Check that all the output files have the same format (required by MSFragger)
      for (const String & output_file : this->output_files)
      {
        const FileTypes::Type current_file_type = FileHandler::getTypeByFileName(output_file);
        if (out_file_type != current_file_type)
        {
          _fatalError("Output files do not agree in format: "
              + FileTypes::typeToName(out_file_type)
          + " vs. "
          + FileTypes::typeToName(current_file_type));
        }
      }

      // Only pepXML and TSV currently supported as output file types.
      const bool out_is_pepxml = out_file_type == FileTypes::PEPXML;
      if (out_is_pepxml == false && out_file_type != FileTypes::TSV)
      {
        _fatalError("Output file type (determined from file extension) invalid. Choose either pepXML or tsv!");
      }
      this->file_type = out_is_pepxml ? "pepXML" : "tsv";

      // tolerance
      const double arg_precursor_mass_tolerance(this->getDoubleOption_(TOPPMSFraggerAdapter::param_precursor_mass_tolerance));
      const String & arg_precursor_mass_unit = this->getStringOption_(TOPPMSFraggerAdapter::param_precursor_mass_unit);
      const double arg_precursor_true_tolerance(this->getDoubleOption_(TOPPMSFraggerAdapter::param_precursor_true_tolerance));
      const String & arg_precursor_true_unit = this->getStringOption_(TOPPMSFraggerAdapter::param_precursor_true_unit);
      const double arg_fragment_mass_tolerance(this->getDoubleOption_(TOPPMSFraggerAdapter::param_fragment_mass_tolerance));
      const String & arg_fragment_mass_unit = this->getStringOption_(TOPPMSFraggerAdapter::param_fragment_mass_unit);
      const String & arg_isotope_error = this->getStringOption_(TOPPMSFraggerAdapter::param_isotope_error);

      // digest
      const String & arg_search_enzyme_name = this->getStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_name);
      const String & arg_search_enzyme_cutafter = this->getStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_cutafter);
      const String & arg_search_enzyme_nocutbefore = this->getStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_nocutbefore);

      std::map< String,int > num_enzyme_termini;
      num_enzyme_termini["non-enzymatic"] = 0;
      num_enzyme_termini["semi"] = 1;
      num_enzyme_termini["fully"] = 2;
      const int arg_num_enzyme_termini = num_enzyme_termini[this->getStringOption_(TOPPMSFraggerAdapter::param_num_enzyme_termini)];

      const String & arg_allowed_missed_cleavage = this->getStringOption_(TOPPMSFraggerAdapter::param_allowed_missed_cleavage);
      const int arg_digest_min_length = this->getIntOption_(TOPPMSFraggerAdapter::param_digest_min_length);
      const int arg_digest_max_length = this->getIntOption_(TOPPMSFraggerAdapter::param_digest_max_length);
      ensureRange(arg_digest_min_length, arg_digest_max_length, "Maximum length of digest is not allowed to be smaller than minimum length of digest");

      const double arg_digest_mass_range_min = this->getDoubleOption_(TOPPMSFraggerAdapter::param_digest_mass_range_min);
      const double arg_digest_mass_range_max = this->getDoubleOption_(TOPPMSFraggerAdapter::param_digest_mass_range_max);
      ensureRange(arg_digest_mass_range_min, arg_digest_mass_range_max, "Maximum digest mass is not allowed to be smaller than minimum digest mass!");

      // varmod
      const bool arg_clip_nterm_m = this->getFlag_(TOPPMSFraggerAdapter::param_clip_nterm_m);
      std::vector< double > arg_varmod_masses = this->getDoubleList_(TOPPMSFraggerAdapter::param_varmod_masses);
      std::vector< String > arg_varmod_syntax = this->getStringList_(TOPPMSFraggerAdapter::param_varmod_syntax);

      // Assignment of mass to syntax is by index, so the vectors have to be the same length
      if (arg_varmod_masses.size() != arg_varmod_syntax.size())
      {
        _fatalError("List of arguments for the parameters 'varmod_masses' and 'varmod_syntax' must have the same length!");
      }
      // Only up to 7 variable modifications are allowed
      if (arg_varmod_masses.size() > 7)
      {
        _fatalError("MSFragger is restricted to at most 7 variable modifications.");
      }

      // Add common variable modifications if requested
      if (this->getFlag_(TOPPMSFraggerAdapter::param_varmod_enable_common))
      {
        // Oxidation on methionine
        this->_addVarMod(arg_varmod_masses, arg_varmod_syntax, 15.9949, "M");

        // N-terminal acetylation
        this->_addVarMod(arg_varmod_masses, arg_varmod_syntax, 42.0106, "[*");
      }

      const bool arg_not_allow_multiple_variable_mods_on_residue = this->getFlag_(TOPPMSFraggerAdapter::param_not_allow_multiple_variable_mods_on_residue);
      const String & arg_max_variable_mods_per_mod  = this->getStringOption_(TOPPMSFraggerAdapter::param_max_variable_mods_per_mod);
      const int arg_max_variable_mods_combinations = this->getIntOption_(TOPPMSFraggerAdapter::param_max_variable_mods_combinations);

      // spectrum
      const int arg_minimum_peaks = this->getIntOption_(TOPPMSFraggerAdapter::param_minimum_peaks);
      const int arg_use_topn_peaks  = this->getIntOption_(TOPPMSFraggerAdapter::param_use_topn_peaks);
      const double arg_minimum_ratio = this->getDoubleOption_(TOPPMSFraggerAdapter::param_minimum_ratio);
      const double arg_clear_mz_range_min = this->getDoubleOption_(TOPPMSFraggerAdapter::param_clear_mz_range_min);
      const double arg_clear_mz_range_max = this->getDoubleOption_(TOPPMSFraggerAdapter::param_clear_mz_range_max);
      ensureRange(arg_clear_mz_range_min, arg_clear_mz_range_max, "Maximum clear mz value is not allowed to be smaller than minimum clear mz value!");
      const String & arg_max_fragment_charge = this->getStringOption_(TOPPMSFraggerAdapter::param_max_fragment_charge);
      const bool arg_override_charge = this->getFlag_(TOPPMSFraggerAdapter::param_override_charge);
      const int arg_precursor_charge_min = this->getIntOption_(TOPPMSFraggerAdapter::param_precursor_charge_min);
      const int arg_precursor_charge_max = this->getIntOption_(TOPPMSFraggerAdapter::param_precursor_charge_max);
      ensureRange(arg_precursor_charge_min, arg_precursor_charge_max, "Maximum precursor charge is not allowed to be smaller than minimum precursor charge!");

      // Ensures that the user is aware of overriding the precursoe charges
      if ((arg_precursor_charge_min != 1 || arg_precursor_charge_max != 4) && arg_override_charge == false)
      {
        _fatalError("If you want to ignore the precursor charge, please also set the -" + TOPPMSFraggerAdapter::param_override_charge + " flag!");
      }

      // search
      const int arg_track_zero_topn = this->getIntOption_(TOPPMSFraggerAdapter::param_track_zero_topn);
      const double arg_zero_bin_accept_expect = this->getDoubleOption_(TOPPMSFraggerAdapter::param_zero_bin_accept_expect);
      const double arg_zero_bin_mult_expect = this->getDoubleOption_(TOPPMSFraggerAdapter::param_zero_bin_mult_expect);
      const int arg_add_topn_complementary = this->getIntOption_(TOPPMSFraggerAdapter::param_add_topn_complementary);
      const int arg_min_fragments_modeling = this->getIntOption_(TOPPMSFraggerAdapter::param_min_fragments_modeling);
      const int arg_min_matched_fragments = this->getIntOption_(TOPPMSFraggerAdapter::param_min_matched_fragments);
      const int arg_output_report_topn = this->getIntOption_(TOPPMSFraggerAdapter::param_output_report_topn);
      const double arg_output_max_expect = this->getDoubleOption_(TOPPMSFraggerAdapter::param_output_max_expect);

      // statmod
      const double arg_add_cterm_peptide = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_cterm_peptide);
      const double arg_add_nterm_peptide = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_nterm_peptide);
      const double arg_add_cterm_protein = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_cterm_protein);
      const double arg_add_nterm_protein = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_nterm_protein);
      const double arg_add_G_glycine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_G_glycine);
      const double arg_add_A_alanine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_A_alanine);
      const double arg_add_S_serine      = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_S_serine);
      const double arg_add_P_proline     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_P_proline);
      const double arg_add_V_valine      = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_V_valine);
      const double arg_add_T_threonine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_T_threonine);
      const double arg_add_C_cysteine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_C_cysteine);
      const double arg_add_L_leucine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_L_leucine);
      const double arg_add_I_isoleucine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_I_isoleucine);
      const double arg_add_N_asparagine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_N_asparagine);
      const double arg_add_D_aspartic_acid     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_D_aspartic_acid);
      const double arg_add_Q_glutamine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_Q_glutamine);
      const double arg_add_K_lysine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_K_lysine);
      const double arg_add_E_glutamic_acid     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_E_glutamic_acid);
      const double arg_add_M_methionine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_M_methionine);
      const double arg_add_H_histidine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_H_histidine);
      const double arg_add_F_phenylalanine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_F_phenylalanine);
      const double arg_add_R_arginine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_R_arginine);
      const double arg_add_Y_tyrosine     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_Y_tyrosine);
      const double arg_add_W_tryptophan     = this->getDoubleOption_(TOPPMSFraggerAdapter::param_add_W_tryptophan);

      // Parameters have been read in and verified, they are now going to be written into the fragger.params file in a temporary directory
      this->working_directory = this->makeTempDirectory_().toQString();
      const QFileInfo tmp_param_file(this->working_directory, "fragger.params");
      this->parameter_file_path =  String(tmp_param_file.absoluteFilePath());

      // Create a link to the FASTA database in the temporary directory (since MSFragger needs to construct the pepindex)
      // in the temp directory
      const QFileInfo database_fileinfo(arg_database.toQString());
      const QString database_path = QFileInfo(this->working_directory, database_fileinfo.fileName()).absoluteFilePath();
      if (QFile::link(database_fileinfo.absoluteFilePath(), database_path) == false)
      {
        _fatalError("Could not create link to database in tmp directory: " + String(this->working_directory));
      }

      // Create link to all the input files in the temp directory
      for (const String & input_file : arg_in)
      {
        const QFileInfo file_info(input_file.toQString());
        const QString link_name = QFileInfo(this->working_directory, file_info.fileName()).absoluteFilePath();

        if (QFile::link(file_info.absoluteFilePath(), link_name) == false)
        {
          _fatalError("Could not create link to input file in tmp directory: " + String(this->working_directory));
        }
        this->input_files.push_back(link_name);
      }

      writeDebug_("Parameter file for MSFragger: '" + this->parameter_file_path + "'", TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE);
      writeDebug_("Working Directory: '" + String(this->working_directory) + "'", TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE);
      writeDebug_("If you want to keep the working directory and the parameter file, set the -debug to 2", 1);
      ofstream os(this->parameter_file_path.c_str());

      // Write all the parameters into the file
      os << "database_name = " << String(database_path)
                               << "\nnum_threads = " << this->getIntOption_("threads")
                               << "\n\nprecursor_mass_tolerance = " << arg_precursor_mass_tolerance
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

      // Write variable modifications (and also write to log)
      writeLog_("Variable Modifications set to:");
      for (Size i = 0; i < arg_varmod_masses.size(); ++i)
      {
        const String varmod = "variable_mod_0" + String(i+1) + " = " + String(arg_varmod_masses[i]) + " " + String(arg_varmod_syntax[i]);
        os << "\n" << varmod;
        writeLog_(varmod);
      }

      os << std::endl
          << "\nallow_multiple_variable_mods_on_residue = " << (arg_not_allow_multiple_variable_mods_on_residue ? 0 : 1)
          << "\nmax_variable_mods_per_mod = " << arg_max_variable_mods_per_mod
          << "\nmax_variable_mods_combinations = " << arg_max_variable_mods_combinations
          << "\n\noutput_file_extension = " << this->file_type
          << "\noutput_format = " << this->file_type
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
    process_params << "-Xmx" + QString::number(this->getIntOption_(TOPPMSFraggerAdapter::param_java_heapmemory)) + "m"
        << "-jar" << this->executable.toQString()
        << this->parameter_file_path.toQString();

    for (const QString & input_file : this->input_files)
    {
      process_params << input_file;
    }

    QProcess process_msfragger;
    process_msfragger.setWorkingDirectory(this->working_directory);

    if (this->debug_level_ >= TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE)
    {
      writeDebug_("COMMAND LINE CALL IS:", 1);
      String command_line = this->java_executable;
      for (const String & process_param : process_params)
      {
        command_line += (" " + process_param);
      }
      writeDebug_(command_line, TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE);
    }

    process_msfragger.start(this->java_executable.toQString(), process_params);

    if (process_msfragger.waitForFinished(-1) == false || process_msfragger.exitCode() != 0)
    {
      LOG_FATAL_ERROR << "FATAL: Invokation of MSFraggerAdapter has failed. Error code was: " + process_msfragger.exitCode() << std::endl;
      return EXTERNAL_PROGRAM_ERROR;
    }

    // Copy the output files of MSFragger to the user location
    for (int i = 0; i < this->input_files.size(); ++i)
    {
      QFile::copy(
          (File::removeExtension(this->input_files[i]) + "." + this->file_type).toQString(),
          QFileInfo(this->output_files[i].toQString()).absoluteFilePath());
    }

    return EXECUTION_OK;
  }


private:

  QString working_directory;

  String java_executable;
  String executable;

  String parameter_file_path;
  QStringList input_files;
  StringList output_files;

  String file_type;


  // Adds variable modification if not already present
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
    LOG_FATAL_ERROR << "FATAL: " << message << std::endl;
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
      LOG_ERROR << "FATAL: " << message << std::endl;
      throw 1;
    }
  }
};

const String TOPPMSFraggerAdapter::param_java_executable = "java_executable";
const String TOPPMSFraggerAdapter::param_java_heapmemory = "java_heapmemory";
const String TOPPMSFraggerAdapter::param_executable = "executable";
const String TOPPMSFraggerAdapter::param_in = "in";
const String TOPPMSFraggerAdapter::param_out = "out";
const String TOPPMSFraggerAdapter::param_database = "database";

// tolerance
const String TOPPMSFraggerAdapter::param_precursor_mass_tolerance = "tolerance:precursor_mass_tolerance";
const String TOPPMSFraggerAdapter::param_precursor_mass_unit = "tolerance:precursor_mass_unit";
const String TOPPMSFraggerAdapter::param_precursor_true_tolerance = "tolerance:precursor_true_tolerance";
const String TOPPMSFraggerAdapter::param_precursor_true_unit = "tolerance:precursor_true_unit";
const String TOPPMSFraggerAdapter::param_fragment_mass_tolerance = "tolerance:fragment_mass_tolerance";
const String TOPPMSFraggerAdapter::param_fragment_mass_unit = "tolerance:fragment_mass_unit";
const String TOPPMSFraggerAdapter::param_isotope_error = "tolerance:isotope_error";

// digest
const String TOPPMSFraggerAdapter::param_search_enzyme_name = "digest:search_enzyme_name";
const String TOPPMSFraggerAdapter::param_search_enzyme_cutafter = "digest:search_enzyme_cutafter";
const String TOPPMSFraggerAdapter::param_search_enzyme_nocutbefore = "digest:search_enzyme_nocutbefore";
const String TOPPMSFraggerAdapter::param_num_enzyme_termini = "digest:num_enzyme_termini";
const String TOPPMSFraggerAdapter::param_allowed_missed_cleavage = "digest:allowed_missed_cleavage";
const String TOPPMSFraggerAdapter::param_digest_min_length = "digest:min_length";
const String TOPPMSFraggerAdapter::param_digest_max_length = "digest:max_length";
const String TOPPMSFraggerAdapter::param_digest_mass_range_min = "digest:mass_range_min";
const String TOPPMSFraggerAdapter::param_digest_mass_range_max = "digest:mass_range_max";

// varmod
const String TOPPMSFraggerAdapter::param_clip_nterm_m = "varmod:clip_nterm_m";
const String TOPPMSFraggerAdapter::param_varmod_masses = "varmod:masses";
const String TOPPMSFraggerAdapter::param_varmod_syntax = "varmod:syntaxes";
const String TOPPMSFraggerAdapter::param_varmod_enable_common = "varmod:enable_common";
const String TOPPMSFraggerAdapter::param_not_allow_multiple_variable_mods_on_residue = "varmod:not_allow_multiple_variable_mods_on_residue";
const String TOPPMSFraggerAdapter::param_max_variable_mods_per_mod = "varmod:max_variable_mods_per_mod";
const String TOPPMSFraggerAdapter::param_max_variable_mods_combinations = "varmod:max_variable_mods_combinations";

// spectrum
const String TOPPMSFraggerAdapter::param_minimum_peaks = "spectrum:minimum_peaks";
const String TOPPMSFraggerAdapter::param_use_topn_peaks = "spectrum:use_topn_peaks";
const String TOPPMSFraggerAdapter::param_minimum_ratio = "spectrum:minimum_ratio";
const String TOPPMSFraggerAdapter::param_clear_mz_range_min = "spectrum:clear_mz_range_min";
const String TOPPMSFraggerAdapter::param_clear_mz_range_max = "spectrum:clear_mz_range_max";
const String TOPPMSFraggerAdapter::param_max_fragment_charge = "spectrum:max_fragment_charge";
const String TOPPMSFraggerAdapter::param_override_charge = "spectrum:override_charge";
const String TOPPMSFraggerAdapter::param_precursor_charge_min = "spectrum:precursor_charge_min";
const String TOPPMSFraggerAdapter::param_precursor_charge_max = "spectrum:precursor_charge_max";

// search
const String TOPPMSFraggerAdapter::param_track_zero_topn = "search:track_zero_topn";
const String TOPPMSFraggerAdapter::param_zero_bin_accept_expect = "search:zero_bin_accept_expect";
const String TOPPMSFraggerAdapter::param_zero_bin_mult_expect = "search:zero_bin_mult_expect";
const String TOPPMSFraggerAdapter::param_add_topn_complementary = "search:add_topn_complementary";
const String TOPPMSFraggerAdapter::param_min_fragments_modeling = "search:min_fragments_modeling";
const String TOPPMSFraggerAdapter::param_min_matched_fragments = "search:min_matched_fragments";
const String TOPPMSFraggerAdapter::param_output_report_topn = "search:output_report_topn";
const String TOPPMSFraggerAdapter::param_output_max_expect = "search:output_max_expect";

// statmod
const String TOPPMSFraggerAdapter::param_add_cterm_peptide = "statmod:add_cterm_peptide";
const String TOPPMSFraggerAdapter::param_add_nterm_peptide = "statmod:add_nterm_peptide";
const String TOPPMSFraggerAdapter::param_add_cterm_protein = "statmod:add_cterm_protein";
const String TOPPMSFraggerAdapter::param_add_nterm_protein = "statmod:add_nterm_protein";
const String TOPPMSFraggerAdapter::param_add_G_glycine = "statmod:add_G_glycine";
const String TOPPMSFraggerAdapter::param_add_A_alanine = "statmod:add_A_alanine";
const String TOPPMSFraggerAdapter::param_add_S_serine = "statmod:add_S_serine";
const String TOPPMSFraggerAdapter::param_add_P_proline = "statmod:add_P_proline";
const String TOPPMSFraggerAdapter::param_add_V_valine = "statmod:add_V_valine";
const String TOPPMSFraggerAdapter::param_add_T_threonine = "statmod:add_T_threonine";
const String TOPPMSFraggerAdapter::param_add_C_cysteine = "statmod:add_C_cysteine";
const String TOPPMSFraggerAdapter::param_add_L_leucine = "statmod:add_L_leucine";
const String TOPPMSFraggerAdapter::param_add_I_isoleucine = "statmod:add_I_isoleucine";
const String TOPPMSFraggerAdapter::param_add_N_asparagine = "statmod:add_N_asparagine";
const String TOPPMSFraggerAdapter::param_add_D_aspartic_acid = "statmod:add_D_aspartic_acid";
const String TOPPMSFraggerAdapter::param_add_Q_glutamine = "statmod:add_Q_glutamine";
const String TOPPMSFraggerAdapter::param_add_K_lysine = "statmod:add_K_lysine";
const String TOPPMSFraggerAdapter::param_add_E_glutamic_acid = "statmod:add_E_glutamic_acid";
const String TOPPMSFraggerAdapter::param_add_M_methionine = "statmod:add_M_methionine";
const String TOPPMSFraggerAdapter::param_add_H_histidine = "statmod:add_H_histidine";
const String TOPPMSFraggerAdapter::param_add_F_phenylalanine = "statmod:add_F_phenylalanine";
const String TOPPMSFraggerAdapter::param_add_R_arginine = "statmod:add_R_arginine";
const String TOPPMSFraggerAdapter::param_add_Y_tyrosine = "statmod:add_Y_tyrosine";
const String TOPPMSFraggerAdapter::param_add_W_tryptophan = "statmod:add_W_tryptophan";

const int TOPPMSFraggerAdapter::LOG_LEVEL_VERBOSE = 1;

int main(int argc, const char** argv)
{
  TOPPMSFraggerAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
