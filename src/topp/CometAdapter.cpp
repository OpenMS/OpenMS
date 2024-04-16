// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Leon Bichmann, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
// TODO remove this once we have handler transform support
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>

#include <QStringList>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_CometAdapter CometAdapter

@brief Identifies peptides in MS/MS spectra via Comet.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; CometAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

@em Comet must be installed/downloaded before this wrapper can be used. OpenMS installers ship with Comet.

@warning We recommend to use Comet 2019.01 rev. 5 or later, due to a serious "empty result" bug in earlier versions (which occurs frequently on Windows; Linux seems not/less affected).

Comet settings not exposed by this adapter can be directly adjusted using a param file, which can be generated using comet -p.
By default, All (!) parameters available explicitly via this param file will take precedence over the wrapper parameters.

Parameter names have been changed to match names found in other search engine adapters, however some are Comet specific.
For a detailed description of all available parameters check the Comet documentation at http://comet-ms.sourceforge.net/parameters/parameters_201601/
The default parameters are set for a high resolution instrument.

@note This adapter supports 15N labeling by specifying the 20 AA modifications 'Label:15N(x)' as fixed modifications.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_CometAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_CometAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPCometAdapter :
  public SearchEngineBase
{
public:
  TOPPCometAdapter() :
    SearchEngineBase("CometAdapter", "Annotates MS/MS spectra using Comet.", true,
             {
                 {"Eng, Jimmy K. and Jahan, Tahmina A. and Hoopmann, Michael R.",
                 "Comet: An open-source MS/MS sequence database search tool",
                 "PROTEOMICS 2013; 13-1: 22--24",
                 "10.1002/pmic.201200439"}
             })
  {
  }

protected:

  map<string,int> num_enzyme_termini {{"semi",1},{"fully",2},{"C-term unspecific", 8},{"N-term unspecific",9}};

  void registerOptionsAndFlags_() override
  {

    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", { "mzML" } );
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", { "idXML"} );
    registerInputFile_("database", "<file>", "", "FASTA file", true, false, {"skipexists"});
    setValidFormats_("database", { "FASTA" } );
    registerInputFile_("comet_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      "comet.exe", // this is the name on ALL platforms currently...
      "The Comet executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", true, false, {"is_executable"});

    //
    // Optional parameters
    //

    //Files
    registerOutputFile_("pin_out", "<file>", "", "Output file - for Percolator input", false);
    setValidFormats_("pin_out", ListUtils::create<String>("tsv"));
    registerInputFile_("default_params_file", "<file>", "", "Default Comet params file. All parameters of this take precedence. A template file can be generated using 'comet.exe -p'", false, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("default_params_file", ListUtils::create<String>("txt"));

    //Masses
    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 10.0, "Precursor monoisotopic mass tolerance (Comet parameter: peptide_mass_tolerance).  See also precursor_error_units to set the unit.",false);
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "Unit of precursor monoisotopic mass tolerance for parameter precursor_mass_tolerance (Comet parameter: peptide_mass_units)", false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("amu,ppm,Da"));
    //registerIntOption_("mass_type_parent", "<num>", 1, "0=average masses, 1=monoisotopic masses", false, true);
    //registerIntOption_("mass_type_fragment", "<num>", 1, "0=average masses, 1=monoisotopic masses", false, true);
    //registerIntOption_("precursor_tolerance_type", "<num>", 0, "0=average masses, 1=monoisotopic masses", false, false);
    registerStringOption_(Constants::UserParam::ISOTOPE_ERROR, "<choice>", "off", "This parameter controls whether the peptide_mass_tolerance takes into account possible isotope errors in the precursor mass measurement. Use -8/-4/0/4/8 only for SILAC.", false, false);
    setValidStrings_(Constants::UserParam::ISOTOPE_ERROR, ListUtils::create<String>("off,0/1,0/1/2,0/1/2/3,-8/-4/0/4/8,-1/0/1/2/3"));

    //Fragment Ions
    registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.01,
                          "This is half the bin size, which is used to segment the MS/MS spectrum. Thus, the value should be a bit higher than for other search engines, since the bin might not be centered around the peak apex (see 'fragment_bin_offset')."
                          "CAUTION: Low tolerances have heavy impact on RAM usage (since Comet uses a lot of bins in this case). Consider using use_sparse_matrix and/or spectrum_batch_size.", false);
    setMinFloat_("fragment_mass_tolerance", 0.0001);
    
    registerStringOption_("fragment_error_units", "<unit>", "Da", "Fragment monoisotopic mass error units", false);
    setValidStrings_("fragment_error_units", { "Da" }); // only Da allowed
    
    registerDoubleOption_("fragment_bin_offset", "<fraction>", 0.0, "Offset of fragment bins. Recommended by Comet: low-res: 0.4, high-res: 0.0", false);
    setMinFloat_("fragment_bin_offset", 0.0);
    setMaxFloat_("fragment_bin_offset", 1.0);

    registerStringOption_("instrument", "<choice>", "high_res", "Comets theoretical_fragment_ions parameter: theoretical fragment ion peak representation, high-res: sum of intensities plus flanking bins, ion trap (low-res) ms/ms: sum of intensities of central M bin only", false);
    setValidStrings_("instrument", ListUtils::create<String>("low_res,high_res"));
    registerStringOption_("use_A_ions", "<num>", "false", "use A ions for PSM", false, true);
    setValidStrings_("use_A_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_B_ions", "<num>", "true", "use B ions for PSM", false, true);
    setValidStrings_("use_B_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_C_ions", "<num>", "false", "use C ions for PSM", false, true);
    setValidStrings_("use_C_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_X_ions", "<num>", "false", "use X ions for PSM", false, true);
    setValidStrings_("use_X_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_Y_ions", "<num>", "true", "use Y ions for PSM", false, true);
    setValidStrings_("use_Y_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_Z_ions", "<num>", "false", "use Z ions for PSM", false, true);
    setValidStrings_("use_Z_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_NL_ions", "<num>", "false", "use neutral loss (NH3, H2O) ions from b/y for PSM", false, true);
    setValidStrings_("use_NL_ions", ListUtils::create<String>("true,false"));

    //Search Enzyme
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllCometNames(all_enzymes);
    registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false, false);
    setValidStrings_("enzyme", all_enzymes);
    registerStringOption_("second_enzyme", "<cleavage site>", "", "Additional enzyme used for peptide digestion.", false, true);
    setValidStrings_("second_enzyme", all_enzymes);

    registerStringOption_("num_enzyme_termini", "<choice>", "fully", "Specify the termini where the cleavage rule has to match", false, false);
    setValidStrings_("num_enzyme_termini", { "semi", "fully", "C-term unspecific", "N-term unspecific" } );
    registerIntOption_("missed_cleavages", "<num>", 1, "Number of possible cleavage sites missed by the enzyme. It has no effect if enzyme is unspecific cleavage.", false, false);
    setMinInt_("missed_cleavages", 0);
    setMaxInt_("missed_cleavages", 5);

    registerIntOption_("min_peptide_length", "<num>", 5, "Minimum peptide length to consider.", false);
    setMinInt_("min_peptide_length", 5);
    setMaxInt_("min_peptide_length", 63);
    registerIntOption_("max_peptide_length", "<num>", 63, "Maximum peptide length to consider.", false);
    setMinInt_("max_peptide_length", 5);
    setMaxInt_("max_peptide_length", 63);

    //Output
    registerIntOption_("num_hits", "<num>", 1, "Number of peptide hits (PSMs) per spectrum in output file", false, false);

    //mzXML/mzML parameters
    registerStringOption_("precursor_charge", "[min]:[max]", "0:0", "Precursor charge range to search (if spectrum is not annotated with a charge or if override_charge!=keep any known): 0:[num] == search all charges, 2:6 == from +2 to +6, 3:3 == +3", false, false);
    registerStringOption_("override_charge", "<choice>", "keep known search unknown", "_keep any known_: keep any precursor charge state (from input), _ignore known_: ignore known precursor charge state and use precursor_charge parameter, _ignore outside range_: ignore precursor charges outside precursor_charge range, _keep known search unknown_: keep any known precursor charge state. For unknown charge states, search as singly charged if there is no signal above the precursor m/z or use the precursor_charge range", false, false);
    setValidStrings_("override_charge", ListUtils::create<String>("keep any known,ignore known,ignore outside range,keep known search unknown"));
    registerIntOption_("ms_level", "<num>", 2, "MS level to analyze, valid are levels 2 (default) or 3", false, false);
    setMinInt_("ms_level", 2);
    setMaxInt_("ms_level", 3);
    registerStringOption_("activation_method", "<method>", "ALL", "If not ALL, only searches spectra of the given method", false, false);
    setValidStrings_("activation_method", ListUtils::create<String>("ALL,CID,ECD,ETD,PQD,HCD,IRMPD"));

    //Misc. parameters
    //scan range
    registerStringOption_("digest_mass_range", "[min]:[max]", "600:5000", "MH+ peptide mass range to analyze", false, true);
    registerIntOption_("max_fragment_charge", "<posnum>", 3, "Set maximum fragment charge state to analyze as long as still lower than precursor charge - 1. (Allowed max 5)", false, false);
    setMinInt_("max_fragment_charge", 1);
    setMaxInt_("max_fragment_charge", 5);
    registerIntOption_("max_precursor_charge", "<posnum>", 5, "set maximum precursor charge state to analyze (allowed max 9)", false, true);
    setMinInt_("max_precursor_charge", 1);
    setMaxInt_("max_precursor_charge", 9);
    registerStringOption_("clip_nterm_methionine", "<bool>", "false", "If set to true, also considers the peptide sequence w/o N-term methionine separately and applies appropriate N-term mods to it", false, false);
    setValidStrings_("clip_nterm_methionine", ListUtils::create<String>("true,false"));
    registerIntOption_("spectrum_batch_size", "<posnum>", 20000, "max. number of spectra to search at a time; use 0 to search the entire scan range in one batch", false, true);
    setMinInt_("spectrum_batch_size", 0);
    registerDoubleList_("mass_offsets", "<doubleoffset1, doubleoffset2,...>", {0.0}, "One or more mass offsets to search (values subtracted from deconvoluted precursor mass). Has to include 0.0 if you want the default mass to be searched.", false, true);

    // spectral processing
    registerIntOption_("minimum_peaks", "<posnum>", 10, "Required minimum number of peaks in spectrum to search (default 10)", false, true);
    registerDoubleOption_("minimum_intensity", "<posfloat>", 0.0, "Minimum intensity value to read in", false, true);
    setMinFloat_("minimum_intensity", 0.0);
    registerStringOption_("remove_precursor_peak", "<choice>", "no", "no = no removal, yes = remove all peaks around precursor m/z, charge_reduced = remove all charge reduced precursor peaks (for ETD/ECD). phosphate_loss = remove the HPO3 (-80) and H3PO4 (-98) precursor phosphate neutral loss peaks. See also remove_precursor_tolerance", false, true);
    setValidStrings_("remove_precursor_peak", ListUtils::create<String>("no,yes,charge_reduced,phosphate_loss"));
    registerDoubleOption_("remove_precursor_tolerance", "<posfloat>", 1.5, "one-sided tolerance for precursor removal in Thompson", false, true);
    registerStringOption_("clear_mz_range", "[minfloatmz]:[maxfloatmz]", "0:0", "for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range, if not 0:0", false, true);

    //Modifications
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>("Carbamidomethyl (C)", ','), "Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>("Oxidation (M)", ','), "Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    registerIntList_("binary_modifications", "<mods>", {}, 
        "List of modification group indices. Indices correspond to the binary modification index used by comet to group individually searched lists of variable modifications.\n" 
        "Note: if set, both variable_modifications and binary_modifications need to have the same number of entries as the N-th entry corresponds to the N-th variable_modification.\n"
        "      if left empty (default), all entries are internally set to 0 generating all permutations of modified and unmodified residues.\n"
        "      For a detailed explanation please see the parameter description in the Comet help.",
        false);

    registerIntOption_("max_variable_mods_in_peptide", "<num>", 5, "Set a maximum number of variable modifications per peptide", false, true);
    registerStringOption_("require_variable_mod", "<bool>", "false", "If true, requires at least one variable modification per peptide", false, true);
    setValidStrings_("require_variable_mod", ListUtils::create<String>("true,false"));

    // register peptide indexing parameter (with defaults for this search engine) TODO: check if search engine defaults are needed
    registerPeptideIndexingParameter_(PeptideIndexing().getParameters()); 
  }

  const vector<const ResidueModification*> getModifications_(const StringList& modNames)
  {
    vector<const ResidueModification*> modifications;

    // iterate over modification names and add to vector
    for (const auto& modification : modNames)
    {
      if (modNames.empty())
      {
        continue;
      }
      modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
    }

    return modifications;
  }

  ExitCodes createParamFile_(ostream& os, const String& comet_version)
  {
    os << comet_version << "\n";              // required as first line in the param file
    os << "# Comet MS/MS search engine parameters file.\n";
    os << "# Everything following the '#' symbol is treated as a comment.\n";
    os << "database_name = " << getStringOption_("database") << "\n";
    os << "decoy_search = " << 0 << "\n";                                               // 0=no (default), 1=concatenated search, 2=separate search
    os << "peff_format = 0\n";                                                          // 0=no (normal fasta, default), 1=PEFF PSI-MOD, 2=PEFF Unimod
    os << "peff_obo =\n";                                                               // path to PSI Mod or Unimod OBO file

    os << "num_threads = " << getIntOption_("threads") << "\n";                         // 0=poll CPU to set num threads; else specify num threads directly (max 64)

    // masses
    map<String,int> precursor_error_units;
    precursor_error_units["amu"] = 0;
    precursor_error_units["mmu"] = 1;
    precursor_error_units["ppm"] = 2;

    map<string,int> isotope_error;
    isotope_error["off"] = 0;
    isotope_error["0/1"] = 1;
    isotope_error["0/1/2"] = 2;
    isotope_error["0/1/2/3"] = 3;
    isotope_error["-8/-4/0/4/8"] = 4;
    isotope_error["-1/0/1/2/3"] = 5;

    os << "peptide_mass_tolerance = " << getDoubleOption_("precursor_mass_tolerance") << "\n";
    os << "peptide_mass_units = " << precursor_error_units[getStringOption_("precursor_error_units")] << "\n";                  // 0=amu, 1=mmu, 2=ppm
    os << "mass_type_parent = " << 1 << "\n";                    // 0=average masses, 1=monoisotopic masses
    os << "mass_type_fragment = " << 1 << "\n";                  // 0=average masses, 1=monoisotopic masses
    os << "precursor_tolerance_type = " << 1 << "\n";            // 0=MH+ (default), 1=precursor m/z; only valid for amu/mmu tolerances
    os << "isotope_error = " << isotope_error[getStringOption_(Constants::UserParam::ISOTOPE_ERROR)] << "\n";                   // 0=off, 1=0/1 (C13 error), 2=0/1/2, 3=0/1/2/3, 4=-8/-4/0/4/8 (for +4/+8 labeling)

    // search enzyme

    String enzyme_name = getStringOption_("enzyme");
    String enzyme_number = String(ProteaseDB::getInstance()->getEnzyme(enzyme_name)->getCometID());
    String second_enzyme_name = getStringOption_("second_enzyme");
    String enzyme2_number = "0";
    if (!second_enzyme_name.empty())
    {
      enzyme2_number = String(ProteaseDB::getInstance()->getEnzyme(second_enzyme_name)->getCometID());
    }

    os << "search_enzyme_number = " << enzyme_number << "\n";                // choose from list at end of this params file
    os << "search_enzyme2_number = " << enzyme2_number << "\n";              // second enzyme; set to 0 if no second enzyme
    os << "num_enzyme_termini = " << num_enzyme_termini[getStringOption_("num_enzyme_termini")] << "\n"; // 1 (semi-digested), 2 (fully digested, default), 8 C-term unspecific , 9 N-term unspecific
    os << "allowed_missed_cleavage = " << getIntOption_("missed_cleavages") << "\n";             // maximum value is 5; for enzyme search

    // Up to 9 variable modifications are supported
    // # format:  <mass> <residues> <0=variable/else binary> <max_mods_per_peptide> <term_distance> <n/c-term> <required> <neutral_loss>
    //     e.g. 79.966331 STY 0 3 -1 0 0 97.976896
    vector<String> variable_modifications_names = getStringList_("variable_modifications");
    const vector<const ResidueModification*> variable_modifications = getModifications_(variable_modifications_names);
    if (variable_modifications.size() > 9)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Comet supports at most 9 variable modifications. " + String(variable_modifications.size()) + " provided.");
    }

    IntList binary_modifications = getIntList_("binary_modifications");
    if (!binary_modifications.empty() && binary_modifications.size() != variable_modifications.size())
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: List of binary modifications needs to have same size as variable modifications.");
    }

    int max_variable_mods_in_peptide = getIntOption_("max_variable_mods_in_peptide");
    Size var_mod_index = 0;

    // write out user specified modifications
    for (; var_mod_index < variable_modifications.size(); ++var_mod_index)
    {
      const ResidueModification* mod = variable_modifications[var_mod_index];
      double mass = mod->getDiffMonoMass();
      String residues = mod->getOrigin();

      // support for binary groups, e.g. for SILAC
      int binary_group{0};
      if (!binary_modifications.empty())
      {
        binary_group = binary_modifications[var_mod_index];
      }

      //TODO support mod-specific limit (default for now is the overall max per peptide)
      int max_current_mod_per_peptide = max_variable_mods_in_peptide;
      //TODO support term-distances?
      int term_distance = -1;
      int nc_term = 0;

      //TODO support agglomeration of Modifications to same AA. Watch out for nc_term value then.
      if (mod->getTermSpecificity() == ResidueModification::C_TERM)
      {
        if (mod->getOrigin() == 'X')
        {
          residues = "c";
        } // else stays mod.getOrigin()
        term_distance = 0;
        // Since users need to specify mods that apply to multiple residues/terms separately
        // 3 and -1 should be equal for now.
        nc_term = 3;
      }
      else if (mod->getTermSpecificity() == ResidueModification::N_TERM)
      {
        if (mod->getOrigin() == 'X')
        {
          residues = "n";
        } // else stays mod.getOrigin()
        term_distance = 0;
        // Since users need to specify mods that apply to multiple residues/terms separately
        // 2 and -1 should be equal for now.
        nc_term = 2;
      }
      else if (mod->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
      {
        if (mod->getOrigin() == 'X')
        {
          residues = "n";
        } // else stays mod.getOrigin()
        term_distance = 0;
        nc_term = 0;
      }
      else if (mod->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
      {
        if (mod->getOrigin() == 'X')
        {
          residues = "c";
        } // else stays mod.getOrigin()
        term_distance = 0;
        nc_term = 1;
      }

      //TODO support required variable mods
      bool required = false;

      os << "variable_mod0" << var_mod_index+1 << " = " 
         << mass << " " << residues << " " 
         << binary_group << " " 
         << max_current_mod_per_peptide << " " 
         << term_distance << " " 
         << nc_term << " " 
         << required << " " 
         << "0.0" // TODO: add neutral losses (from Residue or user defined?)
         << "\n";
    }

    // fill remaining modification slots (if any) in Comet with "no modification"
    for (; var_mod_index < 9; ++var_mod_index)
    {
      os << "variable_mod0" << var_mod_index+1 << " = " << "0.0 X 0 3 -1 0 0 0.0" << "\n";
    }

    os << "max_variable_mods_in_peptide = " << getIntOption_("max_variable_mods_in_peptide") << "\n";
    os << "require_variable_mod = " << (int) (getStringOption_("require_variable_mod") == "true") << "\n";

    // fragment ion defaults
    // ion trap ms/ms:  1.0005 tolerance, 0.4 offset (mono masses), theoretical_fragment_ions = 1
    // high res ms/ms:    0.02 tolerance, 0.0 offset (mono masses), theoretical_fragment_ions = 0

    String instrument = getStringOption_("instrument");
    double bin_tol = getDoubleOption_("fragment_mass_tolerance") * 2; // convert 1-sided tolerance to bin size
    double bin_offset = getDoubleOption_("fragment_bin_offset");
    if (instrument == "low_res" && (bin_tol < 0.8 || bin_offset <= 0.2))
    {
      OPENMS_LOG_ERROR << "Fragment bin size (== 2x 'fragment_mass_tolerance') or offset is quite low for low-res instruments (Comet recommends 1.005 Da bin size & 0.4 Da offset). "
                       << "Current value: fragment bin size = " << bin_tol << "(=2x" << bin_tol/2 << ") and offset = " << bin_offset << ". Use the '-force' flag to continue anyway." << std::endl;
      if (!getFlag_("force"))
      {
        return ExitCodes::ILLEGAL_PARAMETERS;
      }
      OPENMS_LOG_ERROR << "You used the '-force'!" << std::endl;
    }
    else if (instrument == "high_res" && (bin_tol > 0.1 || bin_offset > 0.1))
    {
      OPENMS_LOG_ERROR << "Fragment bin size (== 2x 'fragment_mass_tolerance') or offset is quite high for high-res instruments (Comet recommends 0.02 Da bin size & 0.0 Da offset). "
                       << "Current value: fragment bin size = " << bin_tol << "(=2x" << bin_tol / 2 << ") and offset = " << bin_offset << ". Use the '-force' flag to continue anyway." << std::endl;
      if (!getFlag_("force"))
      {
        return ExitCodes::ILLEGAL_PARAMETERS;
      }
      OPENMS_LOG_ERROR << "You used the '-force'!" << std::endl;
    }

    os << "fragment_bin_tol = " << bin_tol << "\n";               // binning to use on fragment ions
    os << "fragment_bin_offset = " << bin_offset  << "\n";              // offset position to start the binning (0.0 to 1.0)
    os << "theoretical_fragment_ions = " << (int)(instrument == "low_res") << "\n";           // 0=use flanking bins as well; 1=use M bin only
    os << "use_A_ions = " << (int)(getStringOption_("use_A_ions")=="true") << "\n";
    os << "use_B_ions = " << (int)(getStringOption_("use_B_ions")=="true") << "\n";
    os << "use_C_ions = " << (int)(getStringOption_("use_C_ions")=="true") << "\n";
    os << "use_X_ions = " << (int)(getStringOption_("use_X_ions")=="true") << "\n";
    os << "use_Y_ions = " << (int)(getStringOption_("use_Y_ions")=="true") << "\n";
    os << "use_Z_ions = " << (int)(getStringOption_("use_Z_ions")=="true") << "\n";
    os << "use_NL_ions = " << (int)(getStringOption_("use_NL_ions")=="true") << "\n";                         // 0=no, 1=yes to consider NH3/H2O neutral loss peaks

    // output
    os << "output_sqtstream = " << 0 << "\n";                    // 0=no, 1=yes  write sqt to standard output
    os << "output_sqtfile = " << 0 << "\n";                      // 0=no, 1=yes  write sqt file
    os << "output_txtfile = " << 0 << "\n";                     // 0=no, 1=yes  write tab-delimited txt file
    os << "output_pepxmlfile = " << 1 << "\n";                   // 0=no, 1=yes  write pep.xml file
    os << "export_additional_pepxml_scores = " << 1 << "\n";     // Hidden parameter of comet that adds additional comet scores to the pep.xml

    os << "output_percolatorfile = " << !getStringOption_("pin_out").empty() << "\n";              // 0=no, 1=yes  write Percolator tab-delimited input file
    os << "print_expect_score = " << 1 << "\n";                  // 0=no, 1=yes to replace Sp with expect in out & sqt
    os << "num_output_lines = " << getIntOption_("num_hits") << "\n";                    // num peptide results to show
    os << "show_fragment_ions = " << 0 << "\n";                  // 0=no, 1=yes for out files only
    os << "sample_enzyme_number = " << enzyme_number << "\n";                // Sample enzyme which is possibly different than the one applied to the search.

    // mzXML parameters
    map<string,int> override_charge;
    override_charge["keep any known"] = 0;
    override_charge["ignore known"] = 1;
    override_charge["ignore outside range"] = 2;
    override_charge["keep known search unknown"] = 3;

    int precursor_charge_min(0), precursor_charge_max(0);
    if (!parseRange_(getStringOption_("precursor_charge"), precursor_charge_min, precursor_charge_max))
    {
      OPENMS_LOG_INFO << "precursor_charge range not set. Defaulting to 0:0 (disable charge filtering)." << endl;
    }

    os << "scan_range = " << "0 0" << "\n";                        // start and scan scan range to search; 0 as 1st entry ignores parameter
    os << "precursor_charge = " << precursor_charge_min << " " << precursor_charge_max << "\n";                  // precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter
    os << "override_charge = " << override_charge[getStringOption_("override_charge")] << "\n";                     // 0=no, 1=override precursor charge states, 2=ignore precursor charges outside precursor_charge range, 3=see online
    os << "ms_level = " << getIntOption_("ms_level") << "\n";                            // MS level to analyze, valid are levels 2 (default) or 3
    os << "activation_method = " << getStringOption_("activation_method") << "\n";                 // activation method; used if activation method set; allowed ALL, CID, ECD, ETD, PQD, HCD, IRMPD

    // misc parameters
    double digest_mass_range_min(600.0), digest_mass_range_max(5000.0);
    if (!parseRange_(getStringOption_("digest_mass_range"), digest_mass_range_min, digest_mass_range_max))
    {
      OPENMS_LOG_INFO << "digest_mass_range not set. Defaulting to 600.0 5000.0." << endl;
    }

    os << "digest_mass_range = " << digest_mass_range_min << " " << digest_mass_range_max << "\n";        // MH+ peptide mass range to analyze
    os << "num_results = " << 100 << "\n";                       // number of search hits to store internally
    os << "skip_researching = " << 1 << "\n";                    // for '.out' file output only, 0=search everything again (default), 1=don't search if .out exists
    os << "max_fragment_charge = " << getIntOption_("max_fragment_charge") << "\n";                 // set maximum fragment charge state to analyze (allowed max 5)
    os << "max_precursor_charge = " << getIntOption_("max_precursor_charge") << "\n";                // set maximum precursor charge state to analyze (allowed max 9)
    os << "nucleotide_reading_frame = " << 0 << "\n";            // 0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six
    os << "clip_nterm_methionine = " << (int)(getStringOption_("clip_nterm_methionine")=="true") << "\n";              // 0=leave sequences as-is; 1=also consider sequence w/o N-term methionine
    os << "peptide_length_range = " << getIntOption_("min_peptide_length") << " " << getIntOption_("max_peptide_length") << "\n";                       // minimum and maximum peptide length to analyze (default 5 63; max length 63)
    os << "spectrum_batch_size = " << getIntOption_("spectrum_batch_size") << "\n";                 // max. // of spectra to search at a time; 0 to search the entire scan range in one loop
    os << "max_duplicate_proteins = 20\n";                       // maximum number of protein names to report for each peptide identification; -1 reports all duplicates
    os << "equal_I_and_L = 1\n";
    os << "output_suffix = " << "" << "\n";                      // add a suffix to output base names i.e. suffix "-C" generates base-C.pep.xml from base.mzXML input
    os << "mass_offsets = " << ListUtils::concatenate(getDoubleList_("mass_offsets"), " ") << "\n"; // one or more mass offsets to search (values subtracted from deconvoluted precursor mass)
    os << "precursor_NL_ions =\n"; //  one or more precursor neutral loss masses, will be added to xcorr analysis 

    // spectral processing
    map<string,int> remove_precursor_peak;
    remove_precursor_peak["no"] = 0;
    remove_precursor_peak["yes"] = 1;
    remove_precursor_peak["charge_reduced"] = 2;
    remove_precursor_peak["phosphate_loss"] = 3;

    double clear_mz_range_min(0.0), clear_mz_range_max(0.0);
    if (!parseRange_(getStringOption_("clear_mz_range"), clear_mz_range_min, clear_mz_range_max))
    {
      OPENMS_LOG_INFO << "clear_mz_range not set. Defaulting to 0:0 (disable m/z filter)." << endl;
    }

    os << "minimum_peaks = " << getIntOption_("minimum_peaks") << "\n";                      // required minimum number of peaks in spectrum to search (default 10)
    os << "minimum_intensity = " << getDoubleOption_("minimum_intensity") << "\n";                   // minimum intensity value to read in
    os << "remove_precursor_peak = " << remove_precursor_peak[getStringOption_("remove_precursor_peak")] << "\n";               // 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD), 3=phosphate neutral loss peaks
    os << "remove_precursor_tolerance = " << getDoubleOption_("remove_precursor_tolerance") << "\n";        // +- Da tolerance for precursor removal
    os << "clear_mz_range = " << clear_mz_range_min << " " << clear_mz_range_max << "\n";                // for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range


    // write fixed modifications - if not specified residue parameter is zero
    // Aminoacid:
    //      add_AA.OneletterCode_AA.ThreeLetterCode = xxx
    // Terminus:
    //      add_N/Cterm_peptide = xxx       protein not available yet
    vector<String> fixed_modifications_names = getStringList_("fixed_modifications");
    const vector<const ResidueModification*> fixed_modifications = getModifications_(fixed_modifications_names);

    // merge duplicates, targeting the same AA
    std::map<String, double> mods;
    // Comet sets Carbamidometyl (C) as modification as default even if not specified.
    // Therefore there is the need to set it to 0, unless its set as flag (see loop below)
    mods["add_C_cysteine"] = 0;

    for (const auto& fm : fixed_modifications)
    {
      // check modification (amino acid or terminal)
      String AA = fm->getOrigin(); // X (constructor) or amino acid (e.g. K)
      String term_specificity = fm->getTermSpecificityName(); // N-term, C-term, none
      if ((AA != "X") && (term_specificity == "none"))
      {
        const Residue* r = ResidueDB::getInstance()->getResidue(AA);
        String name = r->getName();
        mods["add_" + r->getOneLetterCode() + "_" + name.toLower()] += fm->getDiffMonoMass();
      }
      else if (term_specificity == "N-term" || term_specificity == "C-term")
      {
        mods["add_" + term_specificity.erase(1,1) + "_peptide"] += fm->getDiffMonoMass();
      }
      else if (term_specificity == "Protein N-term" || term_specificity == "Protein C-term")
      {
        term_specificity.erase(0,8); // remove "Protein "
        mods["add_" + term_specificity.erase(1,1) + "_protein"] += fm->getDiffMonoMass();
      }
    }
    for (const auto& mod : mods)
    {
      os << mod.first << " = " << mod.second << "\n";
    }

    //TODO register cut_before and cut_after in Enzymes.xml plus datastructures to add all our Enzymes with our names instead.
    // COMET_ENZYME_INFO _must_ be at the end of this parameters file
    os << "[COMET_ENZYME_INFO]" << "\n";
    os << "0.  No_enzyme              0      -           -" << "\n";
    os << "1.  Trypsin                1      KR          P" << "\n";
    os << "2.  Trypsin/P              1      KR          -" << "\n";
    os << "3.  Lys_C                  1      K           P" << "\n";
    os << "4.  Lys_N                  0      K           -" << "\n";
    os << "5.  Arg_C                  1      R           P" << "\n";
    os << "6.  Asp_N                  0      D           -" << "\n";
    os << "7.  CNBr                   1      M           -" << "\n";
    os << "8.  Glu_C                  1      DE          P" << "\n";
    os << "9.  PepsinA                1      FL          P" << "\n";
    os << "10. Chymotrypsin           1      FWYL        P" << "\n";
    os << "11. No_cut                 1      @           @" << "\n";
    os << "12. Arg-C/P                1.     R           _" << "\n";

    return ExitCodes::EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    // do this early, to see if comet is installed
    String comet_executable = getStringOption_("comet_executable");
    File::TempDir tmp_dir(debug_level_ >= 2);

    writeDebug_("Comet is writing the default parameter file...", 1);
    
    TOPPBase::ExitCodes exit_code = runExternalProcess_(comet_executable.toQString(), QStringList() << "-p", tmp_dir.getPath().toQString());
    if (exit_code != EXECUTION_OK)
    {
      return EXTERNAL_PROGRAM_ERROR;
    }
    // the first line of 'comet.params.new' contains a string like: "# comet_version 2017.01 rev. 1"
    String comet_version; 
    {
      std::ifstream ifs(tmp_dir.getPath() + "/comet.params.new");
      getline(ifs, comet_version);
    }
    writeDebug_("Comet Version extracted is: '" + comet_version + "\n", 2);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------


    int ms_level = getIntOption_("ms_level");
    String inputfile_name = getRawfileName(ms_level);
    String out = getStringOption_("out");
    String db_name = getDBFilename();

    // tmp_dir
    String tmp_pepxml = tmp_dir.getPath() + "result.pep.xml";
    String tmp_pin = tmp_dir.getPath() + "result.pin";
    String default_params = getStringOption_("default_params_file");
    String tmp_file;

    // default params given or to be written
    if (default_params.empty())
    {
        tmp_file = tmp_dir.getPath() + "param.txt";
        ofstream os(tmp_file.c_str());
        auto ret = createParamFile_(os, comet_version);
        os.close();
        if (ret != EXECUTION_OK)
        {
          return ret;
        }
    }
    else
    {
        tmp_file = default_params;
    }

    // check for mzML index (comet requires one)
    MSExperiment exp;
    MzMLFile mzml_file{};
    String input_file_with_index = inputfile_name;
    if (!mzml_file.hasIndex(inputfile_name))
    {
      OPENMS_LOG_WARN << "The mzML file provided to CometAdapter is not indexed, but comet requires one. "
                      << "We will add an index by writing a temporary file. If you run this analysis more often, consider indexing your mzML in advance!" << std::endl;
      // Low memory conversion
      // write mzML with index again
      auto tmp_file = File::getTemporaryFile() + ".mzML";
      PlainMSDataWritingConsumer consumer(tmp_file);
      consumer.getOptions().addMSLevel(ms_level); // only load msLevel 2
      bool skip_full_count = true;
      mzml_file.transform(inputfile_name, &consumer, skip_full_count);
      input_file_with_index = tmp_file;
    }

    mzml_file.getOptions().setMetadataOnly(true);
    mzml_file.load(inputfile_name, exp); // always load metadata for raw file name

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    String paramP = "-P" + tmp_file;
    String paramN = "-N" + FileHandler::stripExtension(FileHandler::stripExtension(tmp_pepxml));
    QStringList arguments;
    arguments << paramP.toQString() << paramN.toQString() << input_file_with_index.toQString();

    //-------------------------------------------------------------
    // run comet
    //-------------------------------------------------------------
    // Comet execution with the executable and the arguments StringList
    exit_code = runExternalProcess_(comet_executable.toQString(), arguments);
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }
    //-------------------------------------------------------------
    // writing IdXML output
    //-------------------------------------------------------------

    // read the pep.xml put of Comet and write it to idXML

    vector<String> fixed_modifications_names = getStringList_("fixed_modifications");
    vector<String> variable_modifications_names = getStringList_("variable_modifications");

    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;

    writeDebug_("load PepXMLFile", 1);
    PepXMLFile pepfile{};
    pepfile.setPreferredFixedModifications(getModifications_(fixed_modifications_names));
    pepfile.setPreferredVariableModifications(getModifications_(variable_modifications_names));
    pepfile.load(tmp_pepxml, protein_identifications, peptide_identifications);
    writeDebug_("write idXMLFile", 1);
    writeDebug_(out, 1);

    //Whatever the pepXML says, overwrite origin as the input mzML
    protein_identifications[0].setPrimaryMSRunPath({inputfile_name}, exp);
    // seems like version is not correctly parsed from pepXML. Overwrite it here.
    protein_identifications[0].setSearchEngineVersion(comet_version);
    // TODO let this be parsed by the pepXML parser if this info is present there.
    protein_identifications[0].getSearchParameters().enzyme_term_specificity =
    static_cast<EnzymaticDigestion::Specificity>(num_enzyme_termini[getStringOption_("num_enzyme_termini")]);
    protein_identifications[0].getSearchParameters().charges = getStringOption_("precursor_charge");
    protein_identifications[0].getSearchParameters().db = getStringOption_("database");

    // write all (!) parameters as metavalues to the search parameters
    if (!protein_identifications.empty())
    {
      DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_identifications[0].getSearchParameters(), this->getToolPrefix());
    }

    // if "reindex" parameter is set to true will perform reindexing
    if (auto ret = reindex_(protein_identifications, peptide_identifications); ret != EXECUTION_OK) return ret;

    FileHandler().storeIdentifications(out, protein_identifications, peptide_identifications, {FileTypes::IDXML});

    //-------------------------------------------------------------
    // create (move) optional pin output
    //-------------------------------------------------------------

    String pin_out = getStringOption_("pin_out");
    if (!pin_out.empty())
    { // move the temporary file to the actual destination:
      if (!File::rename(tmp_pin, pin_out))
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPCometAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
