// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>

#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <regex>

#include <QStringList>
#include <chrono>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

#include <boost/math/distributions/normal.hpp>

using namespace OpenMS;
using namespace std;
using boost::math::normal;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SageAdapter SageAdapter

@brief Identifies peptides in MS/MS spectra via sage.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; SageAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

@em Sage must be installed before this wrapper can be used.

Only the closed-search identification mode of Sage is supported by this adapter.
Currently, also neither "wide window" (= open or DIA) mode, nor "chimeric" mode is supported,
because of limitations in OpenMS' data structures and file formats.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SageAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SageAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


#define CHRONOSET



class TOPPSageAdapter :
  public SearchEngineBase
{
public: 
  TOPPSageAdapter() :
    SearchEngineBase("SageAdapter", "Annotates MS/MS spectra using Sage.", true,
             {
                 {"Michael Lazear",
                 "Sage: An Open-Source Tool for Fast Proteomics Searching and Quantification at Scale",
                 "J. Proteome Res. 2023, 22, 11, 3652–3659",
                 "https://doi.org/10.1021/acs.jproteome.3c00486"}
             })
  {
  }

  // Note: Modification analysis functionality has been moved to OpenSearchModificationAnalysis class


protected:
  // create a template-based configuration file for sage
  // variable values correspond to sage parameter that can be configured via TOPP tool parameter.
  // values will be pasted into the config_template at the corresponding tag. E.g. bucket_size at tag ##bucket_size##
  static constexpr size_t bucket_size = 8192;
  static constexpr size_t min_len = 5; 
  static constexpr size_t max_len = 50; 
  static constexpr size_t missed_cleavages = 2;
  static constexpr double fragment_min_mz = 200.0;
  static constexpr double fragment_max_mz = 2000.0;
  static constexpr double peptide_min_mass = 500.0;
  static constexpr double peptide_max_mass = 5000.0;
  static constexpr size_t min_ion_index = 2;
  static constexpr size_t max_variable_mods = 2;
  const std::string precursor_tol_unit = "ppm";
  static constexpr double precursor_tol_left = -6.0;
  static constexpr double precursor_tol_right = 6.0;
  const std::string fragment_tol_unit = "ppm";
  static constexpr double fragment_tol_left = -10.0;
  static constexpr double fragment_tol_right = 10.0;
  const std::string isotope_errors = "-1, 3";
  const std::string charges_if_not_annotated = "2, 5";
  static constexpr size_t min_matched_peaks = 6;
  static constexpr size_t report_psms = 1;
  static constexpr size_t min_peaks = 15;
  static constexpr size_t max_peaks = 150;

  std::string config_template = R"(
{
  "database": {
    "bucket_size": ##bucket_size##,
    "enzyme": {
      "missed_cleavages": ##missed_cleavages##,
      "min_len": ##min_len##,
      "max_len": ##max_len##,
      ##enzyme_details##
    },
    "fragment_min_mz": ##fragment_min_mz##,
    "fragment_max_mz": ##fragment_max_mz##,
    "peptide_min_mass": ##peptide_min_mass##,
    "peptide_max_mass": ##peptide_max_mass##,
    "ion_kinds": ["b", "y"],
    "min_ion_index": ##min_ion_index##,
    "static_mods": {
      ##static_mods##
    },
    "variable_mods": {
      ##variable_mods##
    },
    "max_variable_mods": ##max_variable_mods##,
    "generate_decoys": false,
    "decoy_tag": "##decoy_prefix##"
  },
  "precursor_tol": {
    "##precursor_tol_unit##": [
      ##precursor_tol_left##,
      ##precursor_tol_right##
    ]
  },
  "fragment_tol": {
    "##fragment_tol_unit##": [
    ##fragment_tol_left##,
    ##fragment_tol_right##
    ]
  },
  "precursor_charge": [
    ##charges_if_not_annotated##
  ],
  "isotope_errors": [
    ##isotope_errors##
  ],
  "deisotope": ##deisotope##,
  "chimera": ##chimera##,
  "predict_rt": ##predict_rt##,
  "min_peaks": ##min_peaks##,
  "max_peaks": ##max_peaks##,
  "min_matched_peaks": ##min_matched_peaks##,
  "report_psms": ##report_psms##, 
  "wide_window": ##wide_window##
}
)";

  // formats a single mod entry as sage json entry
  String getModDetails(const ResidueModification* mod, const Residue* res)
  {
    String origin;
    if (mod->getTermSpecificity() == ResidueModification::N_TERM)
    { 
      origin += "^";
    }
    else if (mod->getTermSpecificity() == ResidueModification::C_TERM)
    {
      origin += "$";
    }
    else if (mod->getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
    {
      origin += "[";
    }
    else if (mod->getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
    {
      origin += "]";
    }
   if (res != nullptr && res->getOneLetterCode() != "X") // omit letter for "any AA"
   {
     origin += res->getOneLetterCode();
   }

    return String("\"") + origin + "\": " + String(mod->getDiffMonoMass());
  }

  // formats all mod entries into a single multi-line json string
  String getModDetailsString(const OpenMS::ModifiedPeptideGenerator::MapToResidueType& mod_map)
  {
    String mod_details;
    for (auto it = mod_map.val.begin(); it != mod_map.val.end(); ++it)
    {
      const auto& mod = it->first;
      const auto& res = it->second;
      mod_details += getModDetails(mod, res);     
      if (std::next(it) != mod_map.val.end())
      {
        mod_details += ",\n";
      }
    }
    return mod_details;
  }

  // impute values into config_template
  // TODO just iterate over all options??
  String imputeConfigIntoTemplate()
  {
    String config_file = config_template;
    config_file.substitute("##bucket_size##", String(getIntOption_("bucket_size")));
    config_file.substitute("##min_len##", String(getIntOption_("min_len")));
    config_file.substitute("##max_len##", String(getIntOption_("max_len")));
    config_file.substitute("##missed_cleavages##", String(getIntOption_("missed_cleavages")));
    config_file.substitute("##fragment_min_mz##", String(getDoubleOption_("fragment_min_mz")));
    config_file.substitute("##fragment_max_mz##", String(getDoubleOption_("fragment_max_mz")));
    config_file.substitute("##peptide_min_mass##", String(getDoubleOption_("peptide_min_mass")));
    config_file.substitute("##peptide_max_mass##", String(getDoubleOption_("peptide_max_mass")));
    config_file.substitute("##min_ion_index##", String(getIntOption_("min_ion_index")));
    config_file.substitute("##max_variable_mods##", String(getIntOption_("max_variable_mods")));
    config_file.substitute("##precursor_tol_unit##", getStringOption_("precursor_tol_unit") == "Da" ? "da" : "ppm"); // sage might expect lower-case "da"
    config_file.substitute("##precursor_tol_left##", String(getDoubleOption_("precursor_tol_left")));
    config_file.substitute("##precursor_tol_right##", String(getDoubleOption_("precursor_tol_right")));
    config_file.substitute("##fragment_tol_unit##", getStringOption_("fragment_tol_unit") == "Da" ? "da" : "ppm"); // sage might expect lower-case "da"
    config_file.substitute("##fragment_tol_left##", String(getDoubleOption_("fragment_tol_left")));
    config_file.substitute("##fragment_tol_right##", String(getDoubleOption_("fragment_tol_right")));
    config_file.substitute("##isotope_errors##", getStringOption_("isotope_error_range"));
    config_file.substitute("##charges_if_not_annotated##", getStringOption_("charges"));
    config_file.substitute("##min_matched_peaks##", String(getIntOption_("min_matched_peaks")));
    config_file.substitute("##min_peaks##", String(getIntOption_("min_peaks")));
    config_file.substitute("##max_peaks##", String(getIntOption_("max_peaks")));
    config_file.substitute("##report_psms##", String(getIntOption_("report_psms")));
    config_file.substitute("##deisotope##", getStringOption_("deisotope")); 
    config_file.substitute("##chimera##", getStringOption_("chimera")); 
    config_file.substitute("##predict_rt##", getStringOption_("predict_rt")); 
    config_file.substitute("##decoy_prefix##", getStringOption_("decoy_prefix")); 
    config_file.substitute("##wide_window##", getStringOption_("wide_window")); 

    
    //Look at decoy handling 

    String enzyme = getStringOption_("enzyme");
    String enzyme_details;
    if (enzyme == "Trypsin")
    {
      enzyme_details = 
   R"("cleave_at": "KR",
      "restrict": "P",
      "c_terminal": true)";
    }
    else if (enzyme == "Trypsin/P")
    {
      enzyme_details = 
   R"("cleave_at": "KR",
      "restrict": null,
      "c_terminal": true)";
    }
    else if (enzyme == "Chymotrypsin")
    {
      enzyme_details = 
   R"("cleave_at": "FWYL",
      "restrict": "P",
      "c_terminal": true)";
    }
    else if (enzyme == "Chymotrypsin/P")
    {
      enzyme_details = 
   R"("cleave_at": "FWYL",
      "restrict": null,
      "c_terminal": true)";
    }
    else if (enzyme == "Arg-C")
    {
      enzyme_details = 
   R"("cleave_at": "R",
      "restrict": "P",
      "c_terminal": true)";
    }
    else if (enzyme == "Arg-C/P")
    {
      enzyme_details = 
   R"("cleave_at": "R",
      "restrict": null,
      "c_terminal": true)";
    }
    else if (enzyme == "Lys-C")
    {
      enzyme_details = 
   R"("cleave_at": "K",
      "restrict": "P",
      "c_terminal": true)";
    }
    else if (enzyme == "Lys-C/P")
    {
      enzyme_details = 
   R"("cleave_at": "K",
      "restrict": null,
      "c_terminal": true)";
    }    
    else if (enzyme == "Lys-N")
    {
      enzyme_details = 
   R"("cleave_at": "K",
      "restrict": null,
      "c_terminal": false)";
    }
    else if (enzyme == "no cleavage")
    {
      enzyme_details = 
   R"("cleave_at": "$")";
    }    
    else if (enzyme == "unspecific cleavage")
    {
      enzyme_details = 
   R"("cleave_at": "")";
    }
    else if (enzyme == "glutamyl endopeptidase")
    {
      enzyme_details =
   R"("cleave_at": "E",
      "restrict": "E",
      "c_terminal":true)";
    }
    else if (enzyme == "leukocyte elastase")
    {
      enzyme_details =
   R"("cleave_at": "ALIV",
      "restrict": null,
      "c_terminal":true)";
    }

    config_file.substitute("##enzyme_details##", enzyme_details);

    
    auto fixed_mods = getStringList_("fixed_modifications");
    set<String> fixed_unique(fixed_mods.begin(), fixed_mods.end());
    fixed_mods.assign(fixed_unique.begin(), fixed_unique.end());   
    ModifiedPeptideGenerator::MapToResidueType fixed_mod_map = ModifiedPeptideGenerator::getModifications(fixed_mods); // std::unordered_map<const ResidueModification*, const Residue*> val;
    String static_mods_details = getModDetailsString(fixed_mod_map);

    auto variable_mods = getStringList_("variable_modifications");
    set<String> variable_unique(variable_mods.begin(), variable_mods.end());
    variable_mods.assign(variable_unique.begin(), variable_unique.end());
    ModifiedPeptideGenerator::MapToResidueType variable_mod_map = ModifiedPeptideGenerator::getModifications(variable_mods);
    String variable_mods_details = getModDetailsString(variable_mod_map);

    //Treat variables as list for sage v0.15 and beyond 
    StringList static_mods_details_list; 
    StringList variable_mods_details_list; 

    String static_mods_details_split = static_mods_details; 
    String variable_mods_details_split = variable_mods_details; 
    static_mods_details_split.split(",", static_mods_details_list); 
    variable_mods_details_split.split(",", variable_mods_details_list); 

    String temp_String_var; 
    for (auto& x : variable_mods_details_list)
    {
      StringList temp_split; 
      x.split(":", temp_split); 
      
      temp_split.insert(temp_split.begin()+1, ":["); 
      temp_split.insert(temp_split.end(), "]"); 
      String temp_split_Str = ""; 

      for (auto& y : temp_split)
      {
        temp_split_Str = temp_split_Str + y; 
      } 
      temp_String_var = temp_String_var + "," + temp_split_Str ; 
    } 
    String temp_String_var_Fin = temp_String_var.substr(1, temp_String_var.size()-1); 
    config_file.substitute("##static_mods##", static_mods_details);
    config_file.substitute("##variable_mods##", temp_String_var_Fin);

    return config_file;
  }

  std::tuple<std::string, std::string, std::string> getVersionNumber_(const std::string& multi_line_input)
  {
      std::regex version_regex("Version ([0-9]+)\\.([0-9]+)\\.([0-9]+)");

      std::sregex_iterator it(multi_line_input.begin(), multi_line_input.end(), version_regex);
      std::smatch match = *it;
      std::cout << "Found Sage version string: " << match.str() << std::endl;      
          
      return make_tuple(it->str(1), it->str(2), it->str(3)); // major, minor, patch
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", { "mzML" } );

    registerOutputFile_("out", "<file>", "", "Single output file containing all search results.", true, false);
    setValidFormats_("out", { "idXML" } );

    registerInputFile_("database", "<file>", "", "FASTA file", true, false, {"skipexists"});
    setValidFormats_("database", { "FASTA" } );

    registerInputFile_("sage_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      #ifdef OPENMS_WINDOWSPLATFORM
        "sage.exe",
      #else
        "sage",
      #endif
      "The Sage executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", true, false, {"is_executable"});

    registerStringOption_("decoy_prefix", "<prefix>", "DECOY_", "Prefix on protein accession used to distinguish decoy from target proteins. Decoy proteins in the FASTA file should have this prefix in their accession. NOTE: Decoy suffix is currently not supported by Sage.", false, false);
    registerIntOption_("batch_size", "<int>", 0, "Number of files to load and search in parallel. Setting this to 0 (default) uses an automatic value (typically number of CPUs/2). Default: 0", false, false);
    
    registerDoubleOption_("precursor_tol_left", "<double>", -6.0, "Start (left side) of the precursor tolerance window w.r.t. precursor location. This value is relative to the experimental precursor mass and used to define the lower bound of the search window. Must be negative (e.g., -6 ppm means 6 ppm below the observed mass).", false, false);
    registerDoubleOption_("precursor_tol_right", "<double>", 6.0, "End (right side) of the precursor tolerance window w.r.t. precursor location. This value is added to the experimental precursor mass to define the upper bound of the search window. Must be positive (e.g., 6 ppm means 6 ppm above the observed mass).", false, false);
    registerStringOption_("precursor_tol_unit", "<unit>", "ppm", "Unit of precursor tolerance (ppm or Da)", false, false);
    setValidStrings_("precursor_tol_unit", ListUtils::create<String>("ppm,Da"));

    registerDoubleOption_("fragment_tol_left", "<double>", -20.0, "Start (left side) of the fragment tolerance window w.r.t. fragment location. This value reduces the experimental fragment mass to define the lower bound of the search window. Must be negative (e.g., -20 ppm means 20 ppm below the observed mass).", false, false);
    registerDoubleOption_("fragment_tol_right", "<double>", 20.0, "End (right side) of the fragment tolerance window w.r.t. fragment location. This value is added to the experimental fragment mass to define the upper bound of the search window. Must be positive (e.g., 20 ppm means 20 ppm above the observed mass).", false, false);
    registerStringOption_("fragment_tol_unit", "<unit>", "ppm", "Unit of fragment tolerance (ppm or Da)", false, false);
    setValidStrings_("fragment_tol_unit", ListUtils::create<String>("ppm,Da"));

    // add advanced options
    registerIntOption_("min_matched_peaks", "<int>", min_matched_peaks, "Minimum number of b+y ions required to match for PSM to be reported. Default: 6", false, true);
    registerIntOption_("min_peaks", "<int>", min_peaks, "Minimum number of peaks required for a spectrum to be considered. Spectra with fewer peaks will be ignored. Default: 15", false, true);
    registerIntOption_("max_peaks", "<int>", max_peaks, "Take the top N most intense MS2 peaks to search. Default: 150", false, true);
    registerIntOption_("report_psms", "<int>", report_psms, "Number of peptide-spectrum matches (PSMs) to report for each spectrum. The top N scoring PSMs will be reported. Values higher than 1 can be useful for chimeric spectra but may affect downstream statistical analysis. Default: 1", false, true);
    registerIntOption_("bucket_size", "<int>", bucket_size, "How many fragments are in each internal mass bucket. Default: 8192 (optimal for high-resolution data). Try increasing it to 32768 or 65536 for low-resolution data. See also: fragment_tol_*", false, true);
    registerIntOption_("min_len", "<int>", min_len, "Minimum peptide length (in amino acids). Default: 5", false, true);
    registerIntOption_("max_len", "<int>", max_len, "Maximum peptide length (in amino acids). Default: 50", false, true);
    registerIntOption_("missed_cleavages", "<int>", missed_cleavages, "Maximum number of missed enzymatic cleavages to allow in peptide generation. Default: 2", false, true);
    registerDoubleOption_("fragment_min_mz", "<double>", fragment_min_mz, "Minimum fragment m/z to consider. Fragment ions below this m/z will be ignored. Default: 200.0", false, true);
    registerDoubleOption_("fragment_max_mz", "<double>", fragment_max_mz, "Maximum fragment m/z to consider. Fragment ions above this m/z will be ignored. Default: 2000.0", false, true);
    registerDoubleOption_("peptide_min_mass", "<double>", peptide_min_mass, "Minimum monoisotopic peptide mass to consider for in silico digestion. Peptides below this mass will be excluded from the search database. Default: 500.0", false, true);
    registerDoubleOption_("peptide_max_mass", "<double>", peptide_max_mass, "Maximum monoisotopic peptide mass to consider for in silico digestion. Peptides above this mass will be excluded from the search database. Default: 5000.0", false, true);
    registerIntOption_("min_ion_index", "<int>", min_ion_index, "Minimum ion index to consider for preliminary scoring. This parameter controls which fragment ions are used in preliminary scoring. Default: 2 (skips b1/b2/y1/y2 ions, which are often missing or unreliable). Setting this to 1 would only skip b1/y1 ions. Does not affect the final scoring of PSMs.", false, true);
    registerIntOption_("max_variable_mods", "<int>", max_variable_mods, "Maximum number of variable modifications allowed per peptide. Default: 2", false, true);  
    registerStringOption_("isotope_error_range", "<start,end>", isotope_errors, "Range of C13 isotope errors to consider for precursor matching, specified as 'start,end' (e.g., '-1,3'). For a range of '-1,3', Sage will consider all isotope errors from -1 to +3 (i.e., -1, 0, 1, 2, 3). This is useful when the monoisotopic peak may not be selected. Can include negative values. Default: '-1,3'. Note: Searching with isotope errors is slower than using a wider precursor tolerance.", false, true);
    registerStringOption_("charges", "<start,end>", charges_if_not_annotated, "Range of precursor charge states to consider if not annotated in the file, specified as 'start,end' (e.g., '2,5'). For a range of '2,5', Sage will consider charge states 2, 3, 4, and 5. This is only used when charge state information is missing from the input file. Default: '2,5'"
      , false, true);
    

    //Search Enzyme
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false, false);
    setValidStrings_("enzyme", all_enzymes);

    //Modifications
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>("Carbamidomethyl (C)", ','), "Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>("Oxidation (M)", ','), "Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    //FDR and misc 

    registerDoubleOption_("q_value_threshold", "<double>", 1, "The FDR (False Discovery Rate) threshold for filtering peptides. PSMs with q-values above this threshold will be excluded. Default: 1 (no filtering)", false, false); 
    registerStringOption_("annotate_matches", "<bool>", "true", "Whether fragment ion matches should be annotated in the output. This provides additional information about which theoretical ions matched experimental peaks. Default: true", false, false); 
    registerStringOption_("deisotope", "<bool>", "false", "Perform deisotoping and charge state deconvolution on MS2 spectra. Recommended for high-resolution MS2 data. May interfere with TMT-MS2 quantification. Default: false", false, false ); 
    registerStringOption_("chimera", "<bool>", "false", "Enable chimeric spectra search mode. When enabled, multiple peptide identifications can be reported for each MS2 scan, useful for co-fragmenting peptides. Default: false", false, false  ); 
    registerStringOption_("predict_rt",  "<bool>", "false", "Use retention time prediction model as a feature for machine learning scoring. Note: This is incompatible with label-free quantification (LFQ). Default: false", false, false ); 
    registerStringOption_("wide_window", "<bool>", "false", "Enable wide-window/DIA search mode. When enabled, the precursor_tol parameter is ignored and a dynamic precursor tolerance is used. Default: false", false, false);
    registerStringOption_("smoothing", "<bool>", "true", "Whether to smooth the PTM (post-translational modification) mass histogram and pick local maxima. If false, uses raw histogram data. Default: true", false, false);  
    registerIntOption_("threads", "<int>", 0, "Number of threads to use. 0 (default) = auto-detect and use all available CPUs via Rayon", false, false); 

    // register peptide indexing parameter (with defaults for this search engine)
    registerPeptideIndexingParameter_(PeptideIndexing().getParameters());
  }


  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    // Validate tolerance parameters
    double precursor_tol_left = getDoubleOption_("precursor_tol_left");
    double precursor_tol_right = getDoubleOption_("precursor_tol_right");
    double fragment_tol_left = getDoubleOption_("fragment_tol_left");
    double fragment_tol_right = getDoubleOption_("fragment_tol_right");

    // Warn if tolerance parameters seem incorrect
    if (precursor_tol_left > 0)
    {
      OPENMS_LOG_WARN << "WARNING: precursor_tol_left is positive (" << precursor_tol_left << "). "
                      << "This parameter is used to reduce the experimental mass, so it should typically be negative. "
                      << "A positive value will likely produce an incorrect search window." << std::endl;
    }
    if (precursor_tol_right < 0)
    {
      OPENMS_LOG_WARN << "WARNING: precursor_tol_right is negative (" << precursor_tol_right << "). "
                      << "This parameter is ADDED to the experimental mass, so it should typically be positive. "
                      << "A negative value will likely produce an incorrect search window." << std::endl;
    }
    if (fragment_tol_left > 0)
    {
      OPENMS_LOG_WARN << "WARNING: fragment_tol_left is positive (" << fragment_tol_left << "). "
                      << "This parameter is used to reduce the experimental mass, so it should typically be negative. "
                      << "A positive value will likely produce an incorrect search window." << std::endl;
    }
    if (fragment_tol_right < 0)
    {
      OPENMS_LOG_WARN << "WARNING: fragment_tol_right is negative (" << fragment_tol_right << "). "
                      << "This parameter is ADDED to the experimental mass, so it should typically be positive. "
                      << "A negative value will likely produce an incorrect search window." << std::endl;
    }

    // do this early, to see if Sage is installed
    String sage_executable = getStringOption_("sage_executable");
    std::cout << sage_executable << " sage executable" << std::endl; 
    String proc_stdout, proc_stderr;
    TOPPBase::ExitCodes exit_code = runExternalProcess_(sage_executable.toQString(), QStringList() << "--help", proc_stdout, proc_stderr, "");
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    auto major_minor_patch = getVersionNumber_(proc_stdout);
    String sage_version = std::get<0>(major_minor_patch) + "." + std::get<1>(major_minor_patch) + "." + std::get<2>(major_minor_patch);
    
    //-------------------------------------------------------------
    // run sage
    //-------------------------------------------------------------
    StringList input_files = getStringList_("in");
    String output_file = getStringOption_("out");
    String output_folder = File::path(output_file);
    String fasta_file = getStringOption_("database");
    int batch = getIntOption_("batch_size");
    int threads = getIntOption_("threads");
    String decoy_prefix = getStringOption_("decoy_prefix");

    // create config
    String config = imputeConfigIntoTemplate();

    // store config in config_file
    OPENMS_LOG_INFO << "Creating temp file name..." << std::endl;
    String config_file = File::getTempDirectory() + "/" + File::getUniqueName() + ".json";
    OPENMS_LOG_INFO << "Creating Sage config file..." << config_file << std::endl;
    ofstream config_stream(config_file.c_str());
    config_stream << config;
    config_stream.close();

    // keep config file if debug mode is set
    if (getIntOption_("debug") > 1)
    {
      String debug_config_file = output_folder + "/" + File::getUniqueName() + ".json";
      ofstream debug_config_stream(debug_config_file.c_str());
      debug_config_stream << config;
      debug_config_stream.close();     
    }

    String annotation_check;    

    QStringList arguments;

  if ( (getStringOption_("annotate_matches").compare("true")) == 0)
  {
    arguments << config_file.toQString() 
              << "-f" << fasta_file.toQString() 
              << "-o" << output_folder.toQString() 
              << "--annotate-matches"
              << "--write-pin"; 
  }
  else
  {
    arguments << config_file.toQString() 
              << "-f" << fasta_file.toQString() 
              << "-o" << output_folder.toQString() 
              << "--write-pin"; 
  }

    if (batch >= 1) arguments << "--batch-size" << String(batch).toQString();
    
    for (auto s : input_files) arguments << s.toQString();

    OPENMS_LOG_INFO << "Sage command line: " << sage_executable << " " << arguments.join(' ').toStdString() << std::endl;
    
    //std::chrono lines for testing/writing purposes only! 

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    // Set RAYON_NUM_THREADS environment variable to control Sage's thread usage
    // Only set if threads > 0; if threads == 0, let Rayon auto-detect (use all CPUs)
    std::map<QString, QString> sage_env;
    if (threads > 0)
    {
      sage_env["RAYON_NUM_THREADS"] = String(threads).toQString();
    }
    
    // Sage execution with the executable and the arguments StringList
    exit_code = runExternalProcess_(sage_executable.toQString(), arguments, "", sage_env);
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    #ifdef CHRONOSET
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    #endif

    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    //-------------------------------------------------------------
    // writing IdXML output
    //-------------------------------------------------------------

    // read the sage output
    OPENMS_LOG_INFO << "Reading sage output..." << std::endl;
    StringList filenames;
    StringList extra_scores = {"ln(-poisson)", "ln(delta_best)", "ln(delta_next)", 
      "ln(matched_intensity_pct)", "longest_b", "longest_y", 
      "longest_y_pct", "matched_peaks", "scored_candidates"}; 
    double FDR_threshhold = getDoubleOption_("q_value_threshold"); 

    PeptideIdentificationList peptide_identifications = PercolatorInfile::load(
      output_folder + "/results.sage.pin",
      true,
      "ln(hyperscore)", //TODO can we get sage's "sage_discriminant_score" out of the pin? Probably not. Suboptimal!
      extra_scores,
      filenames,
      decoy_prefix, 
      FDR_threshhold, 
      true);

    for (auto& id : peptide_identifications)
    {
      auto& hits = id.getHits();
      for (auto& h : hits)
      {
        for (const auto& meta : extra_scores)
        {
          if (h.metaValueExists(meta))
          {
            h.setMetaValue("SAGE:" + meta, h.getMetaValue(meta));
            h.removeMetaValue(meta);    
          }
        }
      }
    }
    
    String smoothing_string = getStringOption_("smoothing"); 
    bool smoothing = !(smoothing_string.compare("true")); 

    // Use shared modification analysis functionality
    OpenSearchModificationAnalysis mod_analyzer;
    auto modification_summaries = mod_analyzer.analyzeModifications(
      peptide_identifications,
      0.01, // precursor_mass_tolerance
      false, // precursor_mass_tolerance_unit_ppm (0.01 Da)
      smoothing,
      output_file
    );
    // remove hits without charge state assigned or charge outside of default range (fix for downstream bugs). TODO: remove if all charges annotated in sage
    IDFilter::filterPeptidesByCharge(peptide_identifications, 2, numeric_limits<int>::max());
    
    if (filenames.empty()) filenames = getStringList_("in");

    // TODO: allow optional split and create multiple idXMLs one per input file
    vector<ProteinIdentification> protein_identifications(1, ProteinIdentification());

    writeDebug_("write idXMLFile", 1);    
    
    protein_identifications[0].setPrimaryMSRunPath(filenames);  
    protein_identifications[0].setDateTime(DateTime::now());
    protein_identifications[0].setSearchEngine("Sage");
    protein_identifications[0].setSearchEngineVersion(sage_version);

    DateTime now = DateTime::now();
    String identifier("Sage_" + now.get());
    protein_identifications[0].setIdentifier(identifier);
    for (auto & pid : peptide_identifications) 
    { 
      pid.setIdentifier(identifier);
      pid.setScoreType("ln(hyperscore)");
      pid.setHigherScoreBetter(true);
    }

    auto& search_parameters = protein_identifications[0].getSearchParameters();
    // protein_identifications[0].getSearchParameters().enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(num_enzyme_termini[getStringOption_("num_enzyme_termini")]);
    protein_identifications[0].getSearchParameters().db = getStringOption_("database");
    
    // add extra scores for percolator rescoring
    vector<String> percolator_features = { "score" };
    for (auto s : extra_scores) percolator_features.push_back("SAGE:" + s);
    search_parameters.setMetaValue("extra_features",  ListUtils::concatenate(percolator_features, ","));
    auto enzyme = *ProteaseDB::getInstance()->getEnzyme(getStringOption_("enzyme"));
    search_parameters.digestion_enzyme = enzyme; // needed for indexing
    search_parameters.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;

    search_parameters.charges = "2:5"; // probably hard-coded in sage https://github.com/lazear/sage/blob/master/crates/sage/src/scoring.rs#L301

    search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
    search_parameters.fixed_modifications = getStringList_("fixed_modifications");
    search_parameters.variable_modifications = getStringList_("variable_modifications");
    search_parameters.missed_cleavages = getIntOption_("missed_cleavages");
    search_parameters.fragment_mass_tolerance = (std::fabs(getDoubleOption_("fragment_tol_left")) + std::fabs(getDoubleOption_("fragment_tol_right"))) * 0.5;
    search_parameters.precursor_mass_tolerance = (std::fabs(getDoubleOption_("precursor_tol_left")) + std::fabs(getDoubleOption_("precursor_tol_right"))) * 0.5;
    search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor_tol_unit") == "ppm";
    search_parameters.fragment_mass_tolerance_ppm = getStringOption_("fragment_tol_unit") == "ppm";

    // write all (!) parameters as metavalues to the search parameters
    if (!protein_identifications.empty())
    {
      DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_identifications[0].getSearchParameters(), this->getToolPrefix());
    }

    // if "reindex" parameter is set to true: will perform reindexing
    if (auto ret = reindex_(protein_identifications, peptide_identifications); ret != EXECUTION_OK) return ret;

    map<String,unordered_map<int,String>> file2specnr2nativeid;
    for (const auto& mzml : input_files)
    {
      // TODO stream mzml?
      MzMLFile m;
      MSExperiment exp;
      auto opts = m.getOptions();
      opts.setMSLevels({2,3});
      opts.setFillData(false);
      //opts.setMetadataOnly(true);
      m.setOptions(opts);
      m.load(mzml, exp);
      String nIDType = "";
      if (!exp.getSourceFiles().empty())
      {
        // TODO we could also guess the regex from the first nativeID if it is not stored here
        //  but I refuse to link to Boost::regex just for this
        //  Someone has to rework the API first!
        nIDType = exp.getSourceFiles()[0].getNativeIDTypeAccession();
      }

      for (const auto& spec : exp)
      {
        const String& nID = spec.getNativeID();
        int nr = SpectrumLookup::extractScanNumber(nID, nIDType);
        if (nr >= 0)
        {
          auto [it, inserted] = file2specnr2nativeid.emplace(File::basename(mzml), unordered_map<int,String>({{nr,nID}}));
          if (!inserted)
          {
            it->second.emplace(nr,nID);
          }
        }
      }
    }

    map<Size, String> idxToFile;
    StringList fnInRun;
    protein_identifications[0].getPrimaryMSRunPath(fnInRun);
    Size cnt = 0;
    for (const auto& f : fnInRun)
    {
      idxToFile.emplace(cnt, f);
      ++cnt;
    }

    for (auto& id : peptide_identifications)
    {
      Int64 scanNrAsInt = 0;

      try
      { // check if spectrum reference is a string that just contains a number
        scanNrAsInt = id.getSpectrumReference().toInt64();
        // no exception -> conversion to int was successful. Now lookup full native ID in corresponding file for given spectrum number.
        // idxToFile values can be full paths but file2specnr2nativeid keys are basenames, so normalize first
        String file_basename = File::basename(idxToFile[id.getMetaValue(Constants::UserParam::ID_MERGE_INDEX)]);
        auto file_it = file2specnr2nativeid.find(file_basename);
        if (file_it != file2specnr2nativeid.end())
        {
          id.setSpectrumReference(file_it->second.at(scanNrAsInt));
        }
      }
      catch (...)
      {
      }
    }

    // Annotate FAIMS compensation voltage if present in any input file
    // Pre-group peptide indices by file for efficient lookup (avoids O(files * peptides))
    std::map<Size, std::vector<Size>> file_to_peptide_indices;
    for (Size i = 0; i < peptide_identifications.size(); ++i)
    {
      const auto& pep = peptide_identifications[i];
      if (pep.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
      {
        file_to_peptide_indices[pep.getMetaValue(Constants::UserParam::ID_MERGE_INDEX)].push_back(i);
      }
    }

    for (const auto& mzml : input_files)
    {
      // Find file index for this mzML
      Size file_idx = 0;
      for (const auto& [idx, fname] : idxToFile)
      {
        if (File::basename(fname) == File::basename(mzml))
        {
          file_idx = idx;
          break;
        }
      }

      // Skip if no peptides for this file
      auto it = file_to_peptide_indices.find(file_idx);
      if (it == file_to_peptide_indices.end() || it->second.empty())
      {
        continue;
      }

      // Load mzML metadata (no peak data needed)
      MzMLFile m;
      MSExperiment exp_full;
      auto opts = m.getOptions();
      opts.setFillData(false);
      m.setOptions(opts);
      m.load(mzml, exp_full);

      // Collect peptide IDs for this file
      PeptideIdentificationList file_peptides;
      file_peptides.reserve(it->second.size());
      for (Size idx : it->second)
      {
        file_peptides.push_back(peptide_identifications[idx]);
      }

      // Annotate FAIMS and copy back
      SpectrumMetaDataLookup::addMissingFAIMSToPeptideIDs(file_peptides, exp_full);
      for (Size i = 0; i < file_peptides.size(); ++i)
      {
        if (file_peptides[i].metaValueExists(Constants::UserParam::FAIMS_CV))
        {
          peptide_identifications[it->second[i]].setMetaValue(
            Constants::UserParam::FAIMS_CV,
            file_peptides[i].getMetaValue(Constants::UserParam::FAIMS_CV));
        }
      }
    }

    IdXMLFile().store(output_file, protein_identifications, peptide_identifications);
    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPPSageAdapter tool;
  return tool.main(argc, argv);
}

/// @endcond
