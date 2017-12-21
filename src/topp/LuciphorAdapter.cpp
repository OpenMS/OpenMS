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
// $Maintainer: Petra Gutenbrunner, Oliver Alka $
// $Authors: Petra Gutenbrunner $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/JavaInfo.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

#include <cstddef>
#include <fstream>
#include <map>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_LuciphorAdapter LuciphorAdapter

   @brief Adapter for the LuciPHOr2: a site localisation tool of generic post-translational modifications from tandem mass spectrometry data.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ LuciphorAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    LuciPHOr2 must be installed before this wrapper can be used. Please make sure that Java and LuciPHOr2 are working.@n
    The following LuciPHOr2 version is required: luciphor2 (JAVA-based version of Luciphor) (1.2014Oct10). At the time of writing, it could be downloaded from http://luciphor2.sourceforge.net.

    Input spectra for LuciPHOr2 have to be in pepXML file format. The input mzML file must be the same as the one used to create the pepXML input file.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_LuciphorAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_LuciphorAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

using namespace OpenMS;
using namespace std;

class LuciphorAdapter :
  public TOPPBase
{
public:
  LuciphorAdapter() :
    TOPPBase("LuciphorAdapter", "Modification site localisation using LuciPHOr2.", true),
    // parameter choices (the order of the values must be the same as in the LuciPHOr2 parameters!):
    fragment_methods_(ListUtils::create<String>("CID,HCD")),
    fragment_error_units_(ListUtils::create<String>("Da,ppm")),
    score_selection_method_(ListUtils::create<String>("Peptide Prophet probability,Mascot Ion Score,-log(E-value),X!Tandem Hyperscore,Sequest Xcorr"))
  {
  }

protected:
  struct LuciphorPSM
  {
    String spec_id;
    int scan_nr;
    int scan_idx;
    int charge;
    String predicted_pep;
    double delta_score;
    double predicted_pep_score;
    
    LuciphorPSM() : scan_nr(-1), scan_idx(-1), charge(-1), delta_score(-1), predicted_pep_score(-1) {}
  };

  // lists of allowed parameter values:
  vector<String> fragment_methods_, fragment_error_units_, score_selection_method_;

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input spectrum file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    
    registerInputFile_("id", "<file>", "", "Protein/peptide identifications file");
    setValidFormats_("id", ListUtils::create<String>("idXML"));

    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerInputFile_("executable", "<file>", "luciphor2.jar", "LuciPHOr2 .jar file, e.g. 'c:\\program files\\luciphor2.jar'", true, false, ListUtils::create<String>("skipexists"));

    registerStringOption_("fragment_method", "<choice>", fragment_methods_[0], "Fragmentation method", false);
    setValidStrings_("fragment_method", fragment_methods_);
    
    registerDoubleOption_("fragment_mass_tolerance", "<value>", 0.5, "Tolerance of the peaks in the fragment spectrum", false);
    registerStringOption_("fragment_error_units", "<choice>", fragment_error_units_[0], "Unit of fragment mass tolerance", false);
    setValidStrings_("fragment_error_units", fragment_error_units_);

    registerDoubleOption_("min_mz", "<value>", 150.0, "Do not consider peaks below this value for matching fragment ions", false);
    
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("target_modifications", "<mods>", ListUtils::create<String>("Phospho (S),Phospho (T),Phospho (Y)"), "List the amino acids to be searched for and their mass modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("target_modifications", all_mods);
    
    registerStringList_("neutral_losses", "<value>", ListUtils::create<String>("sty -H3PO4 -97.97690"), "List the types of neutral losses that you want to consider. The residue field is case sensitive. For example: lower case 'sty' implies that the neutral loss can only occur if the specified modification is present. Syntax: NL = <RESDIUES> -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>", false);
    
    registerDoubleOption_("decoy_mass", "<value>", 79.966331, "How much to add to an amino acid to make it a decoy", false);
    setMinFloat_("decoy_mass", 1.0);
    registerStringList_("decoy_neutral_losses", "<value>", ListUtils::create<String>("X -H3PO4 -97.97690"), "For handling the neutral loss from a decoy sequence. The syntax for this is identical to that of the normal neutral losses given above except that the residue is always 'X'. Syntax: DECOY_NL = X -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>", false);
    
    registerIntOption_("max_charge_state", "<num>", 5, "Do not consider PSMs with a charge state above this value", false);
    setMinInt_("max_charge_state", 1);

    registerIntOption_("max_peptide_length", "<num>", 40, "Restrict scoring to peptides with a length shorter than this value", false);
    setMinInt_("max_peptide_length", 1);

    registerIntOption_("max_num_perm", "<num>", 16384, "Maximum number of permutations a sequence can have", false);
    setMinInt_("max_num_perm", 1);

    registerDoubleOption_("modeling_score_threshold", "<value>", 0.95, "Minimum score a PSM needs to be considered for modeling", false);
    setMinFloat_("modeling_score_threshold", 0.0);
    
    registerDoubleOption_("scoring_threshold", "<value>", 0.0, "PSMs below this value will be discarded", false);
    setMinFloat_("scoring_threshold", 0.0);
    
    registerIntOption_("min_num_psms_model", "<num>", 50, "The minimum number of PSMs you need for any charge state in order to build a model for it", false);
    setMinInt_("min_num_psms_model", 1);

    registerIntOption_("num_threads", "<num>", 6, "For multi-threading, 0 = use all CPU found by JAVA", false);
    setMinInt_("num_threads", 0);

    registerStringOption_("run_mode", "<choice>", "0", "Determines how Luciphor will run: 0 = calculate FLR then rerun scoring without decoys (two iterations), 1 = Report Decoys: calculate FLR but don't rescore PSMs, all decoy hits will be reported", false);
    setValidStrings_("run_mode", ListUtils::create<String>("0,1")); 
    
    registerInputFile_("java_executable", "<file>", "java", "The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java", false, false, ListUtils::create<String>("skipexists"));

    registerIntOption_("java_memory", "<num>", 3500, "Maximum Java heap size (in MB)", false);
    registerIntOption_("java_permgen", "<num>", 0, "Maximum Java permanent generation space (in MB); only for Java 7 and below", false, true);
  }
  
  String makeModString_(const String& mod_name)
  {
    ResidueModification mod = ModificationsDB::getInstance()->getModification(mod_name);
    String residue = mod.getOrigin();
    
    return String(residue +  " " + mod.getDiffMonoMass());    
  }
 
  ExitCodes parseParameters_(map<String, vector<String> >& config_map, const String& id, const String& in,
    const String& out, const vector<String>& target_mods, String selection_method)
  {
    FileHandler fh;
    
    config_map["SPECTRUM_PATH"].push_back(File::path(File::absolutePath(in)));
    config_map["SPECTRUM_SUFFIX"].push_back(FileTypes::typeToName(fh.getTypeByFileName(in)));
    config_map["INPUT_DATA"].push_back(id);
    
    String type = FileTypes::typeToName(fh.getTypeByFileName(id));
    config_map["INPUT_TYPE"].push_back(0);
    
    config_map["ALGORITHM"].push_back(ListUtils::getIndex<String>(fragment_methods_, getStringOption_("fragment_method")));
    config_map["MS2_TOL"].push_back(getDoubleOption_("fragment_mass_tolerance"));
    config_map["MS2_TOL_UNITS"].push_back(ListUtils::getIndex<String>(fragment_error_units_, getStringOption_("fragment_error_units")));
    config_map["MIN_MZ"].push_back(getDoubleOption_("min_mz"));
    config_map["OUTPUT_FILE"].push_back(out);
    config_map["DECOY_MASS"].push_back(getDoubleOption_("decoy_mass"));
    config_map["MAX_CHARGE_STATE"].push_back(getIntOption_("max_charge_state"));
    config_map["MAX_PEP_LEN"].push_back(getIntOption_("max_peptide_length"));
    config_map["MAX_NUM_PERM"].push_back(getIntOption_("max_num_perm"));
    config_map["SELECTION_METHOD"].push_back(ListUtils::getIndex<String>(score_selection_method_, selection_method));    
    config_map["MODELING_SCORE_THRESHOLD"].push_back(getDoubleOption_("modeling_score_threshold"));
    config_map["SCORING_THRESHOLD"].push_back(getDoubleOption_("scoring_threshold"));
    config_map["MIN_NUM_PSMS_MODEL"].push_back(getIntOption_("min_num_psms_model"));
    config_map["NUM_THREADS"].push_back(getIntOption_("num_threads"));
    config_map["RUN_MODE"].push_back(getStringOption_("run_mode"));
    
    for (vector<String>::const_iterator it = target_mods.begin(); it != target_mods.end(); ++it)
    {
      config_map["TARGET_MOD"].push_back(makeModString_(*it));
    }
    
    vector<String> neutral_losses = getStringList_("neutral_losses");
    for (vector<String>::const_iterator it = neutral_losses.begin(); it != neutral_losses.end(); ++it)
    {
      config_map["NL"].push_back(*it);
    }
    
    vector<String> dcy_neutral_losses = getStringList_("decoy_neutral_losses");
    for (vector<String>::const_iterator it = dcy_neutral_losses.begin(); it != dcy_neutral_losses.end(); ++it)
    {
      config_map["DECOY_NL"].push_back(*it);
    }
    
    return EXECUTION_OK;
  }
  
  void writeConfigurationFile_(const String& out_path, map<String, vector<String> >& config_map)
  {
    ofstream output(out_path.c_str());
    output << "## Input file for Luciphor2 (aka: LucXor). (part of OpenMS)\n\n";

    for (std::map<String, vector<String> >::iterator it = config_map.begin(); it != config_map.end(); ++it)
    {
      String key = it->first;
      if (!key.empty())
      {
        for (vector<String>::iterator it_val = it->second.begin(); it_val != it->second.end(); ++it_val)
        {
          output << key << " = " << *it_val << "\n";
        }        
      }
    }

    //------------------------------------------------------------------
    // static parameter definition
    //------------------------------------------------------------------
    // Generate a tab-delimited file of all the matched peaks
    output << "WRITE_MATCHED_PEAKS_FILE = 0\n";
    
    output << "MOD_PEP_REP = 0 ## 0 = show single character modifications, 1 = show TPP-formatted modifications\n";
    
    output << "## This option can be used to help diagnose problems with Luciphor. Multi-threading is disabled in debug mode.\n";
    output << "DEBUG_MODE = 0 ## 0 = default: turn off debugging\n";
    output << "               ## 1 = write peaks selected for modeling to disk\n";
    output << "               ## 2 = write the scores of all permutations for each PSM to disk\n";
    output << "               ## 3 = write the matched peaks for the top scoring permutation to disk\n";
    output << "               ## 4 = write HCD non-parametric models to disk (HCD-mode only option)\n";
  }
  
  void removeTempDir_(const String& temp_dir)
  {
    if (temp_dir.empty()) return; // no temp. dir. created

    if (debug_level_ >= 2)
    {
      writeDebug_("Keeping temporary files in directory '" + temp_dir + "'. Set debug level to 1 or lower to remove them.", 2);
    }
    else
    {
      if (debug_level_ == 1) writeDebug_("Deleting temporary directory '" + temp_dir + "'. Set debug level to 2 or higher to keep it.", 1);
      File::removeDirRecursively(temp_dir);
    }
  }
  
  struct LuciphorPSM splitSpecId_(const String& spec_id)
  {
    struct LuciphorPSM l_psm;
    l_psm.spec_id = spec_id;
    
    vector<String> parts;
    spec_id.split(".", parts);
    l_psm.scan_nr = parts[1].toInt();
    l_psm.charge = parts[3].toInt();
    
    return l_psm;
  }
  
  ExitCodes convertTargetModification_(const vector<String>& target_mods, map<String, String>& modifications)
  {
    modifications.clear();
    for (vector<String>::const_iterator it = target_mods.begin(); it !=target_mods.end(); ++it)
    {
      String mod_param_value = *it;
      String mod;
      
      vector<String> parts;
      mod_param_value.split(' ', parts);
      if (parts.size() != 2)
      {
        writeLog_("Error: cannot parse modification '" + mod_param_value + "'");
        return PARSE_ERROR;
      }
      else
      {
        mod = parts[0];
        String AAs = parts[1];
        
        // LuciPHOr2 discards C-term and N-term modifications in the sequence. The modifications must be added based on the original sequence.
        if (!AAs.hasPrefix("(C-term") && !AAs.hasPrefix("(N-term"))
        {
          AAs.remove(')');
          AAs.remove('(');
          // because origin can be e.g. (STY)
          for (String::iterator aa = AAs.begin(); aa != AAs.end(); ++aa)
          {
            modifications[*aa] = mod;
          }
        }          
      }
    }
    return EXECUTION_OK;
  }
  
  String parseLuciphorOutput_(const String& l_out, map<int, LuciphorPSM>& l_psms, const SpectrumLookup& lookup)
  {
    CsvFile tsvfile(l_out, '\t');
    String spec_id = "";
        
    for (Size row_count = 1; row_count < tsvfile.rowCount(); ++row_count) // skip header line
    {
      vector<String> elements;
      if (!tsvfile.getRow(row_count, elements))
      {
        writeLog_("Error: could not split row " + String(row_count) + " of file '" + l_out + "'");
        return PARSE_ERROR;
      }
      
      spec_id = elements[0];
      struct LuciphorPSM l_psm = splitSpecId_(spec_id);      
      l_psm.scan_idx = lookup.findByScanNumber(l_psm.scan_nr);
      l_psm.predicted_pep = elements[2];
      l_psm.delta_score = elements[7].toDouble();
      l_psm.predicted_pep_score = elements[8].toDouble();
      
      if (l_psms.count(l_psm.scan_idx) > 0)
      {
        return "Duplicate scannr existing " + String(l_psm.scan_nr) + ".";
      }
      l_psms[l_psm.scan_idx] = l_psm;
    }    
    return "";
    
    // String msg = "Spectrum could not be parsed";
    // throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, spec_id, msg);
  }
  
  // remove all modifications which are LuciPHOr2 target modifications,
  // because for these LuciPHOr2 could predict a different position.
  AASequence removeLuciphorTargetMods_(const AASequence& original_seq, const map<String, String>& target_mods_conv)
  {
    if (!original_seq.isModified()) {
      return original_seq;
    }
    
    AASequence seq_converted = AASequence::fromString(original_seq.toUnmodifiedString());
    
    // set C-term/N-term modification
    if (original_seq.hasNTerminalModification())
    {
      seq_converted.setNTerminalModification(original_seq.getNTerminalModificationName());
    }
    if (original_seq.hasCTerminalModification())
    {
      seq_converted.setCTerminalModification(original_seq.getCTerminalModificationName());
    }
    
    // set all modifications, which were not changed by LuciPHOr2
    for (Size i = 0; i < original_seq.size(); ++i)
    {
      if (original_seq.getResidue(i).isModified())
      {
        String mod = original_seq.getResidue(i).getModificationName();
        
        // no target modification, modification can be set
        bool found = false;
        for (map<String, String>::const_iterator iter = target_mods_conv.begin(); iter != target_mods_conv.end() && !found; ++iter)
        {
          if (mod == iter->second)
          {
            found = true;
          }
        }
        if (!found)
        {
          seq_converted.setModification(i, mod);
        }
      }
    }
    return seq_converted;    
  }
  
  // set modifications changed by LuciPHOr2
  ExitCodes setLuciphorTargetMods_(AASequence& seq, String seq_luciphor, const map<String, String>& target_mods_conv)
  {
    for (Size i = 0; i < seq_luciphor.length(); ++i)
    {
      char aa = seq_luciphor[i];
      if (std::islower(aa))
      {
        map<String, String>::const_iterator iter = target_mods_conv.find(String(aa).toUpper());
        if (iter != target_mods_conv.end())
        {
          if (seq.getResidue(i).isModified())
          {
            writeLog_("Error: ambiguous modifications on AA '" + iter->first + "' (" + seq.getResidue(i).getModificationName() + ", " + iter->second + ")");
            return PARSE_ERROR;
          }
          else 
          {
            seq.setModification(i, iter->second);
          }
        }
      }
    }
    return EXECUTION_OK;
  }
  
  void addScoreToMetaValues_(PeptideHit& hit, const String score_type)
  {
    if (!hit.metaValueExists(score_type) && !hit.metaValueExists(score_type + "_score"))
    {
      if (score_type.hasSubstring("score"))
      {
        hit.setMetaValue(score_type, hit.getScore());
      }
      else
      {
        hit.setMetaValue(score_type + "_score", hit.getScore());
      }
    }
  }
  
  String getSelectionMethod_(const PeptideIdentification& pep_id, String search_engine)
  {
    String selection_method = "";
    if (pep_id.getScoreType() == "Posterior Error Probability" || search_engine == "Percolator")
    {
      selection_method = score_selection_method_[0];
    }
    else if (search_engine == "Mascot")
    {
      selection_method = score_selection_method_[1];
    }
    else if (search_engine == "XTandem")
    {
      selection_method = score_selection_method_[3];
    }
    else
    {
      String msg = "SELECTION_METHOD parameter could not be set. Only Mascot, X! Tandem, or Posterior Error Probability score types are supported.";
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg);
    }
    return selection_method;
  }
  
  ExitCodes main_(int, const char**) override
  {
    String java_executable = getStringOption_("java_executable");
    if (!getFlag_("force"))
    {
      if (!JavaInfo::canRun(java_executable))
      {
        writeLog_("Fatal error: Java is needed to run LuciPHOr2!");
        return EXTERNAL_PROGRAM_ERROR;
      }
    }
    else
    {
      writeLog_("The installation of Java was not checked.");
    }

    // create temporary directory
    String temp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString());
    writeDebug_("Creating temporary directory '" + temp_dir + "'", 1);
    QDir d;
    d.mkpath(temp_dir.toQString());

    // create a temporary config file for LuciPHOr2 parameters
    String conf_file = temp_dir + "luciphor2_input_template.txt";
    
    String id = getStringOption_("id");
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    
    FileHandler fh;
    FileTypes::Type in_type = fh.getType(id);

    vector<PeptideIdentification> pep_ids;
    vector<ProteinIdentification> prot_ids;

    PeakMap exp;
    MzMLFile file;
    file.setLogType(log_type_);
    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);

    file.load(in, exp);
    exp.sortSpectra(true);

    // convert input to pepXML if necessary
    if (in_type == FileTypes::IDXML)
    {
      IdXMLFile().load(id, prot_ids, pep_ids);
      IDFilter::keepNBestHits(pep_ids, 1); // LuciPHOR2 only calculates the best hit
      
      // create a temporary pepXML file for LuciPHOR2 input
      String id_file_name = File::removeExtension(File::basename(id));
      id = temp_dir + id_file_name + ".pepXML";
      
      PepXMLFile().store(id, prot_ids, pep_ids, in, "", false);
    }
    else
    {
      writeLog_("Error: Unknown input file type given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    
    vector<String> target_mods = getStringList_("target_modifications");
    if (target_mods.empty())
    {
      writeLog_("Error: No target modification existing.");
      return ILLEGAL_PARAMETERS;
    }
    
    // initialize map
    map<String, vector<String> > config_map;
    String selection_method = getSelectionMethod_(pep_ids[0], prot_ids.begin()->getSearchEngine());
    
    ExitCodes ret = parseParameters_(config_map, id, in, out, target_mods, selection_method);
    if (ret != EXECUTION_OK)
    {
      return ret;
    }
    
    writeConfigurationFile_(conf_file, config_map);    

    // memory for JVM
    QString java_memory = "-Xmx" + QString::number(getIntOption_("java_memory")) + "m";
    int java_permgen = getIntOption_("java_permgen");

    QString executable = getStringOption_("executable").toQString();
    
    // Hack for KNIME. Looks for LUCIPHOR_PATH in the environment which is set in binaries.ini
    QProcessEnvironment env;
    String luciphorpath = "LUCIPHOR_PATH";
    QString qluciphorpath = env.systemEnvironment().value(luciphorpath.toQString());

    if (!qluciphorpath.isEmpty())
    {
      executable = qluciphorpath;
    }

    QStringList process_params; // the actual process is Java, not LuciPHOr2!
    process_params << java_memory;
    
    if (java_permgen > 0)
    {
      process_params << "-XX:MaxPermSize=" + QString::number(java_permgen);
    }

    process_params << "-jar" << executable << conf_file.toQString();                   
    // execute LuciPHOr2    
    int status = QProcess::execute(java_executable.toQString(), process_params);
    if (status != 0)
    {
      writeLog_("Fatal error: Running LuciPHOr2 returned an error code. Does the LuciPHOr2 executable (.jar file) exist?");
      return EXTERNAL_PROGRAM_ERROR;
    }

    SpectrumLookup lookup;
    lookup.rt_tolerance = 0.05;
    lookup.readSpectra(exp.getSpectra());
      
    map<int, LuciphorPSM> l_psms;    
    ProteinIdentification::SearchParameters search_params;
    
    String error = parseLuciphorOutput_(out, l_psms, lookup);
    if (error != "")
    {
      error = "Error: LuciPHOr2 output is not correctly formated. " + error;
      writeLog_(error);
      return PARSE_ERROR;
    }
    
    //-------------------------------------------------------------
    // writing output - merge LuciPHOr2 result to idXML
    //-------------------------------------------------------------
    vector<PeptideIdentification> pep_out;
    map<String, String> target_mods_conv;
    ret = convertTargetModification_(target_mods, target_mods_conv);
    if (ret != EXECUTION_OK)
    {
      return ret;
    }
    
    for (vector<PeptideIdentification>::iterator pep_id = pep_ids.begin(); pep_id != pep_ids.end(); ++pep_id)
    {
      Size scan_idx = lookup.findByRT(pep_id->getRT());
      
      vector<PeptideHit> scored_peptides;
      if (!pep_id->getHits().empty())
      {
        PeptideHit scored_hit = pep_id->getHits()[0];
        addScoreToMetaValues_(scored_hit, pep_id->getScoreType());
        
        struct LuciphorPSM l_psm;
        if (l_psms.count(scan_idx) > 0)
        {
          l_psm = l_psms.at(scan_idx);
          AASequence original_seq = scored_hit.getSequence();
          
          AASequence predicted_seq = removeLuciphorTargetMods_(original_seq, target_mods_conv);          
          ret = setLuciphorTargetMods_(predicted_seq, l_psm.predicted_pep, target_mods_conv);
          if (ret != EXECUTION_OK)
          {
            return ret;
          }
          scored_hit.setMetaValue("search_engine_sequence", scored_hit.getSequence().toString());
          scored_hit.setMetaValue("Luciphor_pep_score", l_psm.predicted_pep_score);
          scored_hit.setScore(l_psm.delta_score);
          scored_hit.setSequence(predicted_seq);
        }
        else
        {
          scored_hit.setScore(-1);
        }
        scored_peptides.push_back(scored_hit);
      }
      else
      {
        writeLog_("Error: LuciPHOr2 output does not match with idXML.");
        return PARSE_ERROR;
      }
      
      PeptideIdentification new_pep_id(*pep_id);
      new_pep_id.setScoreType("Luciphor_delta_score");
      new_pep_id.setHigherScoreBetter(true);
      new_pep_id.setHits(scored_peptides);
      new_pep_id.assignRanks();
      pep_out.push_back(new_pep_id);
    }
    IdXMLFile().store(out, prot_ids, pep_out);

    removeTempDir_(temp_dir);

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  LuciphorAdapter tool;
  return tool.main(argc, argv);
}
