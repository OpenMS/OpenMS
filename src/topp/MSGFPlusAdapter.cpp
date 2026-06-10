// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Dilek Dere, Mathias Walzer, Petra Gutenbrunner, Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/JavaInfo.h>

#include <boost/interprocess/sync/file_lock.hpp>

#include <algorithm>
#include <fstream>
#include <map>
#include <cstddef>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MSGFPlusAdapter MSGFPlusAdapter

@brief Adapter for the MS-GF+ protein identification (database search) engine.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; MSGFPlusAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes @n (or another centroiding tool)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

MS-GF+ must be installed before this wrapper can be used. Please make sure that Java and MS-GF+ are working.@n
At the time of writing, MS-GF+ can be downloaded from https://github.com/MSGFPlus/msgfplus/releases.

The following MS-GF+ version is required: <b>MS-GF+ 2019/07/03</b>. Older versions will not work properly, giving
an error: <em>[Error] Invalid parameter: -maxMissedCleavages.</em>

Input spectra for MS-GF+ have to be centroided; profile spectra will raise an error in the adapter.

The first time MS-GF+ is applied to a database (FASTA file), it will index the file contents and
generate a number of auxiliary files in the same directory as the database (e.g. for "db.fasta": "db.canno", "db.cnlap", "db.csarr" and "db.cseq" will be generated).
It is advisable to keep these files for future MS-GF+ searches, to save the indexing step.@n

@note This Adapter uses an internal locking mechanism (a file lock), to ensure that MSGF+ does not attempt to create the database index
in parallel (which would fail badly) when multiple instances of this Adapter are run concurrently on the same FASTA database.
After the database has been indexed, multiple MS-GF+ processes (even without this Adapters locking) can use it in parallel.

This adapter supports relative database filenames, which (when not found in the current working directory) are looked up in the directories specified 
by 'OpenMS.ini:id_db_dir'.

The adapter works in three steps to generate an idXML file: First MS-GF+ is run on the input MS data and the sequence database, 
producing an mzIdentML (.mzid) output file containing the search results. This file is then converted to a text file (.tsv) using MS-GF+' "MzIDToTsv" tool.
Finally, the .tsv file is parsed and a result in idXML format is generated.

An optional MSGF+ configuration file can be added via '-conf' parameter.
See https://github.com/MSGFPlus/msgfplus/blob/master/docs/examples/MSGFPlus_Params.txt for 
an example and consult the MSGF+ documentation for further details.
Parameters specified in the configuration file are ignored by MS-GF+ if they are also specified on the command line.
This adapter passes all flags which you can set on the command line, so use the configuration file <b>only</b> for parameters which
are not available here (this includes fixed/variable modifications, which are passed on the commandline via <code>-mod &lt;file&gt;</code>).
Thus, be very careful that your settings in '-conf' actually take effect (try running again without '-conf' file and test if the results change).

@note This adapter supports 15N labeling by specifying the 20 AA modifications 'Label:15N(x)' as fixed modifications.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MSGFPlusAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MSGFPlusAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

using namespace OpenMS;
using namespace std;

class MSGFPlusAdapter :
  public SearchEngineBase
{
public:
  MSGFPlusAdapter() :
    SearchEngineBase("MSGFPlusAdapter", "MS/MS database search using MS-GF+.", true),
    // parameter choices (the order of the values must be the same as in the MS-GF+ parameters!):
    fragment_methods_(ListUtils::create<std::string>("from_spectrum,CID,ETD,HCD")),
    instruments_(ListUtils::create<std::string>("low_res,high_res,TOF,Q_Exactive")),
    protocols_(ListUtils::create<std::string>("automatic,phospho,iTRAQ,iTRAQ_phospho,TMT,none")),
    tryptic_(ListUtils::create<std::string>("non,semi,fully"))
  {
    ProteaseDB::getInstance()->getAllMSGFNames(enzymes_);
    std::sort(enzymes_.begin(),enzymes_.end());
  }

protected:
  /// parts of a sequence of the form "K.AAAA.R"
  struct SequenceParts
  {
    char aa_before, aa_after; // may be '\0' if not given
    std::string peptide;

    SequenceParts(): aa_before(0), aa_after(0) {}
  };

  // lists of allowed parameter values:
  vector<std::string> fragment_methods_, instruments_, enzymes_, protocols_, tryptic_;

  // primary MS run referenced in the mzML file
  StringList primary_ms_run_path_;

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (MS-GF+ parameter '-s')");
    setValidFormats_("in", {"mzML", "mzXML", "mgf", "ms2" });
    registerOutputFile_("out", "<file>", "", "Output file (.idXML) or directory bundle (.idparquet) containing the search results.", false);
    setValidFormats_("out", {"idXML", "idparquet"});
    registerOutputFile_("mzid_out", "<file>", "", "Alternative output file (MS-GF+ parameter '-o')\nEither 'out' or 'mzid_out' are required. They can be used together.", false);
    setValidFormats_("mzid_out", ListUtils::create<std::string>("mzid"));
    registerInputFile_("executable", "<file>", "MSGFPlus.jar", "The MSGFPlus Java archive file. Provide a full or relative path, or make sure it can be found in your PATH environment.", true, false, {"is_executable"});
    registerInputFile_("database", "<file>", "", "Protein sequence database (FASTA file; MS-GF+ parameter '-d'). Non-existing relative filenames are looked up via 'OpenMS.ini:id_db_dir'.", true, false, ListUtils::create<std::string>("skipexists"));
    setValidFormats_("database", ListUtils::create<std::string>("FASTA"));

    registerDoubleOption_("precursor_mass_tolerance", "<value>", 10, "Precursor monoisotopic mass tolerance (MS-GF+ parameter '-t')", false);
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "Unit of precursor mass tolerance (MS-GF+ parameter '-t')", false);
    setValidStrings_("precursor_error_units", ListUtils::create<std::string>("Da,ppm"));

    registerStringOption_("isotope_error_range", "<range>", "0,1", "Range of allowed isotope peak errors (MS-GF+ parameter '-ti'). Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation. Combined with 'precursor_mass_tolerance'/'precursor_error_units', this determines the actual precursor mass tolerance. E.g. for experimental mass 'exp' and calculated mass 'calc', '-precursor_mass_tolerance 20 -precursor_error_units ppm -isotope_error_range -1,2' tests '|exp - calc - n * 1.00335 Da| < 20 ppm' for n = -1, 0, 1, 2.", false);

    registerStringOption_("fragment_method", "<choice>", fragment_methods_[0], "Fragmentation method ('from_spectrum' relies on spectrum meta data and uses CID as fallback option; MS-GF+ parameter '-m')", false);
    setValidStrings_("fragment_method", fragment_methods_);

    registerStringOption_("instrument", "<choice>", instruments_[0], "Instrument that generated the data ('low_res'/'high_res' refer to LCQ and LTQ instruments; MS-GF+ parameter '-inst')", false);
    setValidStrings_("instrument", instruments_);

    registerStringOption_("enzyme", "<choice>", enzymes_[6], "Enzyme used for digestion, or type of cleavage. Note: MS-GF+ does not support blocking rules. (MS-GF+ parameter '-e')", false);
    setValidStrings_("enzyme", enzymes_);

    registerStringOption_("protocol", "<choice>", protocols_[0], "Labeling or enrichment protocol used, if any (MS-GF+ parameter '-p')", false);
    setValidStrings_("protocol", protocols_);

    registerStringOption_("tryptic", "<choice>", tryptic_[2], "Level of cleavage specificity required (MS-GF+ parameter '-ntt')", false);
    setValidStrings_("tryptic", tryptic_);

    registerIntOption_("min_precursor_charge", "<num>", 2, "Minimum precursor ion charge (only used for spectra without charge information; MS-GF+ parameter '-minCharge')", false);
    setMinInt_("min_precursor_charge", 1);
    registerIntOption_("max_precursor_charge", "<num>", 3, "Maximum precursor ion charge (only used for spectra without charge information; MS-GF+ parameter '-maxCharge')", false);
    setMinInt_("max_precursor_charge", 1);

    registerIntOption_("min_peptide_length", "<num>", 6, "Minimum peptide length to consider (MS-GF+ parameter '-minLength')", false);
    setMinInt_("min_peptide_length", 1);
    registerIntOption_("max_peptide_length", "<num>", 40, "Maximum peptide length to consider (MS-GF+ parameter '-maxLength')", false);
    setMinInt_("max_peptide_length", 1);

    registerIntOption_("matches_per_spec", "<num>", 1, "Number of matches per spectrum to be reported (MS-GF+ parameter '-n')", false);
    setMinInt_("matches_per_spec", 1);

    registerIntOption_("min_peaks", "<num>", 10, "Minimum number of ions a spectrum must have to be examined", false); 
    setMinInt_("min_peaks", 10); 

    registerStringOption_("add_features", "<true/false>", "true", "Output additional features (MS-GF+ parameter '-addFeatures'). This is required by Percolator and hence by default enabled.", false, false);
    setValidStrings_("add_features", ListUtils::create<std::string>("true,false"));
    
    registerIntOption_("max_mods", "<num>", 2, "Maximum number of modifications per peptide. If this value is large, the search may take very long.", false);
    setMinInt_("max_mods", 0);

    registerIntOption_("max_missed_cleavages", "<num>", -1, "Maximum number of missed cleavages allowed for a peptide to be considered for scoring. (default: -1 meaning unlimited)", false);
    setMinInt_("max_missed_cleavages", -1);

    registerIntOption_("tasks", "<num>", 0, "(Override the number of tasks to use on the threads; Default: (internally calculated based on inputs))\n"
                                             "   More tasks than threads will reduce the memory requirements of the search, but will be slower (how much depends on the inputs).\n"
                                             "   1 <= tasks <= numThreads: will create one task per thread, which is the original behavior.\n"
                                             "   tasks = 0: use default calculation - minimum of: (threads*3) and (numSpectra/250).\n"
                                             "   tasks < 0: multiply number of threads by abs(tasks) to determine number of tasks (i.e., -2 means \"2 * numThreads\" tasks).\n"
                                             "   One task per thread will use the most memory, but will usually finish the fastest.\n"
                                             "   2-3 tasks per thread will use comparably less memory, but may cause the search to take 1.5 to 2 times as long.", false);

    vector<std::string> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", {"Carbamidomethyl (C)"}, "Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", {"Oxidation (M)"}, "Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'",
                        false);
    setValidStrings_("variable_modifications", all_mods);

    registerFlag_("allow_dense_centroided_peaks", "Pass '-allowDenseCentroidedPeaks 1' to MS-GF+ to also search spectra whose median distance of consecutive peaks is below ~50 ppm. The MS-GF+ default treats such (centroided) spectra as profile spectra and skips them.", false);

    registerFlag_("legacy_conversion", "Use the indirect conversion of MS-GF+ results to idXML via export to TSV. Try this only if the default conversion takes too long or uses too much memory.", true);

    registerInputFile_("conf", "<file>", "", "Optional MSGF+ configuration file (passed as -conf <file> to MSGF+). See documentation for examples. Parameters of the adapter take precedence. Use conf file only for settings not available here (for example, any fixed/var modifications, in the conf file will be ignored, since they are provided via -mod flag)", false, false);

    registerInputFile_("java_executable", "<file>", "java", "The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java", false, false, {"is_executable"});
    registerIntOption_("java_memory", "<num>", 3500, "Maximum Java heap size (in MB)", false);
    registerIntOption_("java_permgen", "<num>", 0, "Maximum Java permanent generation space (in MB); only for Java 7 and below", false, true);

    // register peptide indexing parameter (with defaults for this search engine) TODO: check if search engine defaults are needed
    registerPeptideIndexingParameter_(PeptideIndexing().getParameters()); 
  }

  // The following sequence modification methods are used to modify the sequence stored in the TSV such that it can be used by AASequence

  // Method to cut the amino acids before/after the peptide (splice sites) off the sequence.
  // The sequences in the TSV file have the format 'K.XXXR.X' (where XXXR is the actual peptide sequence).
  // This method returns the sequence split into its three parts (e.g. "K", "XXXR", "X").
  struct SequenceParts splitSequence_(const std::string& sequence)
  {
    struct SequenceParts parts;
    size_t len = sequence.size(), start = 0, count = string::npos;
    if (len > 3) // in 'X.Y', which side would we cut off?
    {
      if (sequence[1] == '.')
      {
        start = 2;
        parts.aa_before = sequence[0];
      }
      if (sequence[len - 2] == '.')
      {
        count = len - start - 2;
        parts.aa_after = sequence[len - 1];
      }
    }
    parts.peptide = StringUtils::substr(sequence, start, count);
    return parts;
  }

  std::string modifyNTermAASpecificSequence_(const std::string& seq)
  {
    std::string swap;
    string modifiedSequence = seq;
    vector<pair<std::string, char> > massShiftList;

    massShiftList.push_back(make_pair("-18.011", 'E'));
    massShiftList.push_back(make_pair("-17.027", 'Q'));

    for (vector<pair<std::string, char> >::const_iterator it = massShiftList.begin(); it != massShiftList.end(); ++it)
    {
      string modMassShift = it->first;
      size_t found = modifiedSequence.find(modMassShift);

      if (found != string::npos)
      {
        std::string tmp = StringUtils::substr(modifiedSequence, 0, found + modMassShift.length() + 1);
        size_t foundAA = tmp.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");

        if ((foundAA > found) && (tmp[foundAA] == it->second)) // no AA at the begin
        {
          if (found > 0)
          {
            swap = StringUtils::substr(modifiedSequence, 0, found);
          }
          return swap += *tmp.rbegin() + modMassShift + StringUtils::substr(modifiedSequence, found + modMassShift.length() + 1);
        }
      }
    }
    return  modifiedSequence;
  }

  // Method to replace the mass representation of modifications.
  // Modifications in the TSV file have the format 'M+15.999'
  // After using this method the sequence should look like this: 'M[+15.999]'
  std::string modifySequence_(const std::string& seq)
  {
    std::string modifiedSequence = seq;
    size_t found1 = modifiedSequence.find_first_of("+-");
    while (found1 != string::npos)
    {
      modifiedSequence = modifiedSequence.insert(found1, 1, '[');
      size_t found2 = modifiedSequence.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ", found1);
      if (found2 != string::npos)
      {
        modifiedSequence.insert(found2, 1, ']');
        found1 = modifiedSequence.find_first_of("+-", found2 + 2);
      }
      else // last amino acid is modified
      {
        modifiedSequence = modifiedSequence + ']';
        return modifiedSequence;
      }
    }
    return modifiedSequence;
  }

  // Parse mzML and create RTMapping
  // get RT: it doesn't exist in output from MS-GF+
  // get m/z: it is rounded after converting to TSV
  void generateInputfileMapping_(map<std::string, vector<float> >& rt_mapping)
  {
    std::string exp_name = getStringOption_("in");

    if (!exp_name.empty())
    {
      PeakMap exp;
      // load only MS2 spectra:
      FileHandler f;
      f.getOptions().addMSLevel(2);
      f.getOptions().setFillData(false);
      f.loadExperiment(exp_name, exp, {FileTypes::MZML});
      exp.getPrimaryMSRunPath(primary_ms_run_path_);
      // if no primary run is assigned, the mzML file is the (unprocessed) primary file
      if (primary_ms_run_path_.empty())
      {
        primary_ms_run_path_.push_back(exp_name);
      }

      for (MSSpectrum& ms : exp)
      {
        std::string id = ms.getNativeID(); // expected format: "... scan=#"
        if (!id.empty())
        {
          rt_mapping[id].push_back(ms.getRT());
          rt_mapping[id].push_back(ms.getPrecursors()[0].getMZ());
        }
      }
    }
  }

  std::string makeModString_(const std::string& mod_name, bool fixed=true)
  {
    const ResidueModification* mod = ModificationsDB::getInstance()->getModification(mod_name);
    char residue = mod->getOrigin();
    if (residue == 'X')
    {
      residue = '*'; // terminal mod. without residue specificity
    }
    std::string position = mod->getTermSpecificityName();
    if (position == "Protein N-term")
    {
      position = "Prot-N-term";
    }
    else if (position == "Protein C-term")
    {
      position = "Prot-C-term";
    }
    else if (position == "none")
    {
      position = "any";
    }
    return StringUtils::toStr(mod->getDiffMonoMass()) + ", " + residue + (fixed ? ", fix, " : ", opt, ") + position + ", " + mod->getId() + "    # " + mod_name;
  }

  void writeModificationsFile_(const std::string& out_path, const vector<std::string>& fixed_mods, const vector<std::string>& variable_mods, Size max_mods)
  {
    ofstream output(out_path.c_str());
    if (!output)
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_path);
    }
    output << "# MS-GF+ modifications file written by MSGFPlusAdapter (part of OpenMS)\n"
           << "NumMods=" << max_mods
           << "\n\n# Fixed modifications:\n";
    if (fixed_mods.empty())
    {
      output << "# (none)\n";
    }
    else
    {
      for (vector<std::string>::const_iterator it = fixed_mods.begin(); it != fixed_mods.end(); ++it)
      {
        output << makeModString_(*it) << "\n";
      }
    }
    output << "\n# Variable modifications:\n";
    if (variable_mods.empty())
    {
      output << "# (none)\n";
    }
    else
    {
      for (vector<std::string>::const_iterator it = variable_mods.begin(); it != variable_mods.end(); ++it)
      {
        output << makeModString_(*it, false) << "\n";
      }
    }
  }

  std::string describeHit_(const PeptideHit& hit)
  {
    return "peptide hit with sequence '" + hit.getSequence().toString() +
      "', charge " + StringUtils::toStr(hit.getCharge()) + ", score " +
      StringUtils::toStr(hit.getScore());
  }

  // Set the MS-GF+ e-value (MS:1002052) as new peptide identification score.
  void switchScores_(PeptideIdentification& id)
  {
    for (PeptideHit& hit : id.getHits())
    {
      // MS:1002052 == MS-GF spectral E-value
      if (!hit.metaValueExists("MS:1002052"))
      {
        std::string msg = "Meta value 'MS:1002052' not found for " + describeHit_(hit);
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg);
      }

      hit.setScore(hit.getMetaValue("MS:1002052"));
    }
    id.setScoreType("SpecEValue");
    id.setHigherScoreBetter(false);
  }

  bool createLockedDBIndex(const std::string& db_name, const std::string& java_executable, const std::string& java_memory, const std::string& executable)
  {
    const std::string db_indexfile = FileHandler::stripExtension(db_name) + ".canno";
    const std::string lockfile_path = db_name + ".lock";
    OPENMS_LOG_DEBUG << "Checking for db index, using a lock file ..." << std::endl;

    // Create the lock file if it doesn't exist (needed by boost::interprocess::file_lock)
    {
      std::ofstream touch(lockfile_path.c_str(), std::ios::app);
      if (!touch.is_open())
      {
        // Cannot create lock file - check if index already exists
        std::string msg = "The lock file could not be created, for lack of permissions in the parent directory.";
        if (!File::exists(db_indexfile))
        {
          OPENMS_LOG_ERROR << msg << " Checking index anyway: No database index found! Please make the directory writable or pre-create an DB index." << std::endl;
          return false;
        }
        OPENMS_LOG_DEBUG << msg << " Checking index anyway: found it!" << std::endl;
        return true;
      }
    }

    try
    {
      boost::interprocess::file_lock lock_db(lockfile_path.c_str());
      if (!lock_db.try_lock())
      {
        // Could not acquire lock - another process holds it. Wait for it (blocking).
        lock_db.lock();
      }

      // we have a lock: now check if we need to create a new index (which only one instance should do)
      if (!File::exists(db_indexfile))
      {
        OPENMS_LOG_INFO << "\nNo database index found! Creating index while holding a lock ..." << std::endl;
        std::vector<std::string> process_params = {java_memory, "-cp", executable, "edu.ucsd.msjava.msdbsearch.BuildSA", "-d", db_name, "-tda", "0"};

        // collect all output since MSGF+ might return 'success' even though it did not like the command arguments (e.g. if the version is too old)
        // If no output file is produced, we can print the stderr below.
        std::string proc_stdout, proc_stderr;

        TOPPBase::ExitCodes exit_code = runExternalProcess_(java_executable, process_params, proc_stdout, proc_stderr);
        if (exit_code != EXECUTION_OK)
        {
          // if there was sth like a segfault, runExternalProcess_ will write a warning about the type of error,
          //  but not print the output of the program.
          OPENMS_LOG_ERROR << "The output of MSGF+'s Index Database Creation was:\nSTDOUT:\n" << proc_stdout << "\nSTDERR:\n" << proc_stderr << endl;
          lock_db.unlock();
          return false;
        }
        // Verify index was actually created (BuildSA may return success without producing output)
        if (!File::exists(db_indexfile))
        {
          OPENMS_LOG_ERROR << "MSGF+ BuildSA reported success but index file '" << db_indexfile << "' was not created.\nSTDOUT:\n" << proc_stdout << "\nSTDERR:\n" << proc_stderr << endl;
          lock_db.unlock();
          return false;
        }
        OPENMS_LOG_INFO << " ... done" << std::endl;
      }

      // free lock, since database index exists at this point
      lock_db.unlock();
      OPENMS_LOG_DEBUG << "... releasing DB lock" << std::endl;
      return true;
    }
    catch (const boost::interprocess::interprocess_exception& e)
    {
      OPENMS_LOG_ERROR << "An error occurred while trying to acquire a file lock using the file '" << lockfile_path
                       << "': " << e.what()
                       << ".\nPlease check the previous error message and contact OpenMS support if you cannot solve the problem.";
      return false;
    }
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parse parameters
    //-------------------------------------------------------------

    std::string in = getRawfileName();
    std::string out = getStringOption_("out");
    std::string mzid_out = getStringOption_("mzid_out");
    if (mzid_out.empty() && out.empty())
    {
      writeLogError_("Error:  no output file given (parameter 'out' or 'mzid_out')");
      return ILLEGAL_PARAMETERS;
    }

    const std::string java_executable = getStringOption_("java_executable");
    if (!getFlag_("force"))
    {
      if (!JavaInfo::canRun(java_executable))
      {
        writeLogError_("Fatal error: Java is needed to run MS-GF+!");
        return EXTERNAL_PROGRAM_NOTFOUND;
      }
    }
    else
    {
      writeLogWarn_("The installation of Java was not checked.");
    }
    
    const std::string java_memory = "-Xmx" + StringUtils::toStr(getIntOption_("java_memory")) + "m";
    const std::string executable = getStringOption_("executable");

    const std::string db_name = getDBFilename();
    if (!createLockedDBIndex(db_name, java_executable, java_memory, executable))
    {
      OPENMS_LOG_ERROR << "Could not create/verify database index. Aborting ..." << std::endl;
      return ExitCodes::INTERNAL_ERROR;
    }

    vector<std::string> fixed_mods = getStringList_("fixed_modifications");
    vector<std::string> variable_mods = getStringList_("variable_modifications");
    bool no_mods = fixed_mods.empty() && variable_mods.empty();
    Int max_mods = getIntOption_("max_mods");
    if ((max_mods == 0) && !no_mods)
    {
      writeLogWarn_("Warning: Modifications are defined ('fixed_modifications'/'variable_modifications'), but the number of allowed modifications is zero ('max_mods'). Is that intended?");
    }

    // create temporary directory (and modifications file, if necessary):
    File::TempDir tmp_dir(debug_level_ >= 2);
    std::string mzid_temp, mod_file;
    // always create a temporary mzid file first, even if mzid output is requested via "mzid_out"
    // (reason: TOPPAS may pass a filename with wrong extension to "mzid_out", which would cause an error in MzIDToTSVConverter below,
    // so we make sure that we have a properly named mzid file for the converter; see https://github.com/OpenMS/OpenMS/issues/1251)
    mzid_temp = tmp_dir.getPath() + "msgfplus_output.mzid";
    if (!no_mods)
    {
      mod_file = tmp_dir.getPath() + "msgfplus_mods.txt";
      writeModificationsFile_(mod_file, fixed_mods, variable_mods, max_mods);
    }

    // parameters also used by OpenMS (see idXML creation below):
    std::string enzyme = getStringOption_("enzyme");
    double precursor_mass_tol = getDoubleOption_("precursor_mass_tolerance");
    std::string precursor_error_units = getStringOption_("precursor_error_units");
    Int min_precursor_charge = getIntOption_("min_precursor_charge");
    Int max_precursor_charge = getIntOption_("max_precursor_charge");
    // parameters only needed for MS-GF+:
    // no need to handle "not found" case - would have given error during parameter parsing:
    Int fragment_method_code = ListUtils::getIndex<std::string>(fragment_methods_, getStringOption_("fragment_method"));
    Int instrument_code = ListUtils::getIndex<std::string>(instruments_, getStringOption_("instrument"));
    Int enzyme_code = ProteaseDB::getInstance()->getEnzyme(enzyme)->getMSGFID();
    Int protocol_code = ListUtils::getIndex<std::string>(protocols_, getStringOption_("protocol"));
    // protocol code = 0 corresponds to "automatic" (MS-GF+ docu 2017) and "none" (MS-GF+ docu 2013). We keep 0 = "none" for backward compatibility.
    if (protocol_code == 5)
    {
        protocol_code = 0;
    }
    Int tryptic_code = ListUtils::getIndex<std::string>(tryptic_, getStringOption_("tryptic"));

    std::vector<std::string> process_params = { // the actual process is Java, not MS-GF+!
      java_memory,
      "-jar", executable,
      "-s", in,
      "-o", mzid_temp,
      "-d", db_name,
      "-t",StringUtils::toStr(precursor_mass_tol) + precursor_error_units,
      "-ti", getStringOption_("isotope_error_range"),
      "-m",StringUtils::toStr(fragment_method_code),
      "-inst",StringUtils::toStr(instrument_code),
      "-e",StringUtils::toStr(enzyme_code),
      "-protocol",StringUtils::toStr(protocol_code),
      "-ntt",StringUtils::toStr(tryptic_code),
      "-minLength",StringUtils::toStr(getIntOption_("min_peptide_length")),
      "-maxLength",StringUtils::toStr(getIntOption_("max_peptide_length")),
      "-minNumPeaks",StringUtils::toStr(getIntOption_("min_peaks")),
      "-minCharge",StringUtils::toStr(min_precursor_charge),
      "-maxCharge",StringUtils::toStr(max_precursor_charge),
      "-maxMissedCleavages",StringUtils::toStr(getIntOption_("max_missed_cleavages")),
      "-n",StringUtils::toStr(getIntOption_("matches_per_spec")),
      "-addFeatures",StringUtils::toStr(int((getParam_().getValue("add_features") == "true"))),
      "-tasks",StringUtils::toStr(getIntOption_("tasks")),
      "-thread",StringUtils::toStr(getIntOption_("threads"))
    };
    if (getFlag_("allow_dense_centroided_peaks"))
    {
      // pass flag and value as separate tokens (each vector element becomes one argv entry)
      process_params.push_back("-allowDenseCentroidedPeaks");
      process_params.push_back("1");
    }
    std::string conf = getStringOption_("conf");
    if (!conf.empty())
    {
      process_params.push_back("-conf"); process_params.push_back(conf);
    }

    if (!mod_file.empty())
    {
      process_params.push_back("-mod"); process_params.push_back(mod_file);
    }

    //-------------------------------------------------------------
    // execute MS-GF+
    //-------------------------------------------------------------

    // run MS-GF+ process and create the .mzid file

    writeLogInfo_("Running MSGFPlus search...");
    // collect all output since MSGF+ might return 'success' even though it did not like the command arguments (e.g. if the version is too old)
    // If no output file is produced, we can print the stderr below.
    std::string proc_stdout, proc_stderr;

    TOPPBase::ExitCodes exit_code = runExternalProcess_(java_executable, process_params, proc_stdout, proc_stderr);
    if (exit_code != EXECUTION_OK)
    {
      // if there was sth like a segfault, runExternalProcess_ will write a warning about the type of error,
      //  but not print the output of the program.
      OPENMS_LOG_ERROR << "The output of MSGF+ was:\nSTDOUT:\n" << proc_stdout << "\nSTDERR:\n" << proc_stderr << endl;
      return exit_code;
    }

    //-------------------------------------------------------------
    // create idXML output
    //-------------------------------------------------------------
    if (!out.empty())
    {
      if (!File::exists(mzid_temp))
      {
        OPENMS_LOG_ERROR << "MSGF+ failed. Temporary output file '" << mzid_temp << "' was not created.\n"
                         << "The output of MSGF+ was:\nSTDOUT:\n" << proc_stdout << "\nSTDERR:\n" << proc_stderr << endl;
        return EXTERNAL_PROGRAM_ERROR;
      }

      vector<ProteinIdentification> protein_ids;
      PeptideIdentificationList peptide_ids;

      if (getFlag_("legacy_conversion"))
      {
        // run TSV converter
        std::string tsv_out = tmp_dir.getPath() + "msgfplus_converted.tsv";
        int java_permgen = getIntOption_("java_permgen");
        process_params.clear();
        process_params.push_back(java_memory);
        if (java_permgen > 0)
        {
          process_params.push_back("-XX:MaxPermSize=" + StringUtils::toStr(java_permgen) + "m");
        }
        process_params.insert(process_params.end(), {"-cp", executable, "edu.ucsd.msjava.ui.MzIDToTsv",
                       "-i", mzid_temp,
                       "-o", tsv_out,
                       "-showQValue", "1",
                       "-showDecoy", "1",
                       "-unroll", "1"});
        writeLogInfo_("Running MzIDToTSVConverter...");
        exit_code = runExternalProcess_(java_executable, process_params);
        if (exit_code != EXECUTION_OK)
        {
          return exit_code;
        }

        // initialize map
        map<std::string, vector<float> > rt_mapping;
        generateInputfileMapping_(rt_mapping);

        // handle the search parameters
        ProteinIdentification::SearchParameters search_parameters;
        search_parameters.db = db_name;
        search_parameters.charges = "+" + StringUtils::toStr(min_precursor_charge) + "-+" + StringUtils::toStr(max_precursor_charge);
        search_parameters.mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
        search_parameters.fixed_modifications = fixed_mods;
        search_parameters.variable_modifications = variable_mods;
        search_parameters.precursor_mass_tolerance = precursor_mass_tol;
        search_parameters.precursor_mass_tolerance_ppm = false;
        if (precursor_error_units == "ppm") // convert to Da (at m/z 666: 0.01 Da ~ 15 ppm)
        {
          search_parameters.precursor_mass_tolerance *= 2.0 / 3000.0;
          search_parameters.precursor_mass_tolerance_ppm = true;
        }

        search_parameters.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme));
        search_parameters.enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(tryptic_code);

        // create idXML file
        ProteinIdentification protein_id;
        protein_id.setPrimaryMSRunPath(primary_ms_run_path_);

        DateTime now = DateTime::now();
        std::string date_string = now.getDate();
        std::string identifier = "MS-GF+_" + date_string;

        protein_id.setIdentifier(identifier);
        protein_id.setDateTime(now);
        protein_id.setSearchParameters(search_parameters);
        protein_id.setSearchEngineVersion("");
        protein_id.setSearchEngine("MSGFPlus");
        protein_id.setScoreType(""); // MS-GF+ doesn't assign protein scores
    
        // store all peptide identifications in a map, the key is the scan number
        map<int, PeptideIdentification> peptide_identifications;
        set<std::string> prot_accessions;
    
        // iterate over the rows of the TSV file
        // columns: #SpecFile, SpecID, ScanNum, FragMethod, Precursor, IsotopeError, PrecursorError(ppm), Charge, Peptide, Protein, DeNovoScore, MSGFScore, SpecEValue, EValue, QValue, PepQValue
        // maybe TODO: replace column indexes ("elements[N]") by something more expressive
        CsvFile tsvfile(tsv_out, '\t');
        for (Size row_count = 1; row_count < tsvfile.rowCount(); ++row_count) // skip header line
        {
          vector<std::string> elements;
          if (!tsvfile.getRow(row_count, elements))
          {
            writeLogError_("Error: could not split row " + StringUtils::toStr(row_count) + " of file '" + tsv_out + "'");
            return PARSE_ERROR;
          }

          int scan_number = 0;
          if ((elements[2].empty()) || (elements[2] == "-1"))
          {
            scan_number = StringUtils::toInt32(elements[1]);
          }
          else
          {
            scan_number = StringUtils::toInt32(elements[2]);
          }

          struct SequenceParts parts = splitSequence_(elements[8]);
          StringUtils::substitute(parts.peptide, ',', '.'); // decimal separator should be dot, not comma
          AASequence seq = AASequence::fromString(modifySequence_(modifyNTermAASpecificSequence_(parts.peptide)));

          std::string accession = elements[9];
          // @BUG If there's a space before the protein accession in the FASTA file (e.g. "> accession ..."),
          // the "Protein" field in the TSV file will be empty, leading to an empty accession and no protein
          // reference in the idXML output file! (The mzIdentML output is not affected by this.)
          prot_accessions.insert(accession);

          PeptideEvidence evidence;
          evidence.setProteinAccession(accession);
          if ((parts.aa_before == 0) && (parts.aa_after == 0))
          {
            evidence.setAABefore(PeptideEvidence::UNKNOWN_AA);
            evidence.setAAAfter(PeptideEvidence::UNKNOWN_AA);
          }
          else // if one cleavage site is given, assume the other side is terminal
          {
            if (parts.aa_before != 0)
            {
              evidence.setAABefore(parts.aa_before);
            }
            else
            {
              evidence.setAABefore(PeptideEvidence::N_TERMINAL_AA);
            }
            if (parts.aa_after != 0)
            {
              evidence.setAAAfter(parts.aa_after);
            }
            else
            {
              evidence.setAAAfter(PeptideEvidence::C_TERMINAL_AA);
            }
          }

          bool hit_exists = false;
          // if the PeptideIdentification doesn't exist yet, a new one will be created:
          PeptideIdentification& pep_ident = peptide_identifications[scan_number];
          if (!pep_ident.getHits().empty()) // previously existing PeptideIdentification
          {
            // do we have a peptide hit with this sequence already?
            for (PeptideHit& hit : pep_ident.getHits())
            {
              if (hit.getSequence() == seq) // yes!
              {
                hit_exists = true;
                hit.addPeptideEvidence(evidence);
                break;
              }
            }
          }
          else // new PeptideIdentification
          {
            std::string spec_id = elements[1];
            pep_ident.setRT(rt_mapping[spec_id][0]);
            pep_ident.setMZ(rt_mapping[spec_id][1]);
            pep_ident.setMetaValue("ScanNumber", scan_number);
            pep_ident.setScoreType("SpecEValue");
            pep_ident.setHigherScoreBetter(false);
            pep_ident.setIdentifier(identifier);
          }
          if (!hit_exists) // add new PeptideHit
          {
            double score = StringUtils::toDouble(elements[12]);
            UInt rank = 0; // set to 0 at the moment
            Int charge = StringUtils::toInt32(elements[7]);
            PeptideHit hit(score, rank, charge, std::move(seq));
            hit.addPeptideEvidence(evidence);
            pep_ident.insertHit(hit);
          }
        }

        vector<ProteinHit> prot_hits;
        for (set<std::string>::iterator it = prot_accessions.begin(); it != prot_accessions.end(); ++it)
        {
          if (it->empty())
          {
            continue; // don't write a protein hit without accession (see @BUG above)
          }
          ProteinHit prot_hit = ProteinHit();
          prot_hit.setAccession(*it);
          prot_hits.push_back(prot_hit);
        }
        protein_id.setHits(prot_hits);
        protein_ids.push_back(protein_id);

        // iterate over map and create a vector of peptide identifications
        PeptideIdentification pep;
        for (map<int, PeptideIdentification>::iterator it = peptide_identifications.begin();
             it != peptide_identifications.end(); ++it)
        {
          pep = it->second;
          pep.sort();
          peptide_ids.push_back(pep);
        }
      }
      else // no legacy conversion
      {
        FileHandler().loadIdentifications(mzid_temp, protein_ids, peptide_ids, {FileTypes::MZIDENTML});

        // MzID might contain missed_cleavages set to -1 which leads to a crash in PeptideIndexer
        for (auto& pid : protein_ids)
        {
          pid.getSearchParameters().missed_cleavages = 1000; // use a high value (1000 was used in previous MSGF+ version)
          pid.getSearchParameters().digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme));
        }
        // set the MS-GF+ spectral e-value as new peptide identification score
        for (auto& pep : peptide_ids)
        { 
          switchScores_(pep);
        }

        // add missing RTs and FAIMS CVs to peptide IDs
        MSExperiment exp;
        MzMLFile mzml_file{};
        // Load spectrum metadata (not just file metadata) but skip peak data
        mzml_file.getOptions().setMetadataOnly(false);
        mzml_file.getOptions().setFillData(false);
        mzml_file.load(in, exp);
        SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(peptide_ids, exp);
        // Annotate FAIMS compensation voltage if present
        SpectrumMetaDataLookup::addMissingFAIMSToPeptideIDs(peptide_ids, exp);
      }

      // use OpenMS meta value key
      for (PeptideIdentification& pid : peptide_ids)
      {
        for (PeptideHit& psm : pid.getHits())
        {
          auto v = psm.getMetaValue("IsotopeError");
          // TODO cast to Int!
          psm.setMetaValue(Constants::UserParam::ISOTOPE_ERROR, v);
          psm.removeMetaValue("IsotopeError");
        }
      }

      // write all (!) parameters as metavalues to the search parameters
      if (!protein_ids.empty())
      {
        DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_ids[0].getSearchParameters(), this->getToolPrefix());
      }

      // if "reindex" parameter is set to true will perform reindexing
      if (auto ret = reindex_(protein_ids, peptide_ids); ret != EXECUTION_OK) return ret;

      // get feature set used in percolator
      StringList feature_set;
      PercolatorFeatureSetHelper::addMSGFFeatures(peptide_ids, feature_set);
      protein_ids.front().getSearchParameters().setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));

      FileHandler().storeIdentifications(out, protein_ids, peptide_ids, {FileTypes::IDXML, FileTypes::IDPARQUET});
    }

    //-------------------------------------------------------------
    // create (move) mzid output
    //-------------------------------------------------------------

    if (!mzid_out.empty())
    { // move the temporary file to the actual destination:
      if (!File::rename(mzid_temp, mzid_out))
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
    }

    return EXECUTION_OK;
  }


};


int main(int argc, const char** argv)
{
  MSGFPlusAdapter tool;

  return tool.main(argc, argv);
}

///@endcond
