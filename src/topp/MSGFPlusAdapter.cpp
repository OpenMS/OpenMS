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
// $Maintainer: Hendrik Weisser $
// $Authors: Dilek Dere, Mathias Walzer, Petra Gutenbrunner, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/JavaInfo.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

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
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MSGFPlusAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes @n (or another centroiding tool)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    MS-GF+ must be installed before this wrapper can be used. Please make sure that Java and MS-GF+ are working.@n
    The following MS-GF+ version is required: MS-GF+ Beta (v10089) (7/31/2014). At the time of writing, it could be downloaded from http://omics.pnl.gov/software/ms-gf. Older versions will not work properly.

    Input spectra for MS-GF+ have to be centroided; profile spectra are ignored.

    The first time MS-GF+ is applied to a database (FASTA file), it will index the file contents and generate a number of auxiliary files in the same directory as the database (e.g. for "db.fasta": "db.canno", "db.cnlap", "db.csarr" and "db.cseq" will be generated). It is advisable to keep these files for future MS-GF+ searches, to save the indexing step.@n

    @note When a new database is used for the first time, make sure to run only one MS-GF+ search against it! Otherwise one process will start the indexing and the others will crash due to incomplete index files. After a database has been indexed, multiple MS-GF+ processes can use it in parallel.

    This adapter supports relative database filenames, which (when not found in the current working directory) are looked up in the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

    The adapter works in three steps to generate an idXML file: First MS-GF+ is run on the input MS data and the sequence database, producing an mzIdentML (.mzid) output file containing the search results. This file is then converted to a text file (.tsv) using MS-GF+' "MzIDToTsv" tool. Finally, the .tsv file is parsed and a result in idXML format is generated.

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
  public TOPPBase
{
public:
  MSGFPlusAdapter() :
    TOPPBase("MSGFPlusAdapter", "MS/MS database search using MS-GF+.", true),
    // parameter choices (the order of the values must be the same as in the MS-GF+ parameters!):
    fragment_methods_(ListUtils::create<String>("from_spectrum,CID,ETD,HCD")),
    instruments_(ListUtils::create<String>("low_res,high_res,TOF,Q_Exactive")),
    protocols_(ListUtils::create<String>("none,phospho,iTRAQ,iTRAQ_phospho,TMT")),
    tryptic_(ListUtils::create<String>("non,semi,fully"))
  {
    ProteaseDB::getInstance()->getAllMSGFNames(enzymes_);
    std::sort(enzymes_.begin(),enzymes_.end());
  }

protected:
  /// parts of a sequence of the form "K.AAAA.R"
  struct SequenceParts
  {
    char aa_before, aa_after; // may be '\0' if not given
    String peptide;

    SequenceParts(): aa_before(0), aa_after(0) {}
  };

  // lists of allowed parameter values:
  vector<String> fragment_methods_, instruments_, enzymes_, protocols_, tryptic_;

  // primary MS run referenced in the mzML file
  StringList primary_ms_run_path_;

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file (MS-GF+ parameter '-s')");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML,mgf,ms2"));
    registerOutputFile_("out", "<file>", "", "Output file", false);
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerOutputFile_("mzid_out", "<file>", "", "Alternative output file (MS-GF+ parameter '-o')\nEither 'out' or 'mzid_out' are required. They can be used together.", false);
    setValidFormats_("mzid_out", ListUtils::create<String>("mzid"));
    registerInputFile_("executable", "<file>", "MSGFPlus.jar", "MS-GF+ .jar file, e.g. 'c:\\program files\\MSGFPlus.jar'", true, false, ListUtils::create<String>("skipexists"));
    registerInputFile_("database", "<file>", "", "Protein sequence database (FASTA file; MS-GF+ parameter '-d'). Non-existing relative filenames are looked up via 'OpenMS.ini:id_db_dir'.", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));

    registerFlag_("add_decoys", "Create decoy proteins (reversed sequences) and append them to the database for the search (MS-GF+ parameter '-tda'). This allows the calculation of FDRs, but should only be used if the database does not already contain decoys.");

    registerDoubleOption_("precursor_mass_tolerance", "<value>", 10, "Precursor monoisotopic mass tolerance (MS-GF+ parameter '-t')", false);
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "Unit of precursor mass tolerance (MS-GF+ parameter '-t')", false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("Da,ppm"));

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

    registerFlag_("add_features", "Output additional features - needed e.g. by Percolator (default: basic scores only; MS-GF+ parameter '-addFeatures')", false);

    registerIntOption_("max_mods", "<num>", 2, "Maximum number of modifications per peptide. If this value is large, the search may take very long.", false);
    setMinInt_("max_mods", 0);

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", vector<String>(), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", vector<String>(), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    registerFlag_("legacy_conversion", "Use the indirect conversion of MS-GF+ results to idXML via export to TSV. Try this only if the default conversion takes too long or uses too much memory.", true);

    registerInputFile_("java_executable", "<file>", "java", "The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java", false, false, ListUtils::create<String>("skipexists"));
    registerIntOption_("java_memory", "<num>", 3500, "Maximum Java heap size (in MB)", false);
    registerIntOption_("java_permgen", "<num>", 0, "Maximum Java permanent generation space (in MB); only for Java 7 and below", false, true);
  }

  // The following sequence modification methods are used to modify the sequence stored in the TSV such that it can be used by AASequence

  // Method to cut the amino acids before/after the peptide (splice sites) off the sequence.
  // The sequences in the TSV file have the format 'K.XXXR.X' (where XXXR is the actual peptide sequence).
  // This method returns the sequence split into its three parts (e.g. "K", "XXXR", "X").
  struct SequenceParts splitSequence_(const String& sequence)
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
    parts.peptide = sequence.substr(start, count);
    return parts;
  }

  String modifyNTermAASpecificSequence_(const String& seq)
  {
    String swap;
    string modifiedSequence = seq;
    vector<pair<String, char> > massShiftList;

    massShiftList.push_back(make_pair("-18.011", 'E'));
    massShiftList.push_back(make_pair("-17.027", 'Q'));

    for (vector<pair<String, char> >::const_iterator it = massShiftList.begin(); it != massShiftList.end(); ++it)
    {
      string modMassShift = it->first;
      size_t found = modifiedSequence.find(modMassShift);

      if (found != string::npos)
      {
        String tmp = modifiedSequence.substr(0, found + modMassShift.length() + 1);
        size_t foundAA = tmp.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");

        if ((foundAA > found) && (tmp[foundAA] == it->second)) // no AA at the begin
        {
          if (found > 0)
          {
            swap = modifiedSequence.substr(0, found);
          }
          return swap += *tmp.rbegin() + modMassShift + modifiedSequence.substr(found + modMassShift.length() + 1);
        }
      }
    }
    return  modifiedSequence;
  }

  // Method to replace the mass representation of modifications.
  // Modifications in the TSV file have the format 'M+15.999'
  // After using this method the sequence should look like this: 'M[+15.999]'
  String modifySequence_(const String& seq)
  {
    String modifiedSequence = seq;
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
  void generateInputfileMapping_(map<String, vector<float> >& rt_mapping)
  {
    String exp_name = getStringOption_("in");

    if (!exp_name.empty())
    {
      PeakMap exp;
      // load only MS2 spectra:
      MzMLFile f;
      f.getOptions().addMSLevel(2);
      f.load(exp_name, exp);
      exp.getPrimaryMSRunPath(primary_ms_run_path_);

      if (exp.getSpectra().empty())
      {
        throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS2 spectra in input file.");
      }

      // determine type of spectral data (profile or centroided)
      SpectrumSettings::SpectrumType spectrum_type = exp[0].getType();

      if (spectrum_type == SpectrumSettings::RAWDATA)
      {
        if (!getFlag_("force"))
        {
          throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided MS2 spectra expected. To enforce processing of the data set the -force flag.");
        }
      }

      for (PeakMap::iterator it = exp.begin(); it != exp.end(); ++it)
      {
        String id = it->getNativeID(); // expected format: "... scan=#"
        if (id != "")
        {
          rt_mapping[id].push_back(it->getRT());
          rt_mapping[id].push_back(it->getPrecursors()[0].getMZ());
        }
      }
    }
  }

  String makeModString_(const String& mod_name, bool fixed=true)
  {
    ResidueModification mod = ModificationsDB::getInstance()->getModification(mod_name);
    char residue = mod.getOrigin();
    if (residue == 'X') residue = '*'; // terminal mod. without residue specificity
    String position = mod.getTermSpecificityName(); // "Prot-N-term", "Prot-C-term" not supported by OpenMS
    if (position == "none") position = "any";
    return String(mod.getDiffMonoMass()) + ", " + residue + (fixed ? ", fix, " : ", opt, ") + position + ", " + mod.getId() + "    # " + mod_name;
  }

  void writeModificationsFile_(const String& out_path, const vector<String>& fixed_mods, const vector<String>& variable_mods, Size max_mods)
  {
    ofstream output(out_path.c_str());
    output << "# MS-GF+ modifications file written by MSGFPlusAdapter (part of OpenMS)\n"
           << "NumMods=" << max_mods
           << "\n\n# Fixed modifications:\n";
    if (fixed_mods.empty())
    {
      output << "# (none)\n";
    }
    else
    {
      for (vector<String>::const_iterator it = fixed_mods.begin(); it != fixed_mods.end(); ++it)
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
      for (vector<String>::const_iterator it = variable_mods.begin(); it != variable_mods.end(); ++it)
      {
        output << makeModString_(*it, false) << "\n";
      }
    }
  }

  String describeHit_(const PeptideHit& hit)
  {
    return "peptide hit with sequence '" + hit.getSequence().toString() +
      "', charge " + String(hit.getCharge()) + ", score " +
      String(hit.getScore());
  }

  // Set the MS-GF+ e-value (MS:1002052) as new peptide identification score.
  void switchScores_(PeptideIdentification& id)
  {
    for (vector<PeptideHit>::iterator hit_it = id.getHits().begin(); hit_it != id.getHits().end(); ++hit_it)
    {
      // MS:1002052 == MS-GF spectral E-value
      if (!hit_it->metaValueExists("MS:1002052"))
      {
        String msg = "Meta value 'MS:1002052' not found for " + describeHit_(*hit_it);
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg);
      }

      hit_it->setScore(hit_it->getMetaValue("MS:1002052"));
    }
    id.setScoreType("SpecEValue");
    id.setHigherScoreBetter(false);
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parse parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String mzid_out = getStringOption_("mzid_out");
    if (mzid_out.empty() && out.empty())
    {
      writeLog_("Fatal error: no output file given (parameter 'out' or 'mzid_out')");
      return ILLEGAL_PARAMETERS;
    }

    String java_executable = getStringOption_("java_executable");
    String db_name = getStringOption_("database");
    if (!File::readable(db_name))
    {
      String full_db_name;
      try
      {
        full_db_name = File::findDatabase(db_name);
      }
      catch (...)
      {
        return ILLEGAL_PARAMETERS;
      }
      db_name = full_db_name;
    }

    vector<String> fixed_mods = getStringList_("fixed_modifications");
    vector<String> variable_mods = getStringList_("variable_modifications");
    bool no_mods = fixed_mods.empty() && variable_mods.empty();
    Int max_mods = getIntOption_("max_mods");
    if ((max_mods == 0) && !no_mods)
    {
      writeLog_("Warning: Modifications are defined ('fixed_modifications'/'variable_modifications'), but the number of allowed modifications is zero ('max_mods'). Is that intended?");
    }

    if (!getFlag_("force"))
    {
      if (!JavaInfo::canRun(java_executable))
      {
        writeLog_("Fatal error: Java is needed to run MS-GF+!");
        return EXTERNAL_PROGRAM_ERROR;
      }
    }
    else
    {
      writeLog_("The installation of Java was not checked.");
    }

    // create temporary directory (and modifications file, if necessary):
    String temp_dir, mzid_temp, mod_file;
    temp_dir = makeTempDirectory_();
    // always create a temporary mzid file first, even if mzid output is requested via "mzid_out"
    // (reason: TOPPAS may pass a filename with wrong extension to "mzid_out", which would cause an error in MzIDToTSVConverter below,
    // so we make sure that we have a properly named mzid file for the converter; see https://github.com/OpenMS/OpenMS/issues/1251)
    mzid_temp = temp_dir + "msgfplus_output.mzid";
    if (!no_mods)
    {
      mod_file = temp_dir + "msgfplus_mods.txt";
      writeModificationsFile_(mod_file, fixed_mods, variable_mods, max_mods);
    }

    // parameters also used by OpenMS (see idXML creation below):
    String enzyme = getStringOption_("enzyme");
    double precursor_mass_tol = getDoubleOption_("precursor_mass_tolerance");
    String precursor_error_units = getStringOption_("precursor_error_units");
    Int min_precursor_charge = getIntOption_("min_precursor_charge");
    Int max_precursor_charge = getIntOption_("max_precursor_charge");
    // parameters only needed for MS-GF+:
    QString java_memory = "-Xmx" + QString::number(getIntOption_("java_memory")) + "m";
    QString executable = getStringOption_("executable").toQString();
    // no need to handle "not found" case - would have given error during parameter parsing:
    Int fragment_method_code = ListUtils::getIndex<String>(fragment_methods_, getStringOption_("fragment_method"));
    Int instrument_code = ListUtils::getIndex<String>(instruments_, getStringOption_("instrument"));
    Int enzyme_code = ProteaseDB::getInstance()->getEnzyme(enzyme)->getMSGFID();
    Int protocol_code = ListUtils::getIndex<String>(protocols_, getStringOption_("protocol"));
    Int tryptic_code = ListUtils::getIndex<String>(tryptic_, getStringOption_("tryptic"));

    // Hack for KNIME. Looks for MSGFPLUS_PATH in the environment which is set in binaries.ini
    QProcessEnvironment env;
    String msgfpath = "MSGFPLUS_PATH";
    QString qmsgfpath = env.systemEnvironment().value(msgfpath.toQString());

    if (!qmsgfpath.isEmpty())
    {
      executable = qmsgfpath;
    }

    QStringList process_params; // the actual process is Java, not MS-GF+!
    process_params << java_memory
                   << "-jar" << executable
                   << "-s" << in.toQString()
                   << "-o" << mzid_temp.toQString()
                   << "-d" << db_name.toQString()
                   << "-t" << QString::number(precursor_mass_tol) + precursor_error_units.toQString()
                   << "-ti" << getStringOption_("isotope_error_range").toQString()
                   << "-tda" << QString::number(int(getFlag_("add_decoys")))
                   << "-m" << QString::number(fragment_method_code)
                   << "-inst" << QString::number(instrument_code)
                   << "-e" << QString::number(enzyme_code)
                   << "-protocol" << QString::number(protocol_code)
                   << "-ntt" << QString::number(tryptic_code)
                   << "-minLength" << QString::number(getIntOption_("min_peptide_length"))
                   << "-maxLength" << QString::number(getIntOption_("max_peptide_length"))
                   << "-minCharge" << QString::number(min_precursor_charge)
                   << "-maxCharge" << QString::number(max_precursor_charge)
                   << "-n" << QString::number(getIntOption_("matches_per_spec"))
                   << "-addFeatures" << QString::number(int(getFlag_("add_features")))
                   << "-thread" << QString::number(getIntOption_("threads"));

    if (!mod_file.empty())
    {
      process_params << "-mod" << mod_file.toQString();
    }

    //-------------------------------------------------------------
    // execute MS-GF+
    //-------------------------------------------------------------

    // run MS-GF+ process and create the .mzid file

    int status = QProcess::execute(java_executable.toQString(), process_params);
    if (status != 0)
    {
      writeLog_("Fatal error: Running MS-GF+ returned an error code '" + String(status) + "'. Does the MS-GF+ executable (.jar file) exist?");
      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // create idXML output
    //-------------------------------------------------------------

    if (!out.empty())
    {
      if (getFlag_("legacy_conversion"))
      {
        // run TSV converter
        String tsv_out = temp_dir + "msgfplus_converted.tsv";
        int java_permgen = getIntOption_("java_permgen");
        process_params.clear();
        process_params << java_memory;
        if (java_permgen > 0)
        {
          process_params << "-XX:MaxPermSize=" + QString::number(java_permgen) + "m";
        }
        process_params << "-cp" << executable << "edu.ucsd.msjava.ui.MzIDToTsv"
                       << "-i" << mzid_temp.toQString()
                       << "-o" << tsv_out.toQString()
                       << "-showQValue" << "1"
                       << "-showDecoy" << "1"
                       << "-unroll" << "1";
        status = QProcess::execute(java_executable.toQString(), process_params);
        if (status != 0)
        {
          writeLog_("Fatal error: Running MzIDToTSVConverter returned an error code '" + String(status) + "'.");
          return EXTERNAL_PROGRAM_ERROR;
        }

        // initialize map
        map<String, vector<float> > rt_mapping;
        generateInputfileMapping_(rt_mapping);

        // handle the search parameters
        ProteinIdentification::SearchParameters search_parameters;
        search_parameters.db = getStringOption_("database");
        search_parameters.charges = "+" + String(min_precursor_charge) + "-+" + String(max_precursor_charge);
        search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
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

        // create idXML file
        vector<ProteinIdentification> protein_ids;
        ProteinIdentification protein_id;
        protein_id.setPrimaryMSRunPath(primary_ms_run_path_);

        DateTime now = DateTime::now();
        String date_string = now.getDate();
        String identifier = "MS-GF+_" + date_string;

        protein_id.setIdentifier(identifier);
        protein_id.setDateTime(now);
        protein_id.setSearchParameters(search_parameters);
        protein_id.setSearchEngineVersion("");
        protein_id.setSearchEngine("MSGFPlus");
        protein_id.setScoreType(""); // MS-GF+ doesn't assign protein scores
    
        // store all peptide identifications in a map, the key is the scan number
        map<int, PeptideIdentification> peptide_identifications;
        set<String> prot_accessions;
    
        // iterate over the rows of the TSV file
        // columns: #SpecFile, SpecID, ScanNum, FragMethod, Precursor, IsotopeError, PrecursorError(ppm), Charge, Peptide, Protein, DeNovoScore, MSGFScore, SpecEValue, EValue, QValue, PepQValue
        // maybe TODO: replace column indexes ("elements[N]") by something more expressive
        CsvFile tsvfile(tsv_out, '\t');
        for (Size row_count = 1; row_count < tsvfile.rowCount(); ++row_count) // skip header line
        {
          vector<String> elements;
          if (!tsvfile.getRow(row_count, elements))
          {
            writeLog_("Error: could not split row " + String(row_count) + " of file '" + tsv_out + "'");
            return PARSE_ERROR;
          }

          int scan_number = 0;
          if ((elements[2] == "") || (elements[2] == "-1"))
          {
            scan_number = elements[1].suffix('=').toInt();
          }
          else
          {
            scan_number = elements[2].toInt();
          }

          struct SequenceParts parts = splitSequence_(elements[8]);
          parts.peptide.substitute(',', '.'); // decimal separator should be dot, not comma
          AASequence seq = AASequence::fromString(modifySequence_(modifyNTermAASpecificSequence_(parts.peptide)));

          String accession = elements[9];
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
            for (vector<PeptideHit>::iterator hit_it = pep_ident.getHits().begin();
                 hit_it != pep_ident.getHits().end(); ++hit_it)
            {
              if (hit_it->getSequence() == seq) // yes!
              {
                hit_exists = true;
                hit_it->addPeptideEvidence(evidence);
                break;
              }
            }
          }
          else // new PeptideIdentification
          {
            String spec_id = elements[1];
            pep_ident.setRT(rt_mapping[spec_id][0]);
            pep_ident.setMZ(rt_mapping[spec_id][1]);
            pep_ident.setMetaValue("ScanNumber", scan_number);
            pep_ident.setScoreType("SpecEValue");
            pep_ident.setHigherScoreBetter(false);
            pep_ident.setIdentifier(identifier);
          }
          if (!hit_exists) // add new PeptideHit
          {
            double score = elements[12].toDouble();
            UInt rank = 0; // set to 0 at the moment
            Int charge = elements[7].toInt();
            PeptideHit hit(score, rank, charge, seq);
            hit.addPeptideEvidence(evidence);
            pep_ident.insertHit(hit);
          }
        }

        vector<ProteinHit> prot_hits;
        for (set<String>::iterator it = prot_accessions.begin(); it != prot_accessions.end(); ++it)
        {
          if (it->empty()) continue; // don't write a protein hit without accession (see @BUG above)
          ProteinHit prot_hit = ProteinHit();
          prot_hit.setAccession(*it);
          prot_hits.push_back(prot_hit);
        }
        protein_id.setHits(prot_hits);
        protein_ids.push_back(protein_id);

        // iterate over map and create a vector of peptide identifications
        vector<PeptideIdentification> peptide_ids;
        PeptideIdentification pep;
        for (map<int, PeptideIdentification>::iterator it = peptide_identifications.begin();
             it != peptide_identifications.end(); ++it)
        {
          pep = it->second;
          pep.sort();
          peptide_ids.push_back(pep);
        }

        IdXMLFile().store(out, protein_ids, peptide_ids);
      }
      else
      {
        vector<ProteinIdentification> protein_ids;
        vector<PeptideIdentification> peptide_ids;
        MzIdentMLFile().load(mzid_temp, protein_ids, peptide_ids);
        // set the MS-GF+ spectral e-value as new peptide identification score
        for (vector<PeptideIdentification>::iterator pep_it = peptide_ids.begin(); pep_it != peptide_ids.end(); ++pep_it)
        {
          switchScores_(*pep_it);
        }

        SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(peptide_ids, in, false);

        IdXMLFile().store(out, protein_ids, peptide_ids);
      }
    }

    //-------------------------------------------------------------
    // create (move) mzid output
    //-------------------------------------------------------------

    if (!mzid_out.empty())
    {
      // existing file? Qt won't overwrite, so try to remove it:
      if (QFile::exists(mzid_out.toQString()) && !QFile::remove(mzid_out.toQString()))
      {
        writeLog_("Fatal error: Could not overwrite existing file '" + mzid_out + "'");
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      // move the temporary file to the actual destination:
      if (!QFile::rename(mzid_temp.toQString(), mzid_out.toQString()))
      {
        writeLog_("Fatal error: Could not move temporary mzid file to '" + mzid_out + "'");
        return CANNOT_WRITE_OUTPUT_FILE;
      }
    }

    removeTempDirectory_(temp_dir);

    return EXECUTION_OK;
  }


};


int main(int argc, const char** argv)
{
  MSGFPlusAdapter tool;

  return tool.main(argc, argv);
}
