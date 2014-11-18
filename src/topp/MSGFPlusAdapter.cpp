// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/JavaInfo.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

#include <fstream>
#include <map>
#include <cstddef>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_MSGFPlusAdapter MSGFPlusAdapter

   @brief 
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MSGFplusAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes @n (or another centroiding tool)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    MS-GF+ must be installed before this wrapper can be used. Please make sure that Java and MS-GF+ are working.

    Input spectra for MS-GF+ have to be centroided; profile spectra are ignored.

    This adapter supports relative database filenames, which (when not found in the 
    current working directory) are looked up in the directories specified by 
    'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).
		
    The adapter works in three steps: First MS-GF+ is run on the input MS data and the sequence database, producing an mzIdentML (.mzid) output file containing the search results. This file is then converted to a text file (.tsv) using MS-GF+' "MzIDToTsv" tool. Finally, the .tsv file is parsed and a result in idXML format is generated.

    This adapter has been tested mostly with the following MS-GF+ version: MS-GF+ Beta (v10089) (7/31/2014)

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_MSGFplusAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_MSGFplusAdapter.html
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
    TOPPBase("MSGFPlusAdapter", "MS/MS database search using MS-GF+.", false)
  {
    // parameter choices (the order of the values must be the same as in the MS-GF+ parameters!):
    fragment_methods_ = ListUtils::create<String>("from_spectrum,CID,ETD,HCD");
    instruments_ = ListUtils::create<String>("low_res,high_res,TOF,Q_Exactive");
    enzymes_ = ListUtils::create<String>("unspecific,trypsin,chymotrypsin,LysC,LysN,GluC,ArgC,AspN,alphaLP,no_cleavage");
    protocols_ = ListUtils::create<String>("none,phospho,iTRAQ,iTRAQ_phospho,TMT");
    tryptic_ = ListUtils::create<String>("non,semi,fully");
  }

protected:
  // lists of allowed parameter values:
  vector<String> fragment_methods_, instruments_, enzymes_, protocols_, tryptic_;

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file (MS-GF+ parameter '-s')");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML,mgf,ms2"));
    registerOutputFile_("out", "<file>", "", "Output file", false);
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerOutputFile_("mzid_out", "<file>", "", "Alternative output file (MS-GF+ parameter '-o')\nEither 'out' or 'mzid_out' are required. They can be used together.", false);
    setValidFormats_("mzid_out", ListUtils::create<String>("mzid"));
    registerInputFile_("executable", "<file>", "", "MS-GF+ .jar file, e.g. 'c:\\program files\\MSGFPlus.jar'");
    registerInputFile_("database", "<file>", "", "Protein sequence database (FASTA file; MS-GF+ parameter '-d'). Non-existing relative filenames are looked up via 'OpenMS.ini:id_db_dir'.", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));

    registerFlag_("add_decoys", "Create decoy proteins (reversed sequences) and append them to the database for the search (MS-GF+ parameter '-tda'). This allows the calculation of FDRs, but should only be used if the database does not already contain decoys.");

    registerDoubleOption_("precursor_mass_tolerance", "<value>", 20, "Precursor monoisotopic mass tolerance (MS-GF+ parameter '-t')", false);
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "Unit of precursor mass tolerance (MS-GF+ parameter '-t')", false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("Da,ppm"));

    registerStringOption_("isotope_error_range", "<range>", "0,1", "Range of allowed isotope peak errors (MS-GF+ parameter '-ti'). Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation. Combined with 'precursor_mass_tolerance'/'precursor_error_units', this determines the actual precursor mass tolerance. E.g. for experimental mass 'exp' and calculated mass 'calc', '-precursor_mass_tolerance 20 -precursor_error_units ppm -isotope_error_range -1,2' tests '|exp - calc - n * 1.00335 Da| < 20 ppm' for n = -1, 0, 1, 2.", false);

    registerStringOption_("fragment_method", "<choice>", fragment_methods_[0], "Fragmentation method ('from_spectrum' relies on spectrum meta data and uses CID as fallback option; MS-GF+ parameter '-m')", false);
    setValidStrings_("fragment_method", fragment_methods_);

    registerStringOption_("instrument", "<choice>", instruments_[0], "Instrument that generated the data ('low_res'/'high_res' refer to LCQ and LTQ instruments; MS-GF+ parameter '-inst')", false);
    setValidStrings_("instrument", instruments_);

    registerStringOption_("enzyme", "<choice>", enzymes_[1], "Enzyme used for digestion, or type of cleavage (MS-GF+ parameter '-e')", false);
    setValidStrings_("enzyme", enzymes_);

    registerStringOption_("protocol", "<choice>", protocols_[0], "Labeling or enrichment protocol used, if any (MS-GF+ parameter '-p')", false);
    setValidStrings_("protocol", protocols_);

    registerStringOption_("tryptic", "<choice>", tryptic_[2], "Level of cleavage specificity required (MS-GF+ parameter '-ntt')", false);
    setValidStrings_("tryptic", tryptic_);
    
    registerIntOption_("min_precursor_charge", "<num>", 2, "Minimum precursor ion charge (MS-GF+ parameter '-minCharge')", false);
    setMinInt_("min_precursor_charge", 1);
    registerIntOption_("max_precursor_charge", "<num>", 3, "Maximum precursor ion charge (MS-GF+ parameter '-maxCharge')", false);
    setMinInt_("max_precursor_charge", 1);

    registerIntOption_("min_peptide_length", "<num>", 6, "Minimum peptide length to consider (MS-GF+ parameter '-minLength')", false);
    setMinInt_("min_peptide_length", 1);
    registerIntOption_("max_peptide_length", "<num>", 40, "Maximum peptide length to consider (MS-GF+ parameter '-maxLength')", false);   
    setMinInt_("max_peptide_length", 1);

    registerIntOption_("matches_per_spec", "<num>", 1, "Number of matches per spectrum to be reported (MS-GF+ parameter '-n')", false);
    setMinInt_("matches_per_spec", 1);

    registerFlag_("add_features", "Output additional features (default: basic scores only; MS-GF+ parameter '-addFeatures')?", false);

    registerIntOption_("max_mods", "<num>", 2, "Maximum number of modifications per peptide. If this value is large, the search may take very long.", false);
    setMinInt_("max_mods", 0);

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    registerIntOption_("java_memory", "<num>", 3500, "Maximum Java heap size (in MB)", false);
    registerIntOption_("java_permgen", "<num>", 0, "Maximum Java permanent generation space (in MB); only for Java 7 and below", false, true);
  }

  // The following sequence modification methods are used to modify the sequence stored in the TSV such that it can be used by AASequence

  // Method to cut the amino acids before/after the peptide (splice sites) off the sequence.
  // The sequences in the TSV file have the format 'K.AAAA.R' (where AAAA is the actual peptide sequence).
  // After this method is used the sequence 'AAAA' results
  String cutSequence_(const String& sequence)
  {
    size_t len = sequence.size();
    if (len < 4) return sequence; // in 'X.Y', which side would we cut off?
    size_t start = 0, count = string::npos;
    if (sequence[1] == '.') start = 2;
    if (sequence[len - 2] == '.') count = len - start - 2;
    return sequence.substr(start, count);
  }

  String modifyNTermAASpecificSequence_(const String& seq)
  {
    String swap;
    string modifiedSequence = seq;
    vector<pair<String, char> > massShiftList;

    massShiftList.push_back(make_pair("-18.011", 'E'));
    massShiftList.push_back(make_pair("-17.027", 'Q'));

    for (vector<pair<String, char> >::const_iterator iter = massShiftList.begin(); iter != massShiftList.end(); iter++)
    {
      string modMassShift = iter->first;
      size_t found = modifiedSequence.find(modMassShift);

      if (found != string::npos)
      {
        String tmp = modifiedSequence.substr(0, found + modMassShift.length() + 1);
        size_t foundAA = tmp.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");

        if ((foundAA > found) && (tmp[foundAA] == iter->second)) // no AA at the begin
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
  void generateInputfileMapping_(Map<String, vector<float> >& rt_mapping)
  {
    String exp_name = getStringOption_("in");

    if (!exp_name.empty())
    {
      PeakMap exp;
      // load only MS2 spectra:
      MzMLFile f;
      f.getOptions().addMSLevel(2);
      f.load(exp_name, exp);

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
    String residue = mod.getOrigin();
    if (residue.size() != 1) residue = "*"; // specificity groups, e.g. "Deamidated (NQ)", are not supported by OpenMS
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

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String mzid_out = getStringOption_("mzid_out");
    if (mzid_out.empty() && out.empty())
    {
      writeLog_("Fatal error: no output file given (parameter 'out' or 'mzid_out')");
      return ILLEGAL_PARAMETERS;
    }

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

    if (!JavaInfo::canRun("java"))
    {
      writeLog_("Fatal error: Java not found. Java is needed to run MS-GF+. Make sure that it can be executed by calling 'java', e.g. add the directory containing the Java binary to your PATH variable.");
      return EXTERNAL_PROGRAM_ERROR;
    }

    // create temporary directory and modifications file, if necessary:
    String temp_dir, mod_file;
    if (!(out.empty() && no_mods)) // 'out' and 'mzid_out' can't both be empty
    {
      temp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString());
      writeDebug_("Creating temporary directory '" + temp_dir + "'", 1);
      QDir d;
      d.mkpath(temp_dir.toQString());
      if (mzid_out.empty())
      {
        mzid_out = temp_dir + "msgfplus_output.mzid";
      }
      if (!no_mods)
      {
        mod_file = temp_dir + "msgfplus_mods.txt";
        writeModificationsFile_(mod_file, fixed_mods, variable_mods, max_mods);
      }
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
    Int enzyme_code = ListUtils::getIndex<String>(enzymes_, enzyme);
    Int protocol_code = ListUtils::getIndex<String>(protocols_, getStringOption_("protocol"));
    Int tryptic_code = ListUtils::getIndex<String>(tryptic_, getStringOption_("tryptic"));      

    QStringList process_params; // the actual process is Java, not MS-GF+!
    process_params << java_memory
                   << "-jar" << executable
                   << "-s" << in.toQString()
                   << "-o" << mzid_out.toQString()
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
    int status = QProcess::execute("java", process_params);
    if (status != 0)
    {
      writeLog_("Fatal error: Running MS-GF+ returned an error code. Does the MS-GF+ executable (.jar file) exist?");
      return EXTERNAL_PROGRAM_ERROR;
    }

    if (out.empty()) 
    {
      removeTempDir_(temp_dir);
      return EXECUTION_OK; // no idXML required? -> we're finished now
    }

    //-------------------------------------------------------------
    // execute TSV converter
    //------------------------------------------------------------- 

    String tsv_out = temp_dir + "msgfplus_converted.tsv";
    int java_permgen = getIntOption_("java_permgen");
    process_params.clear();
    process_params << java_memory;
    if (java_permgen > 0) 
    {
      process_params << "-XX:MaxPermSize=" + QString::number(java_permgen) + "m";
    }
    process_params << "-cp" << executable << "edu.ucsd.msjava.ui.MzIDToTsv"
                   << "-i" << mzid_out.toQString()
                   << "-o" << tsv_out.toQString()
                   << "-showQValue" << "1"
                   << "-showDecoy" << "1"
                   << "-unroll" << "1";
    status = QProcess::execute("java", process_params);
    if (status != 0)
    {
      writeLog_("Fatal error: Running MzIDToTSVConverter returned an error code.");
      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // create idXML
    //------------------------------------------------------------- 

    // initialize map
    Map<String, vector<float> > rt_mapping;
    generateInputfileMapping_(rt_mapping);

    // handle the search parameters
    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = getStringOption_("database");
    search_parameters.charges = "+" + String(min_precursor_charge) + "-+" + String(max_precursor_charge);
    search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
    search_parameters.fixed_modifications = fixed_mods;
    search_parameters.variable_modifications = variable_mods;
    search_parameters.precursor_tolerance = precursor_mass_tol;
    if (precursor_error_units == "ppm") // convert to Da (at m/z 666: 0.01 Da ~ 15 ppm)
    {
      search_parameters.precursor_tolerance *= 2.0 / 3000.0;
    }

    ProteinIdentification::DigestionEnzyme enzyme_type = ProteinIdentification::UNKNOWN_ENZYME;
    if (enzyme == "trypsin")
    {
      enzyme_type = ProteinIdentification::TRYPSIN;
    } 
    else if (enzyme == "chymotrypsin")
    {
      enzyme_type = ProteinIdentification::CHYMOTRYPSIN;
    }
    else if (enzyme == "no_cleavage") 
    {
      enzyme_type = ProteinIdentification::NO_ENZYME ;     
    }
    search_parameters.enzyme = enzyme_type;

    // create idXML file
    vector<ProteinIdentification> protein_ids;
    ProteinIdentification protein_id;

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

    double score; // use SpecEValue from the TSV file
    UInt rank;
    Int charge;
    AASequence sequence;
    int scanNumber;

    // iterate over the rows of the TSV file
    CsvFile tsvfile(tsv_out, '\t');
    for (Size row_count = 1; row_count < tsvfile.rowCount(); ++row_count)
    {
      vector<String> elements;
      if (!tsvfile.getRow(row_count, elements))
      {
        writeLog_("Error: could not split row " + String(row_count) + " of file '" + tsv_out + "'");
        return PARSE_ERROR;
      }

      if ((elements[2] == "") || (elements[2] == "-1")) 
      {
        scanNumber = elements[1].suffix('=').toInt();
      } 
      else 
      {
        scanNumber = elements[2].toInt();
      }
      
      String seq = cutSequence_(elements[8]);
      seq.substitute(',', '.'); // decimal separator should be dot, not comma
      sequence = AASequence::fromString(modifySequence_(modifyNTermAASpecificSequence_(seq)));
      vector<PeptideHit> p_hits;
      String prot_accession = elements[9];

      if (prot_accessions.find(prot_accession) == prot_accessions.end()) 
      {
        prot_accessions.insert(prot_accession);
      }

      if (peptide_identifications.find(scanNumber) == peptide_identifications.end()) 
      {
        score = elements[12].toDouble();
        rank = 0; // set to 0 at the moment
        charge = elements[7].toInt();
        
        PeptideHit p_hit(score, rank, charge, sequence);
        p_hit.addProteinAccession(prot_accession);
        p_hits.push_back(p_hit);

        String spec_id = elements[1];
        peptide_identifications[scanNumber].setRT(rt_mapping[spec_id][0]);
        peptide_identifications[scanNumber].setMZ(rt_mapping[spec_id][1]);

        peptide_identifications[scanNumber].setMetaValue("ScanNumber", scanNumber);        
        peptide_identifications[scanNumber].setScoreType("SpecEValue");
        peptide_identifications[scanNumber].setHigherScoreBetter(false);
        peptide_identifications[scanNumber].setIdentifier(identifier);
      } 
      else 
      {
        p_hits = peptide_identifications[scanNumber].getHits();
        for (vector<PeptideHit>::iterator p_it = p_hits.begin(); p_it != p_hits.end(); ++ p_it)
        {
          if (p_it -> getSequence() == sequence) 
          {
            p_it -> addProteinAccession(prot_accession);
          }
        }
      }
      peptide_identifications[scanNumber].setHits(p_hits);
    }

    vector<ProteinHit> prot_hits;
    for (set<String>::iterator it = prot_accessions.begin(); it != prot_accessions.end(); ++ it) 
    {
      ProteinHit prot_hit = ProteinHit();
      prot_hit.setAccession(*it);
      prot_hits.push_back(prot_hit);
    }
    protein_id.setHits(prot_hits);
    protein_ids.push_back(protein_id);

    // iterate over map and create a vector of peptide identifications
    map<int, PeptideIdentification>::iterator it;
    vector<PeptideIdentification> peptide_ids;
    PeptideIdentification pep;
    for (map<int, PeptideIdentification>::iterator it = peptide_identifications.begin(); 
         it != peptide_identifications.end(); ++ it)
    {
      pep = it->second;
      pep.sort();
      peptide_ids.push_back(pep);
    }
    
    IdXMLFile().store(out, protein_ids, peptide_ids);
    
    removeTempDir_(temp_dir);

    return EXECUTION_OK;
  }	
};


int main(int argc, const char** argv)
{
  MSGFPlusAdapter tool;

  return tool.main(argc, argv);
}
