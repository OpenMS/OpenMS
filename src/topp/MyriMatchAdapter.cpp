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
// $Maintainer:  Chris Bielow $
// $Authors: Dimitri Schachmann, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <iostream>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MyriMatchAdapter MyriMatchAdapter

  @brief Identifies peptides in MS/MS spectra via MyriMatch.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MyriMatchAdapter \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
    </tr>
  </table>
</CENTER>
  @em MyriMatch must be installed on the system to be able to use the @em MyriMatchAdapter.
  MyriMatch is currently available as part of the Bumbershoot package. See http://proteowizard.sourceforge.net/downloads.shtml.
  for further information on how to download and install @em MyriMatch on your system.

  This wrapper has been tested successfully with MyriMatch, version 2.1.x. and 2.2.x.

  Use debug level >=1 to keep intermediate pepXML and configuration files for manual inspection.

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MyriMatchAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_MyriMatchAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


using namespace OpenMS;
using namespace std;

class MyriMatchAdapter :
  public TOPPBase
{
public:
  MyriMatchAdapter() :
    TOPPBase("MyriMatchAdapter", "Annotates MS/MS spectra using MyriMatch.")
  {
  }

protected:

  struct MyriMatchVersion
  {
    MyriMatchVersion() :
      myrimatch_major(0), myrimatch_minor(0), myrimatch_patch(0)
    {}

    MyriMatchVersion(Int maj, Int min, Int pat) :
      myrimatch_major(maj), myrimatch_minor(min), myrimatch_patch(pat)
    {}

    Int myrimatch_major;
    Int myrimatch_minor;
    Int myrimatch_patch;

    bool operator<(const MyriMatchVersion& v) const
    {
      if (myrimatch_major > v.myrimatch_major) return false;
      else if (myrimatch_major < v.myrimatch_major) return true;
      else
      {
        if (myrimatch_minor > v.myrimatch_minor) return false;
        else if (myrimatch_minor < v.myrimatch_minor) return true;
        else
        {
          return myrimatch_patch < v.myrimatch_patch;
        }
      }
    }

  };

  bool getVersion_(const String& version, MyriMatchVersion& myrimatch_version_i) const
  {
    // we expect three components
    IntList nums = ListUtils::create<Int>(ListUtils::create<String>(version, '.'));
    if (nums.size() != 3) return false;

    myrimatch_version_i.myrimatch_major = nums[0];
    myrimatch_version_i.myrimatch_minor = nums[1];
    myrimatch_version_i.myrimatch_patch = nums[2];
    return true;
  }

  /// returns false on failure
  void translateModifications(StringList& static_mod_list, StringList& variable_mod_list)
  {
    // translating UNIMOD notation to MyriMatch notation of PTMs.
    ModificationDefinitionsSet mod_set(getStringList_("fixed_modifications"), getStringList_("variable_modifications"));
    if (!getStringList_("fixed_modifications").empty())
    {
      set<String> mod_names = mod_set.getFixedModificationNames();
      for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
      {
        ResidueModification mod = ModificationsDB::getInstance()->getModification(*it);
        String origin = String(mod.getOrigin());
        String mass_diff = String(mod.getDiffMonoMass());
        if (origin == "N-term")
        {
          origin = "(";
        }
        else if (origin == "C-term")
        {
          origin = ")";
        }
        else if (mod.getTermSpecificityName(mod.getTermSpecificity()) == "N-term")
        {
          origin = "(" + origin;
        }
        else if (mod.getTermSpecificityName(mod.getTermSpecificity()) == "C-term")
        {
          origin = ")" + origin;
        }
        static_mod_list.push_back(origin + " " + mod.getDiffMonoMass());
      }
    }

    if (!getStringList_("variable_modifications").empty())
    {
      set<String> mod_names = mod_set.getVariableModificationNames();

      for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
      {
        ResidueModification mod = ModificationsDB::getInstance()->getModification(*it);
        String origin = String(mod.getOrigin());
        String mass_diff = String(mod.getDiffMonoMass());
        if (origin == "N-term")
        {
          origin = "(";
        }
        else if (origin == "C-term")
        {
          origin = ")";
        }
        else if (mod.getTermSpecificityName(mod.getTermSpecificity()) == "N-term")
        {
          origin = "(" + origin;
        }
        else if (mod.getTermSpecificityName(mod.getTermSpecificity()) == "C-term")
        {
          origin = ")" + origin;
        }
        variable_mod_list.push_back(origin + " * " + mass_diff); // use * for all mods (no unique-per-mod symbol should be required)
      }
    }
  }

  void registerOptionsAndFlags_() override
  {
    addEmptyLine_();

    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 10.0, "Precursor monoisotopic mass tolerance.", false);

    registerStringOption_("precursor_mass_tolerance_unit", "<unit>", "ppm", "Unit to be used for precursor mass tolerance.", false);
    setValidStrings_("precursor_mass_tolerance_unit", ListUtils::create<String>("Da,ppm"));

    registerFlag_("precursor_mass_tolerance_avg", "If this flag is set, the average mass is used in the precursor mass tolerance.");
    registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.3, "Fragment mass error in Dalton", false);

    registerStringOption_("fragment_mass_tolerance_unit", "<unit>", "Da", "Unit to be used for fragment mass tolerance.", false);
    setValidStrings_("fragment_mass_tolerance_unit", ListUtils::create<String>("Da,ppm"));

    registerInputFile_("database", "<fasta-file>", "", "FASTA protein database.", true, false);
    setValidFormats_("database", ListUtils::create<String>("FASTA"));

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>(""),
                        "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>(""),
                        "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'.", false);
    setValidStrings_("variable_modifications", all_mods);

    addEmptyLine_();
    registerInputFile_("myrimatch_executable", "<executable>", "myrimatch",
                       "The 'myrimatch' executable of the MyriMatch installation", true, false, ListUtils::create<String>("skipexists"));
    registerIntOption_("NumChargeStates", "<num>", 3, "The number of charge states that MyriMatch will handle during all stages of the program.", false);
    registerDoubleOption_("TicCutoffPercentage", "<percentage>", 0.98, "Noise peaks are filtered out by sorting the original peaks in descending order of intensity, and then picking peaks from that list until the cumulative ion current of the picked peaks divided by the total ion current (TIC) is greater than or equal to this parameter.", false);

    registerIntOption_("MaxDynamicMods", "<num>", 2, "This parameter sets the maximum number of modified residues that may be in any candidate sequence.", false);
    registerIntOption_("MaxResultRank", "<rank>", 5, "This parameter sets the maximum rank of peptide-spectrum-matches to report for each spectrum.", false);
    registerStringOption_("CleavageRules", "<rule>", "", "This parameter allows the user to control the way peptides are generated from the protein database. For more details, see http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MS&termId=MS:1001045&termName=cleavage%20agent%20name .", false);
    setValidStrings_("CleavageRules", ListUtils::create<String>("Trypsin,Trypsin/P,Arg-C,Asp-N,Asp-N_ambic,CNBr,Chymotrypsin,Formic_acid,Lys-C,Lys-C/P,PepsinA,TrypChymo,V8-DE,V8-E,glutamyl endopeptidase,leukocyte elastase,no cleavage,proline endopeptidase,unspecific cleavage")); // NoEnzyme is deprecated according to http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MS&termId=MS:1001045&termName=cleavage%20agent%20name

    registerIntOption_("MinTerminiCleavages", "<num>", 2, "By default, when generating peptides from the protein database, a peptide must start and end at a valid cleavage site. Setting this parameter to 0 or 1 will reduce that requirement, so that neither terminus or only one terminus of the peptide must match one of the cleavage rules specified in the CleavageRules parameter. This parameter is useful to turn a tryptic digest into a semi-tryptic digest.", false); // TODO: Description copied from MM doc
    registerIntOption_("MaxMissedCleavages", "<num>", -1, "By default, when generating peptides from the protein database, a peptide may contain any number of missed cleavages. A missed cleavage is a site within the peptide that matches one of the cleavage rules (refer to CleavageRules). Settings this parameter to some other number will stop generating peptides from a sequence if it contains more than the specified number of missed cleavages.", false); // TODO: Description copied from MM doc

    // advanced options
    registerDoubleOption_("MinPeptideMass", "<mass>", 0.0, "When preprocessing the experimental spectra, any spectrum with a precursor mass that is less than the specified mass will be disqualified.", false, true); // Description copied from MM doc
    registerDoubleOption_("MaxPeptideMass", "<mass>", 10000.0, "When preprocessing the experimental spectra, any spectrum with a precursor mass that exceeds the specified mass will be disqualified.", false, true); // Description copied from MM doc
    registerIntOption_("MinPeptideLength", "<length>", 5, "When digesting proteins, any peptide which does not meet or exceed the specified length will be disqualified.", false, true); // TODO: Description copied from MM doc
    registerIntOption_("MaxPeptideLength", "<length>", 75, "When digesting proteins, any peptide which exceeds this specified length will be disqualified.", false, true); // TODO: Description copied from MM doc
    registerFlag_("UseSmartPlusThreeModel", "When this parameter is set, then for each peptide bond, an internal calculation is done to estimate the basicity of the b and y fragment sequence. The precursors protons are distributed to those ions based on that calculation, with the more basic sequence generally getting more of the protons..", true); // Description copied from MM doc
    registerIntOption_("NumIntensityClasses", "<num>", 3, "Before scoring any candidates, experimental spectra have their peaks stratified into the number of intensity classes specified by this parameter.", false, true); // Description copied from MM doc
    registerDoubleOption_("ClassSizeMultiplier", "<factor>", 2.0, "When stratifying peaks into a specified, fixed number of intensity classes, this parameter controls the size of each class relative to the class above it (where the peaks are more intense). ", false, true); // Description copied from MM doc
    registerStringOption_("MonoisotopeAdjustmentSet", "<set>", "0", "This parameter defines a set of isotopes (0 being the instrument-called monoisotope) to try as the monoisotopic precursor m/z. To disable this technique, set the value to '0'.", false, true);

    registerStringList_("SpectrumListFilters", "<filterList>", StringList(), "Optional set of filters as described in the MyriMatch documentation.", false, true);

    registerFlag_("ignoreConfigErrors", "Ignore wrong parameter names or values. Use with maximum caution!", true);

  }

  ExitCodes main_(int, const char**) override
  {
    String tmp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString()); // body for the tmp files
    {
      QDir d;
      d.mkpath(tmp_dir.toQString());
    }
    String logfile(getStringOption_("log"));
    String myrimatch_executable(getStringOption_("myrimatch_executable"));

    //-------------------------------------------------------------
    // get version of MyriMatch
    //-------------------------------------------------------------

    QProcess qp;
    String myrimatch_version;
    MyriMatchVersion myrimatch_version_i;

    // we invoke myrimatch w/o arguments. that yields a return code != 0. but
    // there is no other way for version 2.1 to get the version number
    qp.start(myrimatch_executable.toQString(), QStringList(), QIODevice::ReadOnly); // does automatic escaping etc...
    qp.waitForFinished();
    String output(QString(qp.readAllStandardOutput()));

    vector<String> lines;
    vector<String> version_split;
    output.split('\n', lines);

    // the version number is expected to be in the second line
    if (lines.size() < 2)
    {
      writeLog_("Warning: MyriMatch version output (" + output + ") not formatted as expected!");
      return EXTERNAL_PROGRAM_ERROR;
    }

    // the version is expected to be something like:
    // MyriMatch 2.1.111 (2011-12-27)
    lines[1].split(' ', version_split);
    if (version_split.size() == 3 && getVersion_(version_split[1], myrimatch_version_i))
    {
      myrimatch_version = version_split[1].removeWhitespaces();
      writeDebug_("Setting MyriMatch version to " + myrimatch_version, 1);
    }
    else
    {
      writeLog_("Warning: MyriMatch version output (" + output + ") not formatted as expected!");
      return EXTERNAL_PROGRAM_ERROR;
    }
    if (! (   (myrimatch_version_i.myrimatch_major == 2) &&  // major must be 2
              (myrimatch_version_i.myrimatch_minor == 1 || myrimatch_version_i.myrimatch_minor == 2) // minor .1 or .2
          ))
    {
      writeLog_("Warning: unsupported MyriMatch version (" + myrimatch_version + "). Tested only for MyriMatch 2.1.x and 2.2.x."
                "\nIf you encounter parameter errors, you can try the flag 'ignoreConfigErrors', but be aware that MyriMatch might be misconfigured.");
    }

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String inputfile_name = File::absolutePath(getStringOption_("in"));
    String outputfile_name = getStringOption_("out");
    String db_name = File::absolutePath(String(getStringOption_("database")));

    // building parameter String
    StringList parameters;

    if (getFlag_("ignoreConfigErrors")) parameters << "-ignoreConfigErrors";

    // Common Identification engine options
    StringList static_mod_list;
    StringList dynamic_mod_list;
    translateModifications(static_mod_list, dynamic_mod_list);
    if (!static_mod_list.empty())
      parameters << "-StaticMods" << ListUtils::concatenate(static_mod_list, " ");
    if (!dynamic_mod_list.empty())
      parameters << "-DynamicMods" << ListUtils::concatenate(dynamic_mod_list, " ");

    parameters << "-ProteinDatabase"  << File::absolutePath(db_name);

    if (getFlag_("precursor_mass_tolerance_avg"))
    {
      parameters << "-AvgPrecursorMzTolerance";
    }
    else
    {
      parameters << "-MonoPrecursorMzTolerance";
    }
    String precursor_mass_tolerance_unit = getStringOption_("precursor_mass_tolerance_unit") == "Da" ? " m/z" : " ppm";
    parameters << String(getDoubleOption_("precursor_mass_tolerance")) + precursor_mass_tolerance_unit;

    String fragment_mass_tolerance_unit = getStringOption_("fragment_mass_tolerance_unit");
    if (fragment_mass_tolerance_unit == "Da")
    {
      fragment_mass_tolerance_unit = "m/z";
    }

    parameters << "-FragmentMzTolerance" << String(getDoubleOption_("fragment_mass_tolerance")) + " " + fragment_mass_tolerance_unit;

    StringList slf = getStringList_("SpectrumListFilters");
    if (slf.size() > 0) 
    {
      if (myrimatch_version_i.myrimatch_minor <= 1)
      { // use quotes around the slf arguments (will be added automatically by Qt during call), i.e. "-SpectrumListFilters" "peakPicking false 2-"
        parameters << "-SpectrumListFilters" << ListUtils::concatenate(slf, ";") << "";
      }
      else
      { // no quotes -- pass a single argument, i.e. "-SpectrumListFilters peakPicking false 2-"
        parameters << "-SpectrumListFilters " + ListUtils::concatenate(slf, ";") << "";
      }
      
    }
    //parameters << "-ThreadCountMultiplier" << String(getIntOption_("threads")); // MyriMatch does not recognise this, even though it's in the manual.

    // MyriMatch specific parameters
    parameters << "-NumChargeStates" << getIntOption_("NumChargeStates");
    parameters << "-TicCutoffPercentage" << String(getDoubleOption_("TicCutoffPercentage"));
    parameters << "-MaxDynamicMods" << getIntOption_("MaxDynamicMods");
    parameters << "-MaxResultRank" << getIntOption_("MaxResultRank");
    parameters << "-MinTerminiCleavages" << getIntOption_("MinTerminiCleavages");
    parameters << "-MaxMissedCleavages" << getIntOption_("MaxMissedCleavages");
    String cleavage_rule = getStringOption_("CleavageRules");
    if (cleavage_rule.empty())
    {
      cleavage_rule = "Trypsin/P";
    }
    parameters << "-CleavageRules" << cleavage_rule;

    // advanced parameters
    parameters << "-MinPeptideMass"   << getDoubleOption_("MinPeptideMass");
    parameters << "-MaxPeptideMass"   << getDoubleOption_("MaxPeptideMass");
    parameters << "-MinPeptideLength" << getIntOption_("MinPeptideLength");
    parameters << "-MaxPeptideLength" << getIntOption_("MaxPeptideLength");
    parameters << "-NumIntensityClasses" << getIntOption_("NumIntensityClasses");
    parameters << "-ClassSizeMultiplier" << getDoubleOption_("ClassSizeMultiplier");
    parameters << "-MonoisotopeAdjustmentSet" << getStringOption_("MonoisotopeAdjustmentSet");
    parameters << "-cpus" << getIntOption_("threads");


    // Constant parameters

    // DecoyPrefix worked only when set through the config file
    String cfg_file = tmp_dir + "myrimatch.cfg";
    ofstream f(cfg_file.c_str());
    f << "DecoyPrefix=\"\"\n";
    f.close();
    parameters << "-cfg" << cfg_file;

    // path to input file must be the last parameter
    parameters << inputfile_name;

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    QStringList qparam;
    writeDebug_("MyriMatch arguments:", 1);
    writeDebug_(String("\"") + ListUtils::concatenate(parameters, "\" \"") + "\"", 1);
    for (Size i = 0; i < parameters.size(); ++i)
    {
      qparam << parameters[i].toQString();
    }
    QProcess process;

    // Bad style, because it breaks relative paths?
    process.setWorkingDirectory(tmp_dir.toQString());

    process.start(myrimatch_executable.toQString(), qparam, QIODevice::ReadOnly);
    bool success = process.waitForFinished(-1);
    String myri_msg(QString(process.readAllStandardOutput()));
    String myri_err(QString(process.readAllStandardError()));
    writeDebug_(myri_msg, 1);
    writeDebug_(myri_err, 0);
    if (!success || process.exitStatus() != 0 || process.exitCode() != 0)
    {
      writeLog_("Error: MyriMatch problem! (Details can be seen in the logfile: \"" + logfile + "\")");
      writeLog_("Note: This message can also be triggered if you run out of space in your tmp directory");
      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // reading MyriMatch output
    //-------------------------------------------------------------

    writeDebug_("Reading output of MyriMatch", 5);
    String exp_name = File::basename(inputfile_name);
    String pep_file = tmp_dir + File::removeExtension(exp_name) + ".pepXML";

    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> peptide_identifications;

    PeakMap exp;
    if (File::exists(pep_file))
    {
      MzMLFile fh;
      fh.load(inputfile_name, exp);

      SpectrumMetaDataLookup lookup;
      lookup.readSpectra(exp.getSpectra());
      PepXMLFile().load(pep_file, protein_identifications,
                        peptide_identifications, exp_name, lookup);
    }
    else
    {
      writeLog_("Error: MyriMatch problem! No pepXML output file (expected as '" + pep_file + "') was generated by MyriMatch.");
      writeLog_("Note: This message can be triggered if no MS2 spectra were found or no identifications were made.");
      writeLog_("      Myrimatch expects MS2 spectra in mzML files to contain the MSn tag. MSSpectrum with MS level 2 is not sufficient. You can use FileConverter to create such an mzML file by converting from mzML --> mzXML --> mzML.");
      return EXTERNAL_PROGRAM_ERROR;
    }

    if (debug_level_ == 0)
    {
      QFile(pep_file.toQString()).remove();
      QFile(cfg_file.toQString()).remove();
    }
    else
    {
      writeDebug_(String("Not removing '") + pep_file + "' for debugging purposes. Please delete manually!", 1);
      writeDebug_(String("Not removing '") + cfg_file + "' for debugging purposes. Please delete manually!", 1);
    }
    //-------------------------------------------------------------
    // writing results
    //-------------------------------------------------------------
    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = getStringOption_("database");
    ProteinIdentification::PeakMassType mass_type = getFlag_("precursor_mass_tolerance_avg") == true ? ProteinIdentification::AVERAGE : ProteinIdentification::MONOISOTOPIC;
    search_parameters.mass_type = mass_type;
    search_parameters.fixed_modifications = getStringList_("fixed_modifications");
    search_parameters.variable_modifications = getStringList_("variable_modifications");
    search_parameters.missed_cleavages = getIntOption_("MaxMissedCleavages");
    search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");
    search_parameters.precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
    search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor_mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.fragment_mass_tolerance_ppm = getStringOption_("fragment_mass_tolerance_unit") == "ppm" ? true : false;
    protein_identifications[0].setSearchParameters(search_parameters);
    protein_identifications[0].setSearchEngineVersion(myrimatch_version);
    protein_identifications[0].setSearchEngine("MyriMatch");

    if (!protein_identifications.empty())
    {
      StringList ms_runs;
      exp.getPrimaryMSRunPath(ms_runs);
      protein_identifications[0].setPrimaryMSRunPath(ms_runs);
    }
    IdXMLFile().store(outputfile_name, protein_identifications, peptide_identifications);
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  MyriMatchAdapter tool;
  return tool.main(argc, argv);
}

/// @endcond
