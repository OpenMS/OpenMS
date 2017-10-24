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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <iostream>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>
#include <qsignalmapper.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_OMSSAAdapter OMSSAAdapter

    @brief Identifies peptides in MS/MS spectra via OMSSA (Open Mass Spectrometry Search Algorithm).

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ OMSSAAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em OMSSA must be installed on the system to be able to use the @em OMSSAAdapter. See pubchem.ncbi.nlm.nih.gov/omssa/
    for further information on how to download and install @em OMSSA on your system. You might find that the latest OMSSA version
    does not run on your system (to test this, run @em omssacl in your OMMSA/bin/ directory and see if it crashes). If you encounter
    an error message, try another OMSSA version.

    @note OMMSA seems to be discontinued by NCBI due to financial restrictions. The above homepage might not work. Try ftp://ftp.ncbi.nih.gov/pub/lewisg/omssa/
    as an alternative for downloading the binaries.

    Sequence databases in FASTA format must be converted into the NCBI format before OMSSA can read them.
    For this, you can use the program 'formatdb' (old releases) or 'makeblastdb' (recent release) of the NCBI-tools suite, which is freely available for download.
    Use @em formatdb @em -i @em SwissProt_TargetAndDecoy.fasta @em -o to create a BLAST database, which actually consists of multiple files.
    The more recent 'makeblastdb' has a similar syntax, e.g., @em makeblastdb @em -dbtype @em prot @em -in @em SwissProt_TargetAndDecoy.fasta .

    Make sure that your FASTA file (which you convert into BLAST database files) is properly formatted, especially
    that it conforms to the FASTA Defline format and that there is a description(!), i.e. '>ID DESCRIPTION'. Otherwise you might get a

@code
  "..\..\..\..\..\..\..\src\algo\ms\omssa\omssacl.cpp", line 282: Fatal: COMSSA::Run() - Exception in COMSSA::Run: NCBI C++ Exception:
  "..\..\..\..\..\src\serial\serialobject.cpp", line 228: Error: NCBI-BlastDL::Blast-def-line.title
@endcode

    As database parameter for the OMSSAAdapter you can either specify the .psq file as generated by formatdb/makeblastdb
    (e.g., 'SwissProt_TargetAndDecoy.fasta.psq') or the original FASTA file (e.g., 'SwissProt_TargetAndDecoy.fasta').
    Allowing fasta format makes it easy to specify a common TOPPAS input node (using only the FASTA suffix) for multiple downstream adapters.
    Just make sure that the BLAST database files (.psq/.pin/.phr) and .fasta file reside in the same directory!

    This adapter supports relative database filenames, which (when not found in the current working directory) is looked up in
    the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

    The options that specify the protease specificity (@em e) are directly taken from OMSSA. A complete list of available
    proteases can be found by executing @em omssacl @em -el.

    Pre-build versions of OMSSA are 32bit for Windows and 64bit for Linux & MacOSX.
    If the input dataset contains many spectra (>30k), then a 32bit version of OMSSA will likely crash due to memory allocation issues.
    To prevent this, the adapter will automatically split the data into chunks of appropriate size (10k spectra by default) and call OMSSA for each chunk.
    Running time is about the same (slightly faster even) for 10k chunks, but deteriorates slightly (15%) if chunk size is too small (1k spectra).
    The disadvantage of chunking is that no protein hits (nor their scores) will be stored in the output, since peptide evidence is split between chunks.
    If you want to disable chunking at the risk of provoking a memory allocation error in OMSSA, set chunk size to '0'.

    This wrapper has been tested successfully with OMSSA, version 2.x.

    @note OMSSA search is much faster when the database (.psq files etc.) is accessed locally, rather than over a network share (we measured 10x speed increase in some cases).

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_OMSSAAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_OMSSAAdapter.html

    @improvement modes to read OMSSA output data and save in idXML format (Andreas)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPOMSSAAdapter :
  public TOPPBase
{
public:
  TOPPOMSSAAdapter() :
    TOPPBase("OMSSAAdapter", "Annotates MS/MS spectra using OMSSA.")
  {
  }

protected:

  struct OMSSAVersion
  {
    OMSSAVersion() :
      omssa_major(0), omssa_minor(0), omssa_patch(0)
    {}

    OMSSAVersion(Int maj, Int min, Int pat) :
      omssa_major(maj), omssa_minor(min), omssa_patch(pat)
    {}

    Int omssa_major;
    Int omssa_minor;
    Int omssa_patch;

    bool operator<(const OMSSAVersion& v) const
    {
      if (omssa_major > v.omssa_major) return false;
      else if (omssa_major < v.omssa_major) return true;
      else // ==
      {
        if (omssa_minor > v.omssa_minor) return false;
        else if (omssa_minor < v.omssa_minor) return true;
        else
        {
          return omssa_patch < v.omssa_patch;
        }
      }

    }

  };

  bool getVersion_(const String& version, OMSSAVersion& omssa_version_i) const
  {
    // we expect three components
    IntList nums = ListUtils::create<Int>(ListUtils::create<String>(version, '.'));
    if (nums.size() != 3) return false;

    omssa_version_i.omssa_major = nums[0];
    omssa_version_i.omssa_minor = nums[1];
    omssa_version_i.omssa_patch = nums[2];
    return true;
  }

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 10.0, "Precursor monoisotopic mass tolerance", false);
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "Unit of precursor mass tolerance", false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("Da,ppm"));
    registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.3, "Fragment mass error in Dalton", false);
    registerInputFile_("database", "<psq or fasta>", "", "NCBI formatted FASTA files. The .psq filename should be given, e.g. 'SwissProt.fasta.psq'. If the filename does not end in '.psq' (e.g., SwissProt.fasta) the psq suffix will be added automatically. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("psq,fasta"));
    registerIntOption_("min_precursor_charge", "<charge>", 1, "Minimum precursor ion charge", false);
    registerIntOption_("max_precursor_charge", "<charge>", 3, "Maximum precursor ion charge", false);
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    //Sequence library
    //-d <String> Blast sequence library to search.  Do not include .p* filename suffixes.
    //-pc <Integer> The number of pseudocounts to add to each precursor mass bin.
    //registerStringOption_("d", "<file>", "", "Blast sequence library to search.  Do not include .p* filename suffixes", true);
    registerInputFile_("omssa_executable", "<executable>", "omssacl", "The 'omssacl' executable of the OMSSA installation", true, false, ListUtils::create<String>("skipexists"));
    registerIntOption_("pc", "<Integer>", 1, "The number of pseudocounts to add to each precursor mass bin", false, true);

    //registerFlag_("omssa_out", "If this flag is set, the parameter 'in' is considered as an output file of OMSSA and will be converted to idXML");
    //registerStringOption_("omssa_out_format", "<type>", "", "Specifies the output format of OMSSA, if not given the format will be estimated", false);


    //Input format and filename
    //-f <String> single dta file to search
    //-fx <String> multiple xml-encapsulated dta files to search
    //-fb <String> multiple dta files separated by blank lines to search
    //-fm <String> mgf formatted file
    //-fp <String> pkl formatted file
    //-hs <Integer> the minimum number of m/z values a spectrum must have to be searched
    //-fxml <String> omssa xml search request file (contains search parameters and spectra. overrides command line)
    //-pm <String> search parameter input in xml format (contains search parameters but no spectra. overrides command line except for name of file containing spectra)
    // input options are not all necessary as TOPP tools only accept mzML
    registerIntOption_("hs", "<Integer>", 4, "the minimum number of m/z values a spectrum must have to be searched", false, true);
    //registerStringOption_("pm", "<file>", "", "search parameter input in xml format", false);

    //Output results
    //-o <String> filename for text asn.1 formatted search results
    //-ob <String> filename for binary asn.1 formatted search results
    //-ox <String> filename for xml formatted search results
    //-oc <String> filename for comma separated value (excel .csv) formatted search results
    // output options of OMSSA are not necessaryOMSSA

    //The following options output the search parameters and search spectra in the output results. This is necessary for viewing result in the OMSSA browser:
    //-w include spectra and search params in search results

    //To turn off informational messages (but not error messages), use:
    //-ni don't print informational messages

    //Mass type and tolerance
    //-to <float> product ion mass tolerance in Da
    //-te <float> precursor ion mass tolerance in Da
    //-tez <Integer> scaling of precursor mass tolerance with charge (0 = none, 1= linear)
    //registerDoubleOption_("to", "<float>", 0.8, "product ion mass tolerance in Da", false);
    //registerDoubleOption_("te", "<float>", 2.0, "precursor ion mass tolerance in Da", false);
    registerIntOption_("tez", "<Integer>", 1, "scaling of precursor mass tolerance with charge (0 = none, 1= linear)", false, true);

    //A precursor ion is the ion before fragmentation and the product ions are the ions generated after fragmentation. These values are specified in Daltons +/- the measured value, e.g. a value of 2.0 means +/- 2.0 Daltons of the measured value.
    //The tez value allows you to specify how the mass tolerance scales with the charge of the precursor. For example, you may search a precursor assuming that it has a charge state of 2+ and 3+. If you set tez to 1, then the mass tolerance for the +2 charge state will be 2 times the precursor mass tolerance, and for the 3+ charge state it will be 3 times the precursor mass tolerance. If you set tez to 0, the mass tolerance will always be equal to the precursor mass tolerance, irrespective of charge state.

    //-tom <Integer> product ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact.
    //-tem <Integer> precursor ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact.
    registerIntOption_("tom", "<Integer>", 0, "product ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact", false, true);
    registerIntOption_("tem", "<Integer>", 0, "precursor ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact", false, true);

    //Monoisotopic searching searches spectral peaks that correspond to peptides consisting entirely of carbon-12. Average mass searching searches on the average natural isotopic mass of peptides. Exact mass searches on the most abundant isotopic peak for a given mass range.

    //-tex <Double> threshold in Da above which the mass of a neutron should be added in an exact mass search.
    registerDoubleOption_("tex", "<float>", 1446.94, "threshold in Da above which the mass of a neutron should be added in an exact mass search", false, true);

    //Preprocessing
    //Preprocessing is the process of eliminating noise from a spectrum. Normally, you do not need to adjust these options as OMSSA automatically adjusts its preprocessing for best results.

    //-cl <float> low intensity cutoff as a fraction of max peak
    //-ch <float> high intensity cutoff as a fraction of max peak
    //-ci <float> intensity cutoff increment as a fraction of max peak
    //-w1 <Integer> single charge window in Da
    //-w2 <Integer> double charge window in Da
    //-h1 <Integer> number of peaks allowed in single charge window
    //-h2 <Integer> number of peaks allowed in double charge window
    //-cp <Integer> eliminate charge reduced precursors in spectra (0=no, 1=yes). Typically turned on for ETD spectra.

    //Charge Handling
    //Determination of precursor charge and product ion charges.  Presently, OMSSA estimates which precursors are 1+.  All other precursors are searched with charge from the minimum to maximum precursor charge specified.
    //-zl <Integer> minimum precursor charge to search when not 1+
    //-zh <Integer> maximum precursor charge to search when not 1+
    //-zt <Integer> minimum precursor charge to start considering multiply charged products
    //-z1 <Double> the fraction of peaks below the precursor used to determine if the spectrum is charge +1
    //-zc <Integer> should charge +1 be determined algorithmically (1=yes)
    //-zcc <Integer> how should precursor charges be determined? (1=believe the input file,2=use the specified range)
    //-zoh <Integer> set the maximum product charge to search

    //registerIntOption_("zl", "<Integer>", 1, "minimum precursor charge to search when not 1+", false);
    //registerIntOption_("zh", "<Integer>", 3, "maximum precursor charge to search when not 1+", false);
    registerIntOption_("zt", "<Integer>", 3, "minimum precursor charge to start considering multiply charged products", false, true);
    registerDoubleOption_("z1", "<float>", 0.95, "the fraction of peaks below the precursor used to determine if the spectrum is charge +1", false, true);
    registerIntOption_("zc", "<Integer>", 1, "should charge +1 be determined algorithmically (1=yes)", false, true);
    registerIntOption_("zcc", "<Integer>", 2, "how should precursor charges be determined? (1=believe the input file,2=use the specified range)", false, true);
    registerIntOption_("zoh", "<Integer>", 2, "set the maximum product charge to search", false, true);

    //Enzyme specification
    //Additional enzymes can be added upon request.
    //-v <Integer> number of missed cleavages allowed
    //-e <Integer> id number of enzyme to use (trypsin is the default)
    //-el print a list of enzymes and their corresponding id number
    //-no <Integer> minimum size of peptides for no-enzyme and semi-tryptic searches
    //-nox <Integer> maximum size of peptides for no-enzyme and semi-tryptic searches
    registerIntOption_("v", "<Integer>", 1, "number of missed cleavages allowed", false);
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllOMSSANames(all_enzymes);
    registerStringOption_("enzyme", "<enzyme>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("enzyme", all_enzymes);
    //registerIntOption_("e", "<Integer>", 0, "id number of enzyme to use (0 (i.e. trypsin) is the default, 17 would be no enzyme (i.e. unspecific digestion). Please refer to 'omssacl -el' for a listing.", false);
    registerIntOption_("no", "<Integer>", 4, "minimum size of peptides for no-enzyme and semi-tryptic searches", false, true);
    registerIntOption_("nox", "<Integer>", 40, "maximum size of peptides for no-enzyme and semi-tryptic searches", false, true);

    //Ions to search
    //OMSSA searches two ions series, both of which can be specified.  Normally one of the ion series specified is a forward ion series and the other is a reverse ion series.
    //-il print a list of ions and their corresponding id number
    //-i comma delimited list of id numbers of ions to search
    //-sp <Integer> number of product ions to search
    //-sb1 <Integer> should first forward (e.g. b1) product ions be searched (1 = no, 0 = yes)
    //-sct <Integer> should c terminus ions (e.g. y1) be searched (1 = no, 0 = yes)
    registerStringOption_("i", "<Num>,<Num>,<Num>", "1,4", "comma delimited list of id numbers of ions to search", false, true);
    registerIntOption_("sp", "<Integer>", 100, "number of product ions to search", false, true);
    registerIntOption_("sb1", "<Integer>", 1, "should first forward (e.g. b1) product ions be searched (1 = no, 0 = yes)", false, true);
    registerIntOption_("sct", "<Integer>", 0, "should c terminus ions (e.g. y1) be searched (1 = no, 0 = yes)", false, true);

    //Taxonomy
    //By default, OMSSA searches without limiting by taxonomy.  By specifying an NCBI taxonomy id, you can limit your search to a particular organism.  The taxonomy id can by found by searching the NCBI tax browser (enter the scientific name of the organism of interest in the search box and then click the correct search result and then the scientific name in the taxonomy browser to get the numeric taxonomy id).
    //-x comma delimited list of NCBI taxonomy ids to search (0 = all.  This is the default)
    registerStringOption_("x", "<Num>,<Num>,<Num>", "0", "comma delimited list of NCBI taxonomy ids to search (0 = all.  This is the default)", false, true);

    //Search heuristic parameters
    //These are options that can speed up the search.  They can result in decreased sensitivity
    //-hm <Integer> the minimum number of m/z matches a sequence library peptide must have for the hit to the peptide to be recorded
    //-ht <Integer> number of m/z values corresponding to the most intense peaks that must include one match to the theoretical peptide
    registerIntOption_("hm", "<Integer>", 2, "the minimum number of m/z matches a sequence library peptide must have for the hit to the peptide to be recorded", false, true);
    registerIntOption_("ht", "<Integer>", 6, "number of m/z values corresponding to the most intense peaks that must include one match to the theoretical peptide", false, true);

    //Results
    //-hl <Integer> maximum number of hits retained for one spectrum
    //-he <Double> the maximum e-value allowed in the hit list
    registerIntOption_("hl", "<Integer>", 30, "maximum number of hits retained for one spectrum. Note: even when set to 1 OMSSA may report multiple hits with different charge states", false);
    registerDoubleOption_("he", "<float>", 1000, "the maximum e-value allowed in the hit list. If you set this parameter too small (e.g., he=1), this will effectively introduce FDR filtering."
                                                 " Thus, allowing a less stringent FDR during post-processing will nevertheless return the (better) FDR introduced here, since mediocre hits are not even reported.", false);

    //Post translational modifications
    //To specify modifications, first type in "omssacl -ml" to see a list of modifications available and their corresponding id number.  Then when running the search, specify the id numbers of the modification you wish to apply, e.g. "omssacl -mf 5 -mv 1,8 ...". Multiple PTMs can be specified by placing commas between the numbers without any spaces.  At the present time, the list of allowed post translational modifications will be expanded over time.
    //-mf  comma delimited list of id numbers for fixed modifications
    //-mv  comma delimited list of id numbers for variable modifications
    //-ml  print a list of modifications and their corresponding id number
    //registerStringOption_("mf", "<Num>,<Num>,<Num>", "", "comma delimited list of id numbers for fixed modifications", false);
    //registerStringOption_("mv", "<Num>,<Num>,<Num>", "", "comma delimited list of id numbers for variable modifications", false);
    // TODO -ml
    //registerStringOption_("mux", "<file>", "", "use the given file which contains user modifications in OMSSA modifications xml format", false);


    //To add your own user defined modifications, edit the usermod0-29 entries in the mods.xml file. If it is common modification, please contact NCBI so that it can be added to the standard list.
    //To reduce the combinatorial expansion that results when specifying multiple variable modifications, you can put an upper bound on the number of mass ladders generated per peptide using the -mm option.  The ladders are generated in the order of the least number of modification to the most number of modifications.
    //-mm <Integer> the maximum number of mass ladders to generate per database peptide
    registerIntOption_("mm", "<Integer>", 128, "the maximum number of mass ladders to generate per database peptide", false, true);

    //There is an upper bound on the number of combinations of variable mods that can be applied to a peptide from the sequence library. The hard upper bound is 1024, which effectively limits the number of variable modification sites per peptide for an exhaustive search to 10. If you set this number too low, you will miss highly modified peptides. If you set it too high, it will make the e-values less significant by searching for too many possible modifications.
    //OMSSA treats cleavage of the initial methionine in each protein record as a variable modification by default. To turn off this behavior use the command line option
    //-mnm n-term methionine should not be cleaved
    registerFlag_("mnm", "n-term methionine should not be cleaved", true);

    //Iterative searching
    //-is <Double> evalue threshold to include a sequence in the iterative search, 0 = all
    //-ir <Double> evalue threshold to replace a hit, 0 = only if better
    //-ii <Double> evalue threshold to iteratively search a spectrum again, 0 = always
    registerDoubleOption_("is", "<float>", 0.0, "evalue threshold to include a sequence in the iterative search, 0 = all", false, true);
    registerDoubleOption_("ir", "<float>", 0.0, "evalue threshold to replace a hit, 0 = only if better", false, true);
    registerDoubleOption_("ii", "<float>", 0.0, "evalue threshold to iteratively search a spectrum again, 0 = always", false, true);


    //-foms <String> read in search result in .oms format (binary asn.1).
    //-fomx <Double> read in search result in .omx format (xml).
    //Iterative searching is the ability to re-search search results in hopes of increasing the number of spectra identified. To accomplish this, an iterative search may change search parameters, such as using a no-enzyme search, or restrict the sequence search library to sequences already hit.

    registerIntOption_("chunk_size", "<Integer>", 0, "Number of spectra to submit in one chunk to OMSSA. Chunks with more than 30k spectra will likely cause memory allocation issues with 32bit OMSSA versions (which is usually the case on Windows). To disable chunking (i.e. submit all spectra in one big chunk), set it to '0'.", false, true);
  }

  ExitCodes main_(int, const char**)
  {
    StringList parameters;
    // path to the log file
    String logfile(getStringOption_("log"));
    String omssa_executable(getStringOption_("omssa_executable"));
    String unique_name = QDir::toNativeSeparators(String(File::getTempDirectory() + "/" + File::getUniqueName()).toQString()); // body for the tmp files
    String unique_input_name = unique_name + "_OMSSA"; // mfg
    String unique_output_name = unique_name + "_OMSSA"; // xml (OMSSA)
    String unique_usermod_name = unique_name + "_OMSSA_user_mod_file.xml";

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    // get version of OMSSA
    QProcess qp;
    qp.start(omssa_executable.toQString(), QStringList() << "-version", QIODevice::ReadOnly); // does automatic escaping etc...
    bool success = qp.waitForFinished();
    String output(QString(qp.readAllStandardOutput()));
    String omssa_version;
    OMSSAVersion omssa_version_i;
    if (!success || qp.exitStatus() != 0 || qp.exitCode() != 0)
    {
      writeLog_("Warning: unable to determine the version of OMSSA - the process returned an error. Call string was: '" + omssa_executable + " -version'. Make sure that OMSSA exists and the path given in '-omssa_executable' is correct!");
      return ILLEGAL_PARAMETERS;
    }
    else
    {
      vector<String> version_split;
      output.split(' ', version_split);
      if (version_split.size() == 2 && getVersion_(version_split[1], omssa_version_i))
      {
        omssa_version = version_split[1].removeWhitespaces();
        writeDebug_("Setting OMSSA version to " + omssa_version, 1);
      }
      else
      {
        writeLog_("Warning: OMSSA version output (" + output + ") not formatted as expected!");
      }
    }
    // parse arguments
    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");
    String db_name = String(getStringOption_("database"));
    // @todo: find DB for OMSSA (if not given) in OpenMS_bin/share/OpenMS/DB/*.fasta|.pin|...

    //-------------------------------------------------------------
    // Validate user parameters
    //-------------------------------------------------------------
    if (getIntOption_("min_precursor_charge") > getIntOption_("max_precursor_charge"))
    {
      LOG_ERROR << "Given charge range is invalid: max_precursor_charge needs to be >= min_precursor_charge." << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    if (db_name.suffix('.') != "psq")
    {
      db_name += ".psq";
    }

    if (!File::readable(db_name))
    {
      String full_db_name;
      try
      {
        full_db_name = File::findDatabase(db_name);
      }
      catch (...)
      {
        LOG_ERROR << "Unable to find database '" << db_name << "' (searched all folders). Did you mistype its name?" << std::endl;
        return ILLEGAL_PARAMETERS;
      }
      db_name = full_db_name;
    }

    db_name = db_name.substr(0, db_name.size() - 4); // OMSSA requires the filename without the .psq part
    // check for .pin and .phr files
    bool has_pin = File::readable(db_name + ".pin");
    bool has_phr = File::readable(db_name + ".phr");
    if (!has_pin || !has_phr)
    {
      LOG_ERROR << "\nThe NCBI psq database '" << db_name << ".psq' was found, but the following associated index file(s) are missing:\n";
      if (!has_pin) LOG_ERROR << "  missing: '" << db_name << ".pin'\n";
      if (!has_phr) LOG_ERROR << "  missing: '" << db_name << ".phr'\n";
      LOG_ERROR << "Please make sure the file(s) are present!\n" << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    bool db_name_contains_space = false;
    if (db_name.hasSubstring(" "))
    {
      db_name_contains_space = true;
    }
    // This is a workaround for a bug in the NCBI libraries.
    // They internally don't support spaces in path or file names.
    if (db_name_contains_space)
    {
#ifdef OPENMS_WINDOWSPLATFORM
      // Windows: use doubly escaped double quotes (and do a system call instead of QProcess later)
      parameters << "-d" << String("\"\\\"") + String(db_name) + String("\\\"\"");
#else
      // Linux/Mac: wrap into singly escaped double quotes
      parameters << "-d" << String("\"") + String(db_name) + String("\"");
#endif
    }
    else
    {
      parameters << "-d" << String(db_name);
    }

    parameters << "-to" << String(getDoubleOption_("fragment_mass_tolerance")); //String(getDoubleOption_("to"));
    parameters << "-hs" << String(getIntOption_("hs"));
    parameters << "-te" << String(getDoubleOption_("precursor_mass_tolerance")); //String(getDoubleOption_("te"));
    if (getStringOption_("precursor_error_units") == "ppm")
    {
      if (omssa_version_i < OMSSAVersion(2, 1, 8))
      {
        writeLog_("This OMSSA version (" + omssa_version + ") does not support ppm tolerances."
                  + " Please disable it and set the precursor tolerance in Da."
                  + " Required version is 2.1.8 and above.\n");
        return ILLEGAL_PARAMETERS;
      }
      parameters << "-teppm"; // only from OMSSA 2.1.8 on
    }

    parameters << "-zl" << String(getIntOption_("min_precursor_charge")); //String(getIntOption_("zl"));
    parameters << "-zh" <<  String(getIntOption_("max_precursor_charge")); //String(getIntOption_("zh"));
    parameters << "-zt" <<  String(getIntOption_("zt"));
    parameters << "-zc" <<  String(getIntOption_("zc"));
    parameters << "-zcc" << String(getIntOption_("zcc"));
    parameters << "-zoh" << String(getIntOption_("zoh"));
    parameters << "-no" << String(getIntOption_("no"));
    parameters << "-nox" << String(getIntOption_("nox"));
    parameters << "-sp" << String(getIntOption_("sp"));
    parameters << "-sb1" << String(getIntOption_("sb1"));
    parameters << "-sct" << String(getIntOption_("sct"));
    parameters << "-x" << getStringOption_("x");
    parameters << "-hl" << String(getIntOption_("hl"));
    parameters << "-hm" << String(getIntOption_("hm"));
    parameters << "-ht" <<  String(getIntOption_("ht"));
    parameters << "-tex" << String(getDoubleOption_("tex"));
    parameters << "-i" << getStringOption_("i");
    parameters << "-z1" << String(getDoubleOption_("z1"));
    parameters << "-v" << String(getIntOption_("v"));
    parameters << "-e" << String(ProteaseDB::getInstance()->getEnzyme(getStringOption_("enzyme"))->getOMSSAID());
    parameters << "-tez" << String(getIntOption_("tez"));


    parameters << "-tom" << String(getIntOption_("tom"));
    parameters << "-tem" << String(getIntOption_("tem"));

    parameters << "-mm" << String(getIntOption_("mm"));
    parameters << "-is" << String(getDoubleOption_("is"));
    parameters << "-ir" << String(getDoubleOption_("ir"));
    parameters << "-ii" << String(getDoubleOption_("ii"));
    parameters << "-nt" << String(getIntOption_("threads"));

    if (getFlag_("mnm"))
    {
      parameters << "-mnm";
    }

    if (getIntOption_("debug") == 0)
    {
      parameters << "-ni";
    }
    parameters << "-he" << String(getDoubleOption_("he"));


    // read mapping for the modifications
    String file = File::find("CHEMISTRY/OMSSA_modification_mapping");

    TextFile infile(file);
    Map<String, UInt> mods_map;
    for (TextFile::ConstIterator it = infile.begin(); it != infile.end(); ++it)
    {
      vector<String> split;
      it->split(',', split);

      if (it->size() > 0 && (*it)[0] != '#')
      {
        if (split.size() < 2)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "parse mapping file line: '" + *it + "'", "");
        }
        for (Size i = 2; i != split.size(); ++i)
        {
          String tmp(split[i].trim());
          if (!tmp.empty())
          {
            mods_map[tmp] = split[0].trim().toInt();
          }
        }
      }
    }

    writeDebug_("Evaluating modifications", 1);
    ModificationDefinitionsSet mod_set(getStringList_("fixed_modifications"), getStringList_("variable_modifications"));
    writeDebug_("Setting modifications", 1);
    UInt user_mod_num(119);
    vector<pair<UInt, String> > user_mods;
    // fixed modifications
    if (!getStringList_("fixed_modifications").empty())
    {
      set<String> mod_names = mod_set.getFixedModificationNames();
      StringList mod_list;
      for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
      {
        if (mods_map.has(*it))
        {
          mod_list.push_back(String(mods_map[*it]));
        }
        else
        {
          mod_list.push_back(String(user_mod_num));
          // add this to the usermods
          user_mods.push_back(make_pair(user_mod_num++, *it));
          writeDebug_("Inserting unknown fixed modification: '" + *it + "' into OMSSA", 1);
        }
      }
      if (mod_list.size() > 0)
      {
        parameters << "-mf" << ListUtils::concatenate(mod_list, ",");
      }
    }

    if (!getStringList_("variable_modifications").empty())
    {
      set<String> mod_names = mod_set.getVariableModificationNames();
      StringList mod_list;

      for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
      {
        if (mods_map.has(*it))
        {
          mod_list.push_back(String(mods_map[*it]));
        }
        else
        {
          mod_list.push_back(String(user_mod_num));
          // add this to the usermods
          user_mods.push_back(make_pair(user_mod_num++, *it));
          writeDebug_("Inserting unknown variable modification: '" + *it + "' into OMSSA", 1);
        }
      }

      if (mod_list.size() > 0)
      {
        parameters << "-mv" << ListUtils::concatenate(mod_list, ",");
      }
    }

    // write unknown modifications to user mods file
    if (!user_mods.empty())
    {
      writeDebug_("Writing usermod file to " + unique_usermod_name, 1);
      parameters << "-mux" << File::absolutePath(unique_usermod_name);
      ofstream out(unique_usermod_name.c_str());
      out << "<?xml version=\"1.0\"?>" << "\n";
      out << "<MSModSpecSet xmlns=\"http://www.ncbi.nlm.nih.gov\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\" xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov OMSSA.xsd\">" << "\n";

      UInt user_mod_count(1);
      for (vector<pair<UInt, String> >::const_iterator it = user_mods.begin(); it != user_mods.end(); ++it)
      {
        writeDebug_("Writing information into user mod file of modification: " + it->second, 1);
        out << "<MSModSpec>" << "\n";
        out << "\t<MSModSpec_mod>" << "\n";
        out << "\t\t<MSMod value=\"usermod" << user_mod_count++ << "\">" << it->first << "</MSMod>" << "\n";
        out << "\t</MSModSpec_mod>" << "\n";
        out << "\t<MSModSpec_type>" << "\n";

        /*
            0 modaa	-  at particular amino acids
            1 modn	-  at the N terminus of a protein
            2 modnaa	-  at the N terminus of a protein at particular amino acids
            3 modc	-  at the C terminus of a protein
            4 modcaa	-  at the C terminus of a protein at particular amino acids
            5 modnp	-  at the N terminus of a peptide
            6 modnpaa	-  at the N terminus of a peptide at particular amino acids
            7 modcp	-  at the C terminus of a peptide
            8 modcpaa	-  at the C terminus of a peptide at particular amino acids
            9 modmax	-  the max number of modification types
        */

        ResidueModification::TermSpecificity ts = ModificationsDB::getInstance()->getModification(it->second).getTermSpecificity();
        String origin = ModificationsDB::getInstance()->getModification(it->second).getOrigin();
        if (ts == ResidueModification::ANYWHERE)
        {
          out << "\t\t<MSModType value=\"modaa\">0</MSModType>" << "\n";
        }
        if (ts == ResidueModification::C_TERM)
        {
          if (origin == "" || origin == "X")
          {
            out << "\t\t<MSModType value=\"modcp\">7</MSModType>" << "\n";
          }
          else
          {
            out << "\t\t<MSModType value=\"modcpaa\">8</MSModType>" << "\n";
          }
        }
        if (ts == ResidueModification::N_TERM)
        {
          if (origin == "" || origin == "X")
          {
            out << "\t\t<MSModType value=\"modnp\">5</MSModType>" << "\n";
          }
          else
          {
            out << "\t\t<MSModType value=\"modnpaa\">6</MSModType>" << "\n";
          }
        }
        out << "\t</MSModSpec_type>" << "\n";

        out << "\t<MSModSpec_name>" << it->second << "</MSModSpec_name>" << "\n";
        out << "\t<MSModSpec_monomass>" << ModificationsDB::getInstance()->getModification(it->second).getDiffMonoMass()  << "</MSModSpec_monomass>" << "\n";
        out << "\t<MSModSpec_averagemass>" << ModificationsDB::getInstance()->getModification(it->second).getDiffAverageMass() << "</MSModSpec_averagemass>" << "\n";
        out << "\t<MSModSpec_n15mass>0</MSModSpec_n15mass>" << "\n";

        if (origin != "")
        {
          out << "\t<MSModSpec_residues>" << "\n";
          out << "\t\t<MSModSpec_residues_E>" << origin << "</MSModSpec_residues_E>" << "\n";
          out << "\t</MSModSpec_residues>" << "\n";

          /* TODO: Check why these are always 0
          double neutral_loss_mono = ModificationsDB::getInstance()->getModification(it->second).getNeutralLossMonoMass();
          double neutral_loss_avg = ModificationsDB::getInstance()->getModification(it->second).getNeutralLossAverageMass();
          */
          double neutral_loss_mono = ModificationsDB::getInstance()->getModification(it->second).getNeutralLossDiffFormula().getMonoWeight();
          double neutral_loss_avg = ModificationsDB::getInstance()->getModification(it->second).getNeutralLossDiffFormula().getAverageWeight();

          if (fabs(neutral_loss_mono) > 0.00001)
          {
            out << "\t<MSModSpec_neutralloss>" << "\n";
            out << "\t\t<MSMassSet>" << "\n";
            out << "\t\t\t<MSMassSet_monomass>" << neutral_loss_mono << "</MSMassSet_monomass>" << "\n";
            out << "\t\t\t<MSMassSet_averagemass>" << neutral_loss_avg << "</MSMassSet_averagemass>" << "\n";
            out << "\t\t\t<MSMassSet_n15mass>0</MSMassSet_n15mass>" << "\n";
            out << "\t\t</MSMassSet>" << "\n";
            out << "\t</MSModSpec_neutralloss>" << "\n";
          }

          out << "</MSModSpec>" << "\n";
        }
      }

      out << "</MSModSpecSet>" << "\n";
      out.close();
    }

    // prepare some datastructures for result annotation
    // OMSSA does not write fixed modifications so we need to add them to the sequences
    set<String> fixed_mod_names = mod_set.getFixedModificationNames();
    vector<String> fixed_nterm_mods, fixed_cterm_mods;
    Map<String, String> fixed_residue_mods;
    writeDebug_("Splitting modification into N-Term, C-Term and anywhere specificity", 1);
    for (set<String>::const_iterator it = fixed_mod_names.begin(); it != fixed_mod_names.end(); ++it)
    {
      ResidueModification::TermSpecificity ts = ModificationsDB::getInstance()->getModification(*it).getTermSpecificity();
      if (ts == ResidueModification::ANYWHERE)
      {
        fixed_residue_mods[ModificationsDB::getInstance()->getModification(*it).getOrigin()] = *it;
      }
      if (ts == ResidueModification::C_TERM)
      {
        fixed_cterm_mods.push_back(*it);
      }
      if (ts == ResidueModification::N_TERM)
      {
        fixed_nterm_mods.push_back(*it);
      }
    }

    // @todo find OMSSA if not given
    // executable is stored in OpenMS_bin/share/OpenMS/3rdParty/OMSSA/omssacl(.exe)
    // or PATH


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    // names of temporary files for data chunks
    StringList file_spectra_chunks_in, file_spectra_chunks_out, primary_ms_runs;
    Size ms2_spec_count(0);
    { // local scope to free memory after conversion to MGF format is done
      FileHandler fh;
      FileTypes::Type in_type = fh.getType(inputfile_name);
      PeakMap peak_map;
      fh.getOptions().addMSLevel(2);
      fh.loadExperiment(inputfile_name, peak_map, in_type, log_type_, false, false);

      peak_map.getPrimaryMSRunPath(primary_ms_runs);
      ms2_spec_count = peak_map.size();
      writeDebug_("Read " + String(ms2_spec_count) + " spectra from file", 5);

      if (peak_map.getSpectra().empty())
      {
        throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS2 spectra in input file.");
      }

      // determine type of spectral data (profile or centroided)
      SpectrumSettings::SpectrumType spectrum_type = peak_map[0].getType();

      if (spectrum_type == SpectrumSettings::RAWDATA)
      {
        if (!getFlag_("force"))
        {
          throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided MS2 spectra expected. To enforce processing of the data set the -force flag.");
        }
      }

      int chunk(0);
      int chunk_size(getIntOption_("chunk_size"));
      if (chunk_size <= 0)
      {
        writeLog_("Chunk size is <=0; disabling chunking of input! If OMSSA crashes due to memory allocation errors, try setting 'chunk_size' to a value below 30000 (e.g., 10000 is usually ok).");
        chunk_size = (int) peak_map.getSpectra().size();
      }

      for (Size i = 0; i < peak_map.size(); i += chunk_size)
      {
        PeakMap map_chunk;
        PeakMap* chunk_ptr = &map_chunk; // points to the current chunk data
        // prepare a chunk
        if (static_cast<int>(peak_map.size()) <= chunk_size)
        { // we have only one chunk; avoid duplicating the whole data (could be a lot)
          // we do not use swap() since someone might want to access 'map' later and would find it empty
          chunk_ptr = &peak_map;
        }
        else
        {
          map_chunk.getSpectra().insert(map_chunk.getSpectra().begin(), peak_map.getSpectra().begin() + i, peak_map.getSpectra().begin() + std::min(
                  peak_map.size(), i + chunk_size));
        }
        MascotGenericFile omssa_infile;
        String filename_chunk = unique_input_name + String(chunk) + ".mgf";
        file_spectra_chunks_in.push_back(filename_chunk);
        writeDebug_("Storing input file: " + filename_chunk, 5);
        omssa_infile.store(filename_chunk, *chunk_ptr);
        file_spectra_chunks_out.push_back(unique_output_name + String(chunk) + ".xml");
        ++chunk;
      }
    }
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    ProteinIdentification protein_identification;
    protein_identification.setPrimaryMSRunPath(primary_ms_runs);
    vector<PeptideIdentification> peptide_ids;

    ProgressLogger pl;
    pl.setLogType(this->log_type_);
    pl.startProgress(0, file_spectra_chunks_in.size(), "OMSSA search");

    for (Size i = 0; i < file_spectra_chunks_in.size(); ++i)
    {
      pl.setProgress(i);
      StringList parameters_chunk = parameters;

      parameters_chunk << "-fm" << file_spectra_chunks_in[i];
      parameters_chunk << "-ox" << file_spectra_chunks_out[i];


      QStringList qparam;
      for (Size ip = 0; ip < parameters_chunk.size(); ++ip)
      {
        qparam << parameters_chunk[ip].toQString();
      }

      Int status = 0;
#ifdef OPENMS_WINDOWSPLATFORM
      if (db_name_contains_space)
      {
        // for some reason QProcess doesn't handle escaped " in arguments properly so we use a system call
        // see http://www.ncbi.nlm.nih.gov/books/NBK1763/ for the format the NCBI library is expecting internally if spaces are in file/path names
        String call_string = omssa_executable + " " + ListUtils::concatenate(parameters_chunk, " ");
        writeDebug_(call_string, 5);
        status = system(call_string.c_str());
      }
      else
      {
        writeDebug_(omssa_executable + " " + ListUtils::concatenate(parameters_chunk, " "), 5);
        status = QProcess::execute(omssa_executable.toQString(), qparam);
      }
#else
      writeDebug_(omssa_executable + " " + ListUtils::concatenate(parameters_chunk, " "), 5);
      status = QProcess::execute(omssa_executable.toQString(), qparam);
#endif

      if (status != 0)
      {
        writeLog_("Error: OMSSA problem! See above for OMSSA error. If this does not help, increase 'debug' level and run again.");
        writeLog_("Note: This message can also be triggered if you run out of space in your tmp directory or (32bit OMSSA only) OMSSA ran out of RAM because chunking was not used (that's the default) or 'chunk_size' was too large (>30k). Look above!");
        if (getIntOption_("debug") < 2)
        {
          QFile(file_spectra_chunks_in[i].toQString()).remove();
          QFile(file_spectra_chunks_out[i].toQString()).remove();
        }
        else
        {
          writeDebug_(String("Not removing intermediate files, but leaving them for inspection at ") + file_spectra_chunks_in[i] + " (OMSSA input) and " + file_spectra_chunks_out[i] + " (OMSSA output).\n", 2);
        }
        if (!user_mods.empty())
        {
          QFile(unique_usermod_name.toQString()).remove();
        }
        return EXTERNAL_PROGRAM_ERROR;
      }

      // read OMSSA output
      writeDebug_(String("Reading output of OMSSA from ") + file_spectra_chunks_out[i], 10);
      ProteinIdentification protein_identification_chunk;
      vector<PeptideIdentification> peptide_ids_chunk;
      OMSSAXMLFile omssa_out_file;
      omssa_out_file.setModificationDefinitionsSet(mod_set);
      // do not load empty hits for efficiency and correct stats report (below)
      omssa_out_file.load(file_spectra_chunks_out[i], protein_identification_chunk, peptide_ids_chunk, true, false);

      // OMSSA does not write fixed modifications so we need to add them to the sequences
      writeDebug_("Assigning modifications to peptides", 1);
      for (vector<PeptideIdentification>::iterator it = peptide_ids_chunk.begin(); it != peptide_ids_chunk.end(); ++it)
      {
        vector<PeptideHit> hits = it->getHits();
        for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
        {
          AASequence seq = pit->getSequence();
          for (vector<String>::const_iterator mit = fixed_nterm_mods.begin(); mit != fixed_nterm_mods.end(); ++mit)
          {
            seq.setNTerminalModification(*mit);
          }
          for (vector<String>::const_iterator mit = fixed_cterm_mods.begin(); mit != fixed_cterm_mods.end(); ++mit)
          {
            seq.setCTerminalModification(*mit);
          }
          UInt pos = 0;
          for (AASequence::Iterator mit = seq.begin(); mit != seq.end(); ++mit, ++pos)
          {
            if (fixed_residue_mods.has(mit->getOneLetterCode()))
            {
              seq.setModification(pos, fixed_residue_mods[mit->getOneLetterCode()]);
            }
          }
          pit->setSequence(seq);
        }
        it->setHits(hits);
      }

      // merge chunk results is not done, since all the statistics associated with a protein hit will be invalidated if peptide evidence is spread
      // across chunks. So we only retain this information if there is a single chunk (no splitting occurred)
      if (file_spectra_chunks_in.size() == 1)
      {
        peptide_ids = peptide_ids_chunk;
        protein_identification = protein_identification_chunk;
      }
      else
      { // add only first prot ID to have a valid id-identifier mapping (but leave hits empty)
        if (i == 0)
        {
          protein_identification = protein_identification_chunk;
          protein_identification.setHits(std::vector<ProteinHit>()); // remove hits
        }
        // ... and remove any refs from peptides
        for (vector<PeptideIdentification>::iterator it_pep = peptide_ids_chunk.begin();
             it_pep != peptide_ids_chunk.end();
             ++it_pep)
        {
          it_pep->setIdentifier(protein_identification.getIdentifier());

          // clear peptide evidences
          vector<PeptideHit> pep_hits = it_pep->getHits();
          for (vector<PeptideHit>::iterator it_pep_hit = pep_hits.begin();
             it_pep_hit != pep_hits.end();
             ++it_pep_hit)
          {
            it_pep_hit->setPeptideEvidences(std::vector<PeptideEvidence>());
          }
          it_pep->setHits(pep_hits);

          peptide_ids.push_back(*it_pep);
        }
      }

      // delete temporary files
      if (getIntOption_("debug") < 2)
      {
        writeDebug_("Removing temporary files", 10);
        QFile(file_spectra_chunks_in[i].toQString()).remove();
        QFile(file_spectra_chunks_out[i].toQString()).remove();
      }
      else
      {
        writeDebug_(String("Not removing intermediate files, but leaving them for inspection at ") + file_spectra_chunks_in[i] + " (OMSSA input) and " + file_spectra_chunks_out[i] + " (OMSSA output).\n", 2);
      }
    } // chunks loop

    pl.endProgress();

    if (getIntOption_("debug") <= 1)
    {
      writeDebug_("Removing temporary files", 10);
      if (!user_mods.empty())
      {
        QFile(unique_usermod_name.toQString()).remove();
      }
    }

    // handle the search parameters
    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = getStringOption_("database");
    //search_parameters.db_version =
    search_parameters.taxonomy = getStringOption_("x");
    search_parameters.charges = "+" + String(getIntOption_("min_precursor_charge")) + "-+" + String(getIntOption_("max_precursor_charge"));
    ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;

    if (getIntOption_("tom") == 1)
    {
      mass_type = ProteinIdentification::AVERAGE;
    }
    else
    {
      if (getIntOption_("tom") != 0)
      {
        writeLog_("Warning: unrecognized mass type: " + String(getIntOption_("tom")));
      }
    }
    search_parameters.mass_type = mass_type;
    search_parameters.fixed_modifications = getStringList_("fixed_modifications");
    search_parameters.variable_modifications = getStringList_("variable_modifications");
    search_parameters.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(getStringOption_("enzyme")));
    search_parameters.missed_cleavages = getIntOption_("v");
    search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");
    search_parameters.precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
    search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor_error_units") == "ppm";
    search_parameters.fragment_mass_tolerance_ppm = false; // OMSSA doesn't support ppm fragment mass tolerance

    protein_identification.setSearchParameters(search_parameters);
    protein_identification.setSearchEngineVersion(omssa_version);
    protein_identification.setSearchEngine("OMSSA");

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    vector<ProteinIdentification> protein_identifications;
    protein_identifications.push_back(protein_identification);
    IdXMLFile().store(outputfile_name, protein_identifications, peptide_ids);

    // some stats
    LOG_INFO << "Statistics:\n"
             << "  identified MS2 spectra: " << peptide_ids.size() << " / " << ms2_spec_count << " = " << int(peptide_ids.size() * 100.0 / ms2_spec_count) << "% (with e-value < " << String(getDoubleOption_("he")) << ")" << std::endl;


    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPOMSSAAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
