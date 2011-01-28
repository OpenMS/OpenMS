// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/MascotInfile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <iostream>

#include <QtCore/QFile>
#include <QtCore/QProcess>

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
  an error message, try another OMSSA version

	Sequence databases in fasta format must be converted into the NCBI format before OMSSA can read them. Therefore, use the program formatdb
	of the NCBI-tools suite (see ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.13/ for a working version).
  The latest NCBI BLAST distribution does not contain the formatdb executable any longer!).
  Use @em formatdb @em -i @em SwissProt_TargetAndDecoy.fasta @em -o to create
  additional files, which	will be used by @em OMSSA. The database option of the @em OMSSAAdapter should contain the name of the psq file
  , e.g., 'SwissProt_TargetAndDecoy.fasta.psq'. The '.psq' suffix can also be omitted, e.g. 'SwissProt_TargetAndDecoy.fasta' and will be added
  automatically.
  This makes it easy to specifiy a common TOPPAS input node (using only the FASTA suffix) for many adapters.
  
  This adapter supports relative database filenames, which (when not found in the current working directory) is looked up in
  the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

	The options that specify the protease specificity (@em e) are directly taken from OMSSA. A complete list of available
	proteases can be found be executing @em omssacl @em -el.

	This wrapper has been tested successfully with OMSSA, version 2.x.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_OMSSAAdapter.cli

	@improvement modes to read OMSSA output data and save in idXML format (Andreas)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPOMSSAAdapter
	: public TOPPBase
{
	public:
		TOPPOMSSAAdapter()
			: TOPPBase("OMSSAAdapter","Annotates MS/MS spectra using OMSSA.")
		{
		}

	protected:

    struct OMSSAVersion
    {
      OMSSAVersion ()
        : omssa_major(0), omssa_minor(0), omssa_patch(0)
      {}

      OMSSAVersion (Int maj, Int min, Int pat)
        : omssa_major(maj), omssa_minor(min), omssa_patch(pat)
      {}

      Int omssa_major;
      Int omssa_minor;
      Int omssa_patch;

      bool operator < (const OMSSAVersion& v) const
      {
        if (omssa_major > v.omssa_major) return false;
        else if (omssa_major < v.omssa_major) return true;
        else // ==
        {
          if (omssa_minor > v.omssa_minor) return false;
          else if (omssa_minor < v.omssa_minor) return true;
          else
          {
            return (omssa_patch < v.omssa_patch);
          }
        }

      }
    };

    bool getVersion_(const String& version, OMSSAVersion& omssa_version_i) const
    {
      // we expect three components 
      IntList nums = IntList::create(StringList::create(version,'.'));
      if (nums.size()!=3) return false;

      omssa_version_i.omssa_major =nums[0];
      omssa_version_i.omssa_minor =nums[1];
      omssa_version_i.omssa_patch =nums[2];
      return true;
    }

		void registerOptionsAndFlags_()
		{

			addEmptyLine_();
			addText_("Common Identification engine options");

			registerInputFile_("in", "<file>", "", "input file ");
			setValidFormats_("in",StringList::create("mzML"));
      registerOutputFile_("out", "<file>", "", "output file ");
	  	setValidFormats_("out",StringList::create("idXML"));

      registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 1.5, "precursor mass tolerance (Default: Dalton)", false);
      registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.3, "fragment mass error in Dalton", false);
      registerFlag_("precursor_mass_tolerance_unit_ppm", "If this flag is set, ppm is used as precursor mass tolerance unit");
      registerInputFile_("database", "<psq-file>", "", "NCBI formatted fasta files. Only the psq filename should be given, e.g. 'SwissProt.fasta.psq'. If the filename does not end in '.psq' the suffix will be added automatically. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, StringList::create("skipexists"));
			registerIntOption_("min_precursor_charge", "<charge>", 1, "minimum precursor ion charge", false);
      registerIntOption_("max_precursor_charge", "<charge>", 3, "maximum precursor ion charge", false);
			vector<String> all_mods;
			ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
      registerStringList_("fixed_modifications", "<mods>", StringList::create(""), "fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
			setValidStrings_("fixed_modifications", all_mods);
      registerStringList_("variable_modifications", "<mods>", StringList::create(""), "variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
			setValidStrings_("variable_modifications", all_mods);

			addEmptyLine_();
			addText_("OMSSA specific input options");

			//Sequence library
			//-d <String> Blast sequence library to search.  Do not include .p* filename suffixes.
			//-pc <Integer> The number of pseudocounts to add to each precursor mass bin.
			//registerStringOption_("d", "<file>", "", "Blast sequence library to search.  Do not include .p* filename suffixes", true);
      registerInputFile_("omssa_executable", "<executable>", "omssacl", "The 'omssacl' executable of the OMSSA installation", true, false, StringList::create("skipexists"));
      registerInputFile_("omssa_user_mods", "<file>", "", "additional <MSModSpec> subtrees of user modifications.\nSubtrees will be pasted into OMSSAAdapter generated user mod files.\nSee http://www.ncbi.nlm.nih.gov/data_specs/schema/OMSSA.mod.xsd for details about user mod file definition.", false, true, StringList::create("input file"));
			registerIntOption_("pc", "<Integer>", 1, "The number of pseudocounts to add to each precursor mass bin", false, true);

			//registerFlag_("omssa_out", "If this flag is set, the parameter 'in' is considered as an output file of OMSSA and will be converted to IdXML");
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
			//-to <Real> product ion mass tolerance in Da
			//-te <Real> precursor ion mass tolerance in Da
			//-tez <Integer> scaling of precursor mass tolerance with charge (0 = none, 1= linear)
			//registerDoubleOption_("to", "<Real>", 0.8, "product ion mass tolerance in Da", false);
			//registerDoubleOption_("te", "<Real>", 2.0, "precursor ion mass tolerance in Da", false);
			registerIntOption_("tez", "<Integer>", 1, "scaling of precursor mass tolerance with charge (0 = none, 1= linear)", false, true);

			//A precursor ion is the ion before fragmentation and the product ions are the ions generated after fragmentation. These values are specified in Daltons +/- the measured value, e.g. a value of 2.0 means +/- 2.0 Daltons of the measured value.
			//The tez value allows you to specify how the mass tolerance scales with the charge of the precursor. For example, you may search a precursor assuming that it has a charge state of 2+ and 3+. If you set tez to 1, then the mass tolerance for the +2 charge state will be 2 times the precursor mass tolerance, and for the 3+ charge state it will be 3 times the precursor mass tolerance. If you set tez to 0, the mass tolerance will always be equal to the precursor mass tolerance, irrespective of charge state.

			//-tom <Integer> product ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact.
			//-tem <Integer> precursor ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact.
			registerIntOption_("tom", "<Integer>", 0, "product ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact", false, true);
			registerIntOption_("tem", "<Integer>", 0, "precursor ion search type, with 0 = monoisotopic, 1 = average, 2 = monoisotopic N15, 3 = exact", false, true);

			//Monoisotopic searching searches spectral peaks that correspond to peptides consisting entirely of carbon-12. Average mass searching searches on the average natural isotopic mass of peptides. Exact mass searches on the most abundant isotopic peak for a given mass range.

			//-tex <Double> threshold in Da above which the mass of a neutron should be added in an exact mass search.
			registerDoubleOption_("tex", "<Real>", 1446.94, "threshold in Da above which the mass of a neutron should be added in an exact mass search", false, true);

			//Preprocessing
			//Preprocessing is the process of eliminating noise from a spectrum. Normally, you do not need to adjust these options as OMSSA automatically adjusts its preprocessing for best results.

			//-cl <Real> low intensity cutoff as a fraction of max peak
			//-ch <Real> high intensity cutoff as a fraction of max peak
			//-ci <Real> intensity cutoff increment as a fraction of max peak
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
			registerDoubleOption_("z1", "<Real>", 0.95, "the fraction of peaks below the precursor used to determine if the spectrum is charge +1", false, true);
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
			registerIntOption_("e", "<Integer>", 0, "id number of enzyme to use (trypsin is the default)", false);
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
			registerIntOption_("hl", "<Integer>", 30, "maximum number of hits retained for one spectrum", false);
			registerDoubleOption_("he", "<Real>", 1, "the maximum e-value allowed in the hit list", false);

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
			registerDoubleOption_("is", "<Real>", 0.0, "evalue threshold to include a sequence in the iterative search, 0 = all", false, true);
			registerDoubleOption_("ir", "<Real>", 0.0, "evalue threshold to replace a hit, 0 = only if better", false, true);
			registerDoubleOption_("ii", "<Real>", 0.0, "evalue threshold to iteratively search a spectrum again, 0 = always", false, true);


			//-foms <String> read in search result in .oms format (binary asn.1).
			//-fomx <Double> read in search result in .omx format (xml).
			//Iterative searching is the ability to re-search search results in hopes of increasing the number of spectra identified. To accomplish this, an iterative search may change search parameters, such as using a no-enzyme search, or restrict the sequence search library to sequences already hit.

		}

		ExitCodes main_(int , const char**)
		{
			// instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
			String ini_location;
			// path to the log file
			String logfile(getStringOption_("log"));
			String omssa_executable(getStringOption_("omssa_executable"));
			String inputfile_name;
			String outputfile_name;
			PeakMap map;

			String parameters;
			String unique_name = File::getUniqueName(); // body for the tmp files
			String unique_input_name = unique_name + "_OMSSA.mgf";
			String unique_output_name = unique_name + "_OMSSA.xml";
			String unique_version_name = unique_name + "_OMSSA_version";
			String unique_usermod_name = unique_name + "_OMSSA_user_mod_file.xml";

			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------

			// get version of OMSSA
      QProcess qp;
      String call = omssa_executable + " -version";
      qp.start(call.toQString(), QIODevice::ReadOnly); // does automatic escaping etc...
      bool success = qp.waitForFinished();
      String output (QString(qp.readAllStandardOutput ()));
			String omssa_version;
      OMSSAVersion omssa_version_i;
      if (!success || qp.exitStatus() != 0 || qp.exitCode()!=0)
			{
        writeLog_("Warning: unable to determine the version of OMSSA - the process returned an error. Call string was: '" + call + "'. Make sure that the path to the OMSSA executable is correct!");
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
			inputfile_name = getStringOption_("in");
			outputfile_name = getStringOption_("out");
      String db_name = String(getStringOption_("database"));
      // @todo: find DB for OMSSA (if not given) in OpenMS_bin/share/OpenMS/DB/*.fasta|.pin|...


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
			    printUsage_();
			    return ILLEGAL_PARAMETERS;
        }
        db_name = full_db_name;
      }
      
      db_name = db_name.substr(0,db_name.size()-4); // OMSSA requires the filename without the .psq part

			parameters += " -d "  +  db_name;
			parameters += " -to " +  String(getDoubleOption_("fragment_mass_tolerance")); //String(getDoubleOption_("to"));
			parameters += " -hs " + String(getIntOption_("hs"));
			parameters += " -te " +  String(getDoubleOption_("precursor_mass_tolerance")); //String(getDoubleOption_("te"));
      if (getFlag_("precursor_mass_tolerance_unit_ppm"))
      {
        if (omssa_version_i < OMSSAVersion(2,1,8))
        {
          writeLog_("This OMSSA version (" + omssa_version + ") does not support the 'precursor_mass_tolerance_unit_ppm' flag."
                   +" Please disable it and set the precursor tolerance in Da."
                   +" Required version is 2.1.8 and above.\n");
          return ILLEGAL_PARAMETERS;
        }
        parameters += " -teppm "; // only from OMSSA 2.1.8 on
      }
			parameters += " -zl " +  String(getIntOption_("min_precursor_charge")); //String(getIntOption_("zl"));
			parameters += " -zh " +  String(getIntOption_("max_precursor_charge")); //String(getIntOption_("zh"));
			parameters += " -zt " +  String(getIntOption_("zt"));
			parameters += " -zc " +  String(getIntOption_("zc"));
			parameters += " -zcc " + String(getIntOption_("zcc"));
			parameters += " -zoh " + String(getIntOption_("zoh"));
			parameters += " -no " + String(getIntOption_("no"));
			parameters += " -nox " + String(getIntOption_("nox"));
			parameters += " -sp " + String(getIntOption_("sp"));
			parameters += " -sb1 " + String(getIntOption_("sb1"));
			parameters += " -sct " + String(getIntOption_("sct"));
			parameters += " -x " + getStringOption_("x");
			parameters += " -hl " + String(getIntOption_("hl"));
			parameters += " -hm " + String(getIntOption_("hm"));
			parameters += " -ht " +  String(getIntOption_("ht"));
			parameters += " -tex " + String(getDoubleOption_("tex"));
			parameters += " -i " + getStringOption_("i");
			parameters += " -z1 " + String(getDoubleOption_("z1"));
			parameters += " -v " +   String(getIntOption_("v"));
			parameters += " -e " +   String(getIntOption_("e"));
			parameters += " -tez " + String(getIntOption_("tez"));


			parameters += " -tom " + String(getIntOption_("tom"));
			parameters += " -tem " + String(getIntOption_("tem"));

			parameters += " -mm " + String(getIntOption_("mm"));
			parameters += " -is " + String(getDoubleOption_("is"));
			parameters += " -ir " + String(getDoubleOption_("ir"));
			parameters += " -ii " + String(getDoubleOption_("ii"));
			parameters += " -nt " + String(getIntOption_("threads"));

			if (getFlag_("mnm"))
			{
				parameters += " -mnm ";
			}

			parameters += " -fm " + unique_input_name;
			parameters += " -ox " + unique_output_name;

			if (getIntOption_("debug") == 0)
			{
				parameters += " -ni ";
			}
			parameters += " -he " + String(getDoubleOption_("he"));


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
          	throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "parse mapping file line: '" + *it + "'", "");
        	}
        	vector<ResidueModification> mods;
        	for (Size i = 2; i != split.size(); ++i)
        	{
          	String tmp(split[i].trim());
          	if (tmp.size() != 0)
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
			if (getStringList_("fixed_modifications").size() != 0)
			{
				set<String> mod_names = mod_set.getFixedModificationNames();
				String mod_list;
				for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
				{
					if (mods_map.has(*it))
					{
						if (mod_list != "")
						{
							mod_list += ",";
						}
						mod_list += String(mods_map[*it]);
					}
					else
					{
						if (mod_list != "")
						{
							mod_list += ",";
						}
						mod_list += String(user_mod_num);

						// add this to the usermods
						user_mods.push_back(make_pair(user_mod_num++, *it));
            writeDebug_("Inserting unknown fixed modification: '" + *it + "' into OMSSA", 1);
					}
				}
				if (mod_list != "")
				{
					parameters += " -mf " + mod_list;
				}
			}

			if (getStringList_("variable_modifications").size() != 0)
			{
				set<String> mod_names = mod_set.getVariableModificationNames();
				String mod_list;

        for (set<String>::const_iterator it = mod_names.begin(); it != mod_names.end(); ++it)
        {
          if (mods_map.has(*it))
          {
            if (mod_list != "")
            {
              mod_list += ",";
            }
            mod_list += String(mods_map[*it]);
          }
          else
          {
            if (mod_list != "")
            {
              mod_list += ",";
            }
            mod_list += String(user_mod_num);

            // add this to the usermods
            user_mods.push_back(make_pair(user_mod_num++, *it));
            writeDebug_("Inserting unknown variable modification: '" + *it + "' into OMSSA", 1);
          }
        }

				if (mod_list != "")
				{
					parameters += " -mv " + mod_list;
				}
			}

      String additional_user_mods_filename = getStringOption_("omssa_user_mods");
			// write unknown modifications to user mods file
      if (user_mods.size() != 0 || additional_user_mods_filename != "")
			{
				writeDebug_("Writing usermod file to " + unique_usermod_name, 1);
				parameters += " -mux " + File::absolutePath(unique_usermod_name);
				ofstream out(unique_usermod_name.c_str());
				out << "<?xml version=\"1.0\"?>" << endl;
				out << "<MSModSpecSet xmlns=\"http://www.ncbi.nlm.nih.gov\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\" xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov OMSSA.xsd\">" << endl;

				UInt user_mod_count(1);
				for (vector<pair<UInt, String> >::const_iterator it = user_mods.begin(); it != user_mods.end(); ++it)
				{
					writeDebug_("Writing information into user mod file of modification: " + it->second, 1);
					out << "<MSModSpec>" << endl;
					out << "\t<MSModSpec_mod>" << endl;
					out << "\t\t<MSMod value=\"usermod" << user_mod_count++ << "\">" << it->first << "</MSMod>" << endl;
					out << "\t</MSModSpec_mod>" << endl;
					out << "\t<MSModSpec_type>" << endl;

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

					ResidueModification::Term_Specificity ts = ModificationsDB::getInstance()->getModification(it->second).getTermSpecificity();
					String origin = ModificationsDB::getInstance()->getModification(it->second).getOrigin();
					if (ts == ResidueModification::ANYWHERE)
					{
						out << "\t\t<MSModType value=\"modaa\">0</MSModType>" << endl;
					}
					if (ts == ResidueModification::C_TERM)
					{
						if (origin == "" || origin == "X")
						{
							out << "\t\t<MSModType value=\"modcp\">7</MSModType>" << endl;
						}
						else
						{
							out << "\t\t<MSModType value=\"modcpaa\">8</MSModType>" << endl;
						}
					}
					if (ts == ResidueModification::N_TERM)
					{
						if (origin == "" || origin == "X")
						{
							out << "\t\t<MSModType value=\"modnp\">5</MSModType>" << endl;
						}
						else
						{
							out << "\t\t<MSModType value=\"modnpaa\">6</MSModType>" << endl;
						}
					}
					out << "\t</MSModSpec_type>" << endl;

					out << "\t<MSModSpec_name>" << it->second << "</MSModSpec_name>" << endl;
					out << "\t<MSModSpec_monomass>" << ModificationsDB::getInstance()->getModification(it->second).getDiffMonoMass()  << "</MSModSpec_monomass>" << endl;
					out << "\t<MSModSpec_averagemass>" << ModificationsDB::getInstance()->getModification(it->second).getDiffAverageMass() << "</MSModSpec_averagemass>" << endl;
					out << "\t<MSModSpec_n15mass>0</MSModSpec_n15mass>" << endl;

					if (origin != "")
					{
						out << "\t<MSModSpec_residues>" << endl;
						out << "\t\t<MSModSpec_residues_E>" << origin << "</MSModSpec_residues_E>" << endl;
						out << "\t</MSModSpec_residues>" << endl;
						out << "</MSModSpec>" << endl;
					}
				}

        // Add additional MSModSPec subtrees to generated user mods
        ifstream additional_user_mods_file(additional_user_mods_filename.c_str());
        String line;
        if(additional_user_mods_file.is_open())
        {
          while (additional_user_mods_file.good())
          {
            getline(additional_user_mods_file, line);
            out << line << endl;
          }
          additional_user_mods_file.close();
        }
				out << "</MSModSpecSet>" << endl;
				out.close();
			}

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			MzMLFile mzml_infile;
			mzml_infile.setLogType(log_type_);
			ProteinIdentification protein_identification;
			vector<PeptideIdentification> peptide_ids;
			mzml_infile.load(inputfile_name, map);

			writeDebug_("Read " + String(map.size()) + " spectra from file", 5);

			vector<ProteinIdentification> protein_identifications;
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------

			writeDebug_("Storing input file: " + unique_input_name, 5);
			MascotInfile omssa_infile;
			omssa_infile.store(unique_input_name, map, "OMSSA search tmp file");

      // @todo find OMSSA if not given
      // executable is stored in OpenMS_bin/share/OpenMS/3rdParty/OMSSA/omssacl(.exe)
      // or PATH

			writeDebug_("omssa_executable " + parameters, 5);
      Int status = QProcess::execute(omssa_executable.toQString(), QStringList(parameters.toQString().split(" ", QString::SkipEmptyParts))); // does automatic escaping etc...
			if (status != 0)
			{
				writeLog_("Error: OMSSA problem! (Details can be seen in the logfile: \"" + logfile + "\")");

				QFile(unique_input_name.toQString()).remove();
				QFile(unique_output_name.toQString()).remove();
        if (user_mods.size() != 0 || additional_user_mods_filename!="")
				{
					QFile(unique_usermod_name.toQString()).remove();
				}
				return EXTERNAL_PROGRAM_ERROR;
			}

			// read OMSSA output
			writeDebug_("Reading output of OMSSA", 10);
			OMSSAXMLFile omssa_out_file;
			omssa_out_file.setModificationDefinitionsSet(mod_set);
			omssa_out_file.load(unique_output_name, protein_identification, peptide_ids);

			// OMSSA does not write fixed modifications so we need to add them to the sequences
			set<String> fixed_mod_names = mod_set.getFixedModificationNames();
			vector<String> fixed_nterm_mods, fixed_cterm_mods;
			Map<String, String> fixed_residue_mods;
			writeDebug_("Splitting modification into N-Term, C-Term and anywhere specificity", 1);
			for (set<String>::const_iterator it = fixed_mod_names.begin(); it != fixed_mod_names.end(); ++it)
			{
				ResidueModification::Term_Specificity ts = ModificationsDB::getInstance()->getModification(*it).getTermSpecificity();
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
			writeDebug_("Assigning modifications to peptides", 1);
			for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
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

			// delete temporary files
			writeDebug_("Removing temporary files", 10);
			QFile(unique_input_name.toQString()).remove();
			QFile(unique_output_name.toQString()).remove();
			if (user_mods.size() != 0)
			{
				QFile(unique_usermod_name.toQString()).remove();
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
			ProteinIdentification::DigestionEnzyme enzyme = ProteinIdentification::TRYPSIN;

			UInt e(getIntOption_("e"));
			if (e != 0)
			{
				writeLog_("Warning: cannot handle enzyme: " + getIntOption_("e"));
			}

			search_parameters.enzyme = enzyme;
			search_parameters.missed_cleavages = getIntOption_("v");
			search_parameters.peak_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");
			search_parameters.precursor_tolerance = getDoubleOption_("precursor_mass_tolerance");


			protein_identification.setSearchParameters(search_parameters);
			protein_identification.setSearchEngineVersion(omssa_version);
			protein_identification.setSearchEngine("OMSSA");

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------

			protein_identifications.push_back(protein_identification);
			IdXMLFile().store(outputfile_name, protein_identifications, peptide_ids);

			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPOMSSAAdapter tool;

	return tool.main(argc,argv);
}

/// @endcond
