// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Authors: Leon Bichmann, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>
#include <QDebug>
#include <iostream>
#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_MSFraggerAdapter MSFraggerAdapter

    @brief TODO 

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ CometAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em Comet must be installed before this wrapper can be used. This wrapper
    has been successfully tested with version 2016.01.2 of Comet.

    Comet settings not exposed by this adapter can be directly adjusted using a param file, which can be generated using comet -p.
    By default, All (!) parameters available explicitly via this param file will take precedence over the wrapper parameters.

    Parameter names have been changed to match names found in other search engine adapters, however some are Comet specific.
    For a detailed description of all available parameters check the Comet documentation at http://comet-ms.sourceforge.net/parameters/parameters_201601/
    The default parameters are set for a high resolution instrument.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_CometAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_CometAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPMSFraggerAdapter :
  public TOPPBase
{
public:
	  static const String param_executable;
	  static const String param_in;
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
	  static const String param_




	  static const String param_allow_multiple_variable_mods_on_residue;
	  static const String param_max_variable_mods_per_mod;
	  static const String param_max_variable_mods_combinations;
	  static const String param_precursor_charge_min;
	  static const String param_precursor_charge_max;
	  static const String param_override_charge;
	  static const String param_max_fragment_charge;
	  static const String param_track_zero_topn;
	  static const String param_zero_bin_accept_expect;
	  static const String param_zero_bin_mult_expect;
	  static const String param_add_topn_complementary;



  TOPPMSFraggerAdapter() :
    TOPPBase("MSFraggerAdapter", "Annotates MS/MS spectra using MSFragger.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
	 const StringList empty;
	 const StringList validUnits = ListUtils::create<String>("Da,ppm");
	 const StringList isotope_error_and_enzyme_termini = ListUtils::create<String>("0,1,2");
	 const StringList zero_to_five = ListUtils::create<String>("0,1,2,3,4,5");

	 // Handle executable
	 registerInputFile_(TOPPMSFraggerAdapter::param_executable, "<path_to_executable>", "MSFragger.jar", "Path to the MSFragger executable to use; may be empty if the executable is globally available.", true, false, ListUtils::create<String>("skipexists"));

	 // input spectra
	 registerInputFileList_(TOPPMSFraggerAdapter::param_in, "<spectra_files>", empty, "Spectra files to search with MSFragger", true, false);
	 setValidFormats_(TOPPMSFraggerAdapter::param_in, ListUtils::create<String>("mzML,mzXML"));

	 // Path to database to search
	 registerInputFile_(TOPPMSFraggerAdapter::param_database, "<path_to_fasta>", "", "FASTA file", true, false);
	 setValidFormats_(TOPPMSFraggerAdapter::param_database, ListUtils::create<String>("FASTA,fasta,fa,fas"), false);

	 // TOPP tolerance
	 registerTOPPSubsection_("tolerance", "Search Tolerances");

	 // Precursor mass tolerance and unit
	 registerDoubleOption_(TOPPMSFraggerAdapter::param_precursor_mass_tolerance, "<precursor_mass_tolerance>", 20.0, "Precursor mass tolerance", false, false);
	 registerStringOption_(TOPPMSFraggerAdapter::param_precursor_mass_unit, "<precursor_mass_unit>", "ppm", "Unit of precursor mass tolerance", false, false);
	 setValidStrings_(TOPPMSFraggerAdapter::param_precursor_mass_unit, validUnits);

	 // Precursor true tolerance
	 registerDoubleOption_(TOPPMSFraggerAdapter::param_precursor_true_tolerance, "<precursor_true_tolerance>", 0.0, "Precursor true tolerance. It is strongly recommended to set this feature!", false, false);
	 registerStringOption_(TOPPMSFraggerAdapter::param_precursor_true_unit, "<precursor_true_unit>", "ppm", "Unit of precursor true tolerance", false, false);
	 setValidStrings_(TOPPMSFraggerAdapter::param_precursor_true_unit, validUnits);

	 // Fragment mass tolerance
	 registerDoubleOption_(TOPPMSFraggerAdapter::param_fragment_mass_tolerance, "<fragment_mass_tolerance>", 20.0, "Fragment mass tolerance", false, false);
	 registerStringOption_(TOPPMSFraggerAdapter::param_fragment_mass_unit, "<fragment_mass_unit>", "ppm", "Fragment mass unit", false, false);
	 setValidStrings_(TOPPMSFraggerAdapter::param_fragment_mass_unit, validUnits);

	 // Isotope error
	 registerStringOption_(TOPPMSFraggerAdapter::param_isotope_error, "<isotope_error>", "0", "Isotope correction for MS/MS events triggered on isotopic peaks. Should be set to 0 for open search or 0/1/2 for correction of narrow window searches", false, false);
	 setValidStrings_(TOPPMSFraggerAdapter::param_isotope_error, isotope_error_and_enzyme_termini);

	 // TOPP digest
	 registerTOPPSubsection_("digest", "In-Silico Digestion Parameters");

	 // Enzyme
	 StringList enzyme_names;
	 EnzymesDB::getInstance()->getAllNames(enzyme_names);
	 registerStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_name, "<search_enzyme_name>", "Trypsin", "Name of the enzyme to be written to the pepXML file", false, false);
	 setValidStrings_(TOPPMSFraggerAdapter::param_search_enzyme_name, enzyme_names);

	 // Cut after
	registerStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_cutafter, "<search_enzyme_cutafter>", "KR", "Residues after which the enzyme cuts", false , false);

	// No cut before
	registerStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_nocutbefore, "<search_enzyme_nocutbefore>", "P", "Residues that the enzyme will not cut before", false, false);

	// Number of enzyme termini
	registerStringOption_(TOPPMSFraggerAdapter::param_num_enzyme_termini, "<num_enzyme_termini>", "2", "Number of enzyme termini (0, 1, or 2 for non-enzymatic, semi-enzymatic, fully-enzymatic)", false, false);
	setValidStrings_(TOPPMSFraggerAdapter::param_num_enzyme_termini, isotope_error_and_enzyme_termini);

	// Allowed missed cleavages
	registerStringOption_(TOPPMSFraggerAdapter::param_allowed_missed_cleavage, "<allowed_missed_cleavage>", "2", "Allowed number of missed cleavages", false, false);
	setValidStrings_(TOPPMSFraggerAdapter::param_allowed_missed_cleavage, zero_to_five); // 5 is the max. allowed value according to MSFragger

	// Digest min length
	registerIntOption_(TOPPMSFraggerAdapter::param_digest_min_length, "<digest_min_length>", 7, "Minimum length of peptides to be generated during in-silico digestion", false, false);

	// Digest max length
	registerIntOption_(TOPPMSFraggerAdapter::param_digest_max_length, "<digest_max_length>", 64, "Maximum length of peptides to be generated during in-silico digestion", false, false);

	// Digest min mass range
	registerDoubleOption_(TOPPMSFraggerAdapter::param_digest_mass_range_min, "<digest_mass_range_min>", 500.0, "Min mass of peptides to be generated (Da)", false, false);

	// Digest max mass range
	registerDoubleOption_(TOPPMSFraggerAdapter::param_digest_mass_range_max, "<digest_mass_range_max>", 5000.0, "Max mass of peptides to be generated (Da)", false, false);


	// TOPP varmod
	registerTOPPSubsection_("varmod", "Variable Modification Parameters");

	// Clip nterm M
	registerFlag_(TOPPMSFraggerAdapter::param_clip_nterm_m, "Specifies the trimming of a protein N-terminal methionine as a variable modification", false);

	// TODO Modifications



	// allow_multiple_variable_mods_on_residue
	registerFlag_(TOPPMSFraggerAdapter::param_allow_multiple_variable_mods_on_residue, "Allow each amino acid to be modified by multiple variable modifications", false);

	// Max variable mods per mod
	registerStringOption_(TOPPMSFraggerAdapter::param_max_variable_mods_per_mod, "<max_variable_mods_per_mod>", "2", "Maximum number of residues that can be occupied by each variable modification", false, true);
	setValidStrings_(TOPPMSFraggerAdapter::param_max_variable_mods_per_mod, zero_to_five);

	// Max variable mods combinations
	registerIntOption_(TOPPMSFraggerAdapter::param_max_variable_mods_combinations, "<max_variable_mods_combinations>", 5000, "Maximum allowed number of modified variably modified peptides from each peptide sequence, (maximum of 65534). If a greater number than the maximum is generated, only the unmodified peptide is considered.", false, true);

	// Minimum charge of precursor
	registerIntOption_(TOPPMSFraggerAdapter::param_precursor_charge_min, "<precursor_charge_min>", 1, "Min charge of precursor charge range to consider" , false, false);

	// Maximum charge of precursor
	registerIntOption_(TOPPMSFraggerAdapter::param_precursor_charge_max, "<precursor_charge_max>", 4, "Max charge of precursor charge range to consider" , false, false);

	// Override charge
	registerFlag_(TOPPMSFraggerAdapter::param_override_charge, "Ignores precursor charge and uses charge state specified in precursor_charge range" , false);


	// Max fragment charge
	registerStringOption_(TOPPMSFraggerAdapter::param_max_fragment_charge, "<max_fragment_charge>", "2", "Maximum charge state for theoretical fragments to match", false, false);
	setValidStrings_(TOPPMSFraggerAdapter::param_max_fragment_charge, ListUtils::create<String>("1,2,3,4"));

	// track zero top N
	registerIntOption_(TOPPMSFraggerAdapter::param_track_zero_topn, "<track_zero_topn>", 0, "Track top N unmodified peptide results separately from main results internally for boosting features. Should be set to a number greater than output_report_topN if zero bin boosting is desired.", false, true);

	// zero_bin_accept_expect
	registerDoubleOption_(TOPPMSFraggerAdapter::param_zero_bin_accept_expect, "<zero_bin_accept_expect>", 0.0, "Ranks a zero-bin hit above all non-zero-bin hit if it has expectation less than this value", false, false);



	// Output file
	// TODO


   // registerOutputFile_("out", "<file>", "", "Output file");
   // setValidFormats_("out", ListUtils::create<String>("idXML"));
    //
    // Optional parameters //
    //

    //Files
    //registerInputFile_("default_params_file", "<file>", "", "Default Comet params file. All parameters of this take precedence. A template file can be generated using comet.exe -p", false, false, ListUtils::create<String>("skipexists"));
    //setValidFormats_("default_params_file", ListUtils::create<String>("txt"));
    //Masses
   //registerStringOption_("isotope_error", "<choice>", "off", "Isotope correction for MS/MS events triggered on isotopic peaks. Should be set to off for open search or 0/1/2 for correction of narrow window searches. Shifts the precursor mass window to multiples of this value multiplied by the mass of C13-C12.", false, false);
  //  setValidStrings_("isotope_error", ListUtils::create<String>("off,0/1/2"));

    //Search Enzyme
    //vector<String> all_enzymes;
    // EnzymesDB::getInstance()->getAllMSFraggerNames(all_enzymes);
  //  registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false, false);
  //  setValidStrings_("enzyme", all_enzymes);
    //search_enzyme_cutafter = KR
    //search_enzyme_butnotafter = P             
  //  registerStringOption_("num_enzyme_termini", "<choice>", "fully", "semi-digested, fully digested (default), non-specific", false, false); // 2 for enzymatic, 1 for semi-enzymatic, 0 for nonspecific digestion
    //setValidStrings_("num_enzyme_termini", ListUtils::create<String>("semi,fully,C-term unspecific,N-term unspecific"));
   // registerIntOption_("allowed_missed_cleavages", "<num>", 1, "Number of possible cleavage sites missed by the enzyme, maximum value is 5; for enzyme search", false, false);
    //registerIntOption_("min_peptide_length", "<num>", 6, "Minimum peptide length to consider (MSFragger parameter 'digest_min_length')", false);
   // setMinInt_("min_peptide_length", 1);
   // registerIntOption_("max_peptide_length", "<num>", 40, "Maximum peptide length to consider (MSFragger parameter 'digest_max_length')", false);
   // setMinInt_("max_peptide_length", 1);
   // registerStringOption_("digest_mass_range", "600:5000", ":", "MH+ peptide mass range to analyze", false, true);

    //Output
    //registerIntOption_("num_hits", "<num>", 5, "Number of peptide hits in output file", false, false);

    //mzXML/mzML parameters
   // registerStringOption_("precursor_charge", "0:0", ":", "charge range to search: 0 0 == search all charges, 2 6 == from +2 to +6, 3 3 == +3", false, false);
   // registerStringOption_("override_charge", "<choice>", "keep any known", "0 = keep any known precursor charge state, 1 = ignore known precursor charge state and use precursor_charge parameter, 2 = ignore precursor charges outside precursor_charge range, 3 = keep any known precursor charge state. For unknown charge states, search as singly charged if there is no signal above the precursor m/z or use the precursor_charge range", false, false);
   // setValidStrings_("override_charge", ListUtils::create<String>("keep any known,ignore known,ignore outside range,keep known search unknown"));
   // registerIntOption_("ms_level", "<num>", 2, "MS level to analyze, valid are levels 2 (default) or 3", false, false);
   // setMinInt_("ms_level",2);
   // setMaxInt_("ms_level",3);
   // registerStringOption_("activation_method", "<method>", "ALL", "activation method; used if activation method set; allowed ALL, CID, ECD, ETD, PQD, HCD, IRMPD", false, false);
    //setValidStrings_("activation_method", ListUtils::create<String>("ALL,CID,ECD,ETD,PQD,HCD,IRMPD"));

    //Misc. parameters
   // registerIntOption_("max_fragment_charge", "<num>", 3, "set maximum fragment charge state to analyze (allowed max 5)", false, false);
   // registerStringOption_("max_precursor_charge", "<num>", "0+", "set maximum precursor charge state to analyze (allowed max 9)", false, true);
   // registerStringOption_("clip_nterm_methionine", "<num>", "false", "0=leave sequences as-is; 1=also consider sequence w/o N-term methionine", false, false);
   // setValidStrings_("clip_nterm_methionine", ListUtils::create<String>("true,false"));
   // registerDoubleOption_("mass_offsets", "<offset>", 0, "one or more mass offsets to search (values substracted from deconvoluted precursor mass)", false, true);

    // spectral processing
   // registerIntOption_("minimum_peaks", "<num>", 10, "required minimum number of peaks in spectrum to search (default 10)", false, true);
   // registerIntOption_("minimum_intensity", "<num>", 0, "minimum intensity value to read in", false, true);
   // registerStringOption_("remove_precursor_peak", "<choice>", "no", "0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)", false, true);
   // setValidStrings_("remove_precursor_peak", ListUtils::create<String>("no,yes,all"));
   // registerIntOption_("remove_precursor_tolerance", "<num>", 1.5, "+- Da tolerance for precursor removal", false, true);
   // registerStringOption_("clear_mz_range", "0:0", ":", "for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range", false, true);

    //Modifications
   // registerStringList_("fixed_modifications", "<mods>", vector<String>(), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false, false);
   // vector<String> all_mods;
   // ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
   // setValidStrings_("fixed_modifications", all_mods);
   // registerStringList_("variable_modifications", "<mods>", vector<String>(), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false, false);
   // setValidStrings_("variable_modifications", all_mods);
   // registerIntOption_("max_variable_mods_in_peptide", "<num>", 5, "Maximum number of residues that can be occupied by each variable modification (maximum of 5)", false, true);
   // registerStringOption_("allow_multiple_variable_mods_on_residue", "<choice>", "no", "Allow each amino acid to be modified by multiple variable modifications", false, true);
   // setValidStrings_("allow_multiple_variable_mods_on_residue", ListUtils::create<String>("no,yes"));
  //  registerIntOption_("max_variable_mods_combinations", "<num>", 5000, "maximum of 65534, limits number of modified peptides generated from sequence", false, true);
  }


  ExitCodes main_(int, const char**)
  {
   try
   {
	   // Get and validate parameters

	   // tolerance
	   const double arg_precursor_mass_tolerance(this->_getNonNegativeDouble(TOPPMSFraggerAdapter::param_precursor_mass_tolerance));
	   const String & arg_precursor_mass_unit = this->getStringOption_(TOPPMSFraggerAdapter::param_precursor_mass_unit);
	   const double arg_precursor_true_tolerance(this->_getNonNegativeDouble(TOPPMSFraggerAdapter::param_precursor_true_tolerance));
	   const String & arg_precursor_true_unit = this->getStringOption_(TOPPMSFraggerAdapter::param_precursor_true_unit);
	   const double arg_fragment_mass_tolerance(this->_getNonNegativeDouble(TOPPMSFraggerAdapter::param_fragment_mass_tolerance));
	   const String & arg_fragment_mass_unit = this->getStringOption_(TOPPMSFraggerAdapter::param_fragment_mass_unit);
	   const String & arg_isotope_error = this->getStringOption_(TOPPMSFraggerAdapter::param_isotope_error);

	   // digest
	   const String & arg_search_enzyme_name = this->getStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_name);
	   const String & arg_search_enzyme_cutafter = this->getStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_cutafter);
	   const String & arg_search_enzyme_nocutbefore = this->getStringOption_(TOPPMSFraggerAdapter::param_search_enzyme_nocutbefore);
	   const String & arg_num_enzyme_termini = this->getStringOption_(TOPPMSFraggerAdapter::param_num_enzyme_termini);
	   const String & arg_allowed_missed_cleavage = this>getStringOption_(TOPPMSFraggerAdapter::param_allowed_missed_cleavage);
	   const int arg_digest_min_length = this->_getNonNegativeInt(TOPPMSFraggerAdapter::param_digest_min_length);
	   const int arg_digest_max_length = this->_getNonNegativeInt(TOPPMSFraggerAdapter::param_digest_max_length);
	   if (arg_digest_max_length < arg_digest_min_length)
	   {
	     LOG_ERROR << "FATAL: Maximum length of digest is not allowed to be smaller than minimum length of digest" << std::endl;
		 throw 1;
	   }
	   const double arg_digest_mass_range_min = this->_getNonNegativeDouble(TOPPMSFraggerAdapter::param_digest_mass_range_min);
	   const double arg_digest_mass_range_max = this->_getNonNegativeDouble(TOPPMSFraggerAdapter::param_digest_mass_range_max);

	   if (arg_digest_mass_range_max < arg_digest_mass_range_min)
	   {
	     LOG_ERROR << "FATAL: Maximum digest mass is not allowed to be smaller than minimum digest mass!" << std::endl;
	     throw 1;
	   }

	   // varmod
	   const bool arg_clip_nterm_m = this->getFlag_(TOPPMSFraggerAdapter::param_clip_nterm_m);



	   const bool arg_allow_multiple_variable_mods_on_residue = this->getFlag_(TOPPMSFraggerAdapter::param_allow_multiple_variable_mods_on_residue);
	   const String & arg_max_variable_mods_per_mod  = this->getStringOption_(TOPPMSFraggerAdapter::param_max_variable_mods_per_mod);
	   const int arg_max_variable_mods_combinations = this->_getInRangeInt(TOPPMSFraggerAdapter::param_max_variable_mods_combinations, 0, 65534);
	   const int arg_precursor_charge_min = this->_getNonNegativeInt(TOPPMSFraggerAdapter::param_precursor_charge_min);
	   const int arg_precursor_charge_max = this->_getNonNegativeInt(TOPPMSFraggerAdapter::param_precursor_charge_max);
	   if (arg_precursor_charge_max < arg_precursor_charge_min)
	   {
		   LOG_ERROR << "FATAL: Maximum precursor charge is not allowed to be smaller than minimum precursor charge!" << std::endl;
		   throw 1;
	   }
	   const bool arg_override_charge = this->getFlag_(TOPPMSFraggerAdapter::param_override_charge);






	   const String & arg_max_fragment_charge = this->getStringOption_(TOPPMSFraggerAdapter::param_max_fragment_charge);
	   const int arg_track_zero_topn = this->_getNonNegativeInt(TOPPMSFraggerAdapter::param_track_zero_topn);
	   const double arg_zero_bin_accept_expect = this->_getNonNegativeDouble(TOPPMSFraggerAdapter::param_zero_bin_accept_expect);
	   const double arg_zero_bin_mult_expect = this->_getNonNegativeDouble(TOPPMSFraggerAdapter::param_zero_bin_mult_expect);
	   const int arg_add_topn_complementary = this->_getNonNegativeInt(TOPPMSFraggerAdapter::param_add_topn_complementary);

   }
   catch (...)
   {
	   return ILLEGAL_PARAMETERS;
   }







  return EXECUTION_OK;
  }


private:

  double _getNonNegativeDouble(const String & param_name) const
  {
	 const double arg_value(this->getDoubleOption_(param_name));
	 if (arg_value < 0)
	 {
		 LOG_ERROR << "FATAL: Parameter " << param_name << " is not allowed to be negative!" << std::endl;
		 throw 1;
	 }
	 return arg_value;
  }

  int _getNonNegativeInt(const String & param_name) const
  {
	  const int arg_value(this->getIntOption_(param_name));
	  if (arg_value < 0)
	  {
		  LOG_ERROR << "FATAL: Parameter " << param_name << " is not allowed to be negative!" << std::endl;
		  throw 1;
	  }
	  return arg_value;
  }

  int _getInRangeInt(const String & param_name, const int min_value, const int max_value) const
  {
	  const int arg_value(this->getIntOption_(param_name));
	  if (arg_value < min_value || arg_value > max_value)
	  {
		  LOG_ERROR << "FATAL: Parameter " << param_name << " is out of range: " << min_value << "..." << max_value << std::endl;
		  throw 1;
	  }
	  return arg_value;
  }


};

const String TOPPMSFraggerAdapter::param_executable = "executable";
const String TOPPMSFraggerAdapter::param_in = "in";
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
const String TOPPMSFraggerAdapter::param_digest_min_length = "digest:digest_min_length";
const String TOPPMSFraggerAdapter::param_digest_max_length = "digest:digest_max_length";
const String TOPPMSFraggerAdapter::param_digest_mass_range_min = "digest:digest_mass_range_min";
const String TOPPMSFraggerAdapter::param_digest_mass_range_max = "digest:digest_mass_range_max";

// varmod
const String TOPPMSFraggerAdapter::param_clip_nterm_m = "varmod:clip_nterm_m";





const String TOPPMSFraggerAdapter::param_allow_multiple_variable_mods_on_residue = "allow_multiple_variable_mods_on_residue";
const String TOPPMSFraggerAdapter::param_max_variable_mods_per_mod = "max_variable_mods_per_mod";
const String TOPPMSFraggerAdapter::param_max_variable_mods_combinations = "max_variable_mods_combinations";
const String TOPPMSFraggerAdapter::param_precursor_charge_min = "precursor_charge_min";
const String TOPPMSFraggerAdapter::param_precursor_charge_max = "precursor_charge_max";
const String TOPPMSFraggerAdapter::param_override_charge = "override_charge";

const String TOPPMSFraggerAdapter::param_max_fragment_charge = "max_fragment_charge";
const String TOPPMSFraggerAdapter::param_track_zero_topn = "track_zero_topn";
const String TOPPMSFraggerAdapter::param_zero_bin_accept_expect = "zero_bin_accept_expect";
const String TOPPMSFraggerAdapter::param_zero_bin_mult_expect = "zero_bin_mult_expect";
const String TOPPMSFraggerAdapter::param_add_topn_complementary = "add_topn_complementary";


int main(int argc, const char** argv)
{
  TOPPMSFraggerAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
