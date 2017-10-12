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

  TOPPMSFraggerAdapter() :
    TOPPBase("MSFraggerAdapter", "Annotates MS/MS spectra using MSFragger.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
	 const StringList empty;

	 // Handle executable
	 registerInputFile_(TOPPMSFraggerAdapter::param_executable, "<path_to_executable>", "MSFragger.jar", "Path to the MSFragger executable to use; may be empty if the executable is globally available.", true, false, ListUtils::create<String>("skipexists"));

	 // input spectra
	 registerInputFileList_(TOPPMSFraggerAdapter::param_in, "<spectra_files>", empty, "Spectra files to search with MSFragger", true, false);
	 setValidFormats_(TOPPMSFraggerAdapter::param_in, ListUtils::create<String>("mzML,mzXML"));


   // registerOutputFile_("out", "<file>", "", "Output file");
   // setValidFormats_("out", ListUtils::create<String>("idXML"));
    //registerInputFile_("database", "<file>", "", "FASTA file", true, false, ListUtils::create<String>("skipexists"));
    //setValidFormats_("database", ListUtils::create<String>("FASTA"));
    //registerStringOption_("MSFragger_version","<choice>", "2016.01 rev. 2","comet version: (year,version,revision)",false,false);               //required as first line in the param file
    //setValidStrings_("MSFragger_version", ListUtils::create<String>("2016.01 rev. 2"));
    //
    // Optional parameters //
    //

    //Files
    //registerInputFile_("default_params_file", "<file>", "", "Default Comet params file. All parameters of this take precedence. A template file can be generated using comet.exe -p", false, false, ListUtils::create<String>("skipexists"));
    //setValidFormats_("default_params_file", ListUtils::create<String>("txt"));
    //registerIntOption_("threads", "<num>", 1, "number of threads", false, true);

    //Masses
    //registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 20.0, "Precursor mass tolerance (MSFragger)", false, false);
    //registerStringOption_("precursor_error_units", "<choice>", "ppm", "precursor_mass_units 0=Dalton, 1=ppm", false, false);
    //setValidStrings_("precursor_error_units", ListUtils::create<String>("Da,ppm"));
    //registerDoubleOption_("precursor_true_tolerance", "<tolerance>", 0, "Used for tie breaker of results (in spectrally ambiguous cases) and zero bin boosting in open searches (0 disables these features). This option is STRONGLY recommended for open searches.", false, false);
    //registerStringOption_("precursor_true_units", "<choice>", "ppm", "precursor_mass_units 0=Dalton, 1=ppm", false, false);
    //setValidStrings_("precursor_true_units", ListUtils::create<String>("Da,ppm"));
   // registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 6, "fragment_mass_tolerance (MSGF+), fragment_bin_tol (Comet)", false, true);
    //registerDoubleOption_("fragment_mass_units", "<tolerance>", 0.25, "fragment_mass_units 0=Dalton, 1=ppm", false, true);
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

    return EXECUTION_OK;
  }

};


const String TOPPMSFraggerAdapter::param_executable = "executable";
const String TOPPMSFraggerAdapter::param_in = "in";

int main(int argc, const char** argv)
{
  TOPPMSFraggerAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
