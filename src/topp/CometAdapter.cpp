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
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

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
    @page TOPP_CometAdapter CometAdapter

    @brief Identifies peptides in MS/MS spectra via Comet.

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


class TOPPCometAdapter :
  public TOPPBase
{
public:
  TOPPCometAdapter() :
    TOPPBase("CometAdapter", "Annotates MS/MS spectra using Comet.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {

    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerInputFile_("database", "<file>", "", "FASTA file", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));
    registerInputFile_("comet_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      "comet.exe",
      "Comet executable of the installation e.g. 'comet.exe'", true, false, ListUtils::create<String>("skipexists"));
    registerStringOption_("comet_version","<choice>", "2016.01 rev. 2","comet version: (year,version,revision)",false,false);               //required as first line in the param file
    setValidStrings_("comet_version", ListUtils::create<String>("2016.01 rev. 2,2016.01 rev. 3,2017.01 rev. 0beta"));
    //
    // Optional parameters //
    //

    //Files
    registerOutputFile_("pin_out", "<file>", "", "Output file - for Percolator input", false);
    setValidFormats_("pin_out", ListUtils::create<String>("csv"));
    registerInputFile_("default_params_file", "<file>", "", "Default Comet params file. All parameters of this take precedence. A template file can be generated using comet.exe -p", false, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("default_params_file", ListUtils::create<String>("txt"));
    //registerIntOption_("threads", "<num>", 1, "number of threads", false, true);

    //Masses
    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 10.0, "Precursor monoisotopic mass tolerance (Comet parameter: peptide_mass_tolerance)", false, false);
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "peptide_mass_units 0=amu, 1=mmu, 2=ppm", false, false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("amu,ppm,Da"));
    //registerIntOption_("mass_type_parent", "<num>", 1, "0=average masses, 1=monoisotopic masses", false, true);
    //registerIntOption_("mass_type_fragment", "<num>", 1, "0=average masses, 1=monoisotopic masses", false, true);
    //registerIntOption_("precursor_tolerance_type", "<num>", 0, "0=average masses, 1=monoisotopic masses", false, false);
    registerStringOption_("isotope_error", "<choice>", "off", "0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)", false, false);
    setValidStrings_("isotope_error", ListUtils::create<String>("off,-1/0/1/2/3,-8/-4/0/4/8"));

    //Search Enzyme
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllCometNames(all_enzymes);
    registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false, false);
    setValidStrings_("enzyme", all_enzymes);
    registerStringOption_("num_enzyme_termini", "<choice>", "fully", "1 semi-digested, 2 fully digested, (default), 8 C-term unspecific, 9 N-term unspecific", false, false);
    setValidStrings_("num_enzyme_termini", ListUtils::create<String>("semi,fully,C-term unspecific,N-term unspecific"));
    registerIntOption_("allowed_missed_cleavages", "<num>", 1, "Number of possible cleavage sites missed by the enzyme, maximum value is 5; for enzyme search", false, false);

    //Fragment Ions
    registerDoubleOption_("fragment_bin_tolerance", "<tolerance>", 1.0005, "fragment_mass_tolerance (MSGF+), fragment_bin_tol (Comet)", false, true);
    registerDoubleOption_("fragment_bin_offset", "<tolerance>", 0.25, "fragment_bin_offset (Comet)", false, true);
    registerStringOption_("instrument", "<choice>", "high_res", "comets theoretical_fragment_ions parameter: theoretical fragment ion peak representation, high res ms/ms: sum of intensities plus flanking bins, ion trap ms/ms: sum of intensities of central bin only", false, true);
    setValidStrings_("instrument", ListUtils::create<String>("low_res,high_res"));
    registerStringOption_("use_A_ions", "<num>", "false", "use A ions for PSM", false, true);
    setValidStrings_("use_A_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_B_ions", "<num>", "true", "use B ions for PSM, 0 == no, 1 == yes", false, true);
    setValidStrings_("use_B_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_C_ions", "<num>", "false", "use C ions for PSM, 0 == no, 1 == yes", false, true);
    setValidStrings_("use_C_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_X_ions", "<num>", "false", "use X ions for PSM, 0 == no, 1 == yes", false, true);
    setValidStrings_("use_X_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_Y_ions", "<num>", "true", "use Y ions for PSM, 0 == no, 1 == yes", false, true);
    setValidStrings_("use_Y_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_Z_ions", "<num>", "false", "use Z ions for PSM, 0 == no, 1 == yes", false, true);
    setValidStrings_("use_Z_ions", ListUtils::create<String>("true,false"));
    registerStringOption_("use_NL_ions", "<num>", "false", "use Neutral Loss ions for PSM, 0 == no, 1 == yes", false, true);
    setValidStrings_("use_NL_ions", ListUtils::create<String>("true,false"));

    //Output
    registerIntOption_("num_hits", "<num>", 5, "Number of peptide hits in output file", false, false);

    //mzXML/mzML parameters
    registerStringOption_("precursor_charge", "0:0", ":", "charge range to search: 0 0 == search all charges, 2 6 == from +2 to +6, 3 3 == +3", false, false);
    registerStringOption_("override_charge", "<choice>", "keep any known", "0 = keep any known precursor charge state, 1 = ignore known precursor charge state and use precursor_charge parameter, 2 = ignore precursor charges outside precursor_charge range, 3 = keep any known precursor charge state. For unknown charge states, search as singly charged if there is no signal above the precursor m/z or use the precursor_charge range", false, false);
    setValidStrings_("override_charge", ListUtils::create<String>("keep any known,ignore known,ignore outside range,keep known search unknown"));
    registerIntOption_("ms_level", "<num>", 2, "MS level to analyze, valid are levels 2 (default) or 3", false, false);
    setMinInt_("ms_level",2);
    setMaxInt_("ms_level",3);
    registerStringOption_("activation_method", "<method>", "ALL", "activation method; used if activation method set; allowed ALL, CID, ECD, ETD, PQD, HCD, IRMPD", false, false);
    setValidStrings_("activation_method", ListUtils::create<String>("ALL,CID,ECD,ETD,PQD,HCD,IRMPD"));

    //Misc. parameters
    registerStringOption_("digest_mass_range", "600:5000", ":", "MH+ peptide mass range to analyze", false, true);
    registerIntOption_("max_fragment_charge", "<num>", 3, "set maximum fragment charge state to analyze (allowed max 5)", false, false);
    registerStringOption_("max_precursor_charge", "<num>", "0+", "set maximum precursor charge state to analyze (allowed max 9)", false, true);
    registerStringOption_("clip_nterm_methionine", "<num>", "false", "0=leave sequences as-is; 1=also consider sequence w/o N-term methionine", false, false);
    setValidStrings_("clip_nterm_methionine", ListUtils::create<String>("true,false"));
    registerIntOption_("spectrum_batch_size", "<num>", 1000, "max. // of spectra to search at a time; 0 to search the entire scan range in one loop", false, true);
    registerDoubleOption_("mass_offsets", "<offset>", 0, "one or more mass offsets to search (values subtracted from deconvoluted precursor mass)", false, true);

    // spectral processing
    registerIntOption_("minimum_peaks", "<num>", 10, "required minimum number of peaks in spectrum to search (default 10)", false, true);
    registerIntOption_("minimum_intensity", "<num>", 0, "minimum intensity value to read in", false, true);
    registerStringOption_("remove_precursor_peak", "<choice>", "no", "0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)", false, true);
    setValidStrings_("remove_precursor_peak", ListUtils::create<String>("no,yes,all"));
    registerIntOption_("remove_precursor_tolerance", "<num>", 1.5, "+- Da tolerance for precursor removal", false, true);
    registerStringOption_("clear_mz_range", "0:0", ":", "for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range", false, true);

    //Modifications
    registerStringList_("fixed_modifications", "<mods>", vector<String>(), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false, false);
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", vector<String>(), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false, false);
    setValidStrings_("variable_modifications", all_mods);
    registerIntOption_("max_variable_mods_in_peptide", "<num>", 5, "", false, true);
  }

  vector<ResidueModification> getModifications_(StringList modNames)
  {
    vector<ResidueModification> modifications;

    // iterate over modification names and add to vector
    for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
    {
      String modification(*mod_it);
      modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
    }

    return modifications;
  }

  void removeTempDir_(const String& tmp_dir)
  {
    if (tmp_dir.empty()) {return;} // no temporary directory created

    if (debug_level_ >= 2)
    {
      writeDebug_("Keeping temporary files in directory '" + tmp_dir + "'. Set debug level to 1 or lower to remove them.", 2);
    }
    else
    {
      if (debug_level_ == 1) 
      {
        writeDebug_("Deleting temporary directory '" + tmp_dir + "'. Set debug level to 2 or higher to keep it.", 1);
      }
      File::removeDirRecursively(tmp_dir);
    }
  }

  void createParamFile_(ostream& os)
  {
    os << "# comet_version " << getStringOption_("comet_version") << "\n";               //required as first line in the param file
    os << "# Comet MS/MS search engine parameters file.\n";
    os << "# Everything following the '#' symbol is treated as a comment.\n";
    os << "database_name = " << getStringOption_("database") << "\n";
    os << "decoy_search = " << 0 << "\n"; // 0=no (default), 1=concatenated search, 2=separate search
    os << "num_threads = " << getIntOption_("threads") << "\n";  // 0=poll CPU to set num threads; else specify num threads directly (max 64)

    // masses
    map<String,int> precursor_error_units;
    precursor_error_units["amu"] = 0;
    precursor_error_units["mmu"] = 1;
    precursor_error_units["ppm"] = 2;

    map<string,int> isotope_error;
    isotope_error["off"] = 0;
    isotope_error["-1/0/1/2/3"] = 1;
    isotope_error["-8/-4/0/4/8"] = 2;

    os << "peptide_mass_tolerance = " << getDoubleOption_("precursor_mass_tolerance") << "\n";
    os << "peptide_mass_units = " << precursor_error_units[getStringOption_("precursor_error_units")] << "\n";                  // 0=amu, 1=mmu, 2=ppm
    os << "mass_type_parent = " << 1 << "\n";                    // 0=average masses, 1=monoisotopic masses
    os << "mass_type_fragment = " << 1 << "\n";                  // 0=average masses, 1=monoisotopic masses
    os << "precursor_tolerance_type = " << 0 << "\n";            // 0=MH+ (default), 1=precursor m/z; only valid for amu/mmu tolerances
    os << "isotope_error = " << isotope_error[getStringOption_("isotope_error")] << "\n";                      // 0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)

    // search enzyme

    String enzyme_name = getStringOption_("enzyme");
    String enzyme_number = String(ProteaseDB::getInstance()->getEnzyme(enzyme_name)->getCometID());

    map<string,int> num_enzyme_termini;
    num_enzyme_termini["semi"] = 1;
    num_enzyme_termini["fully"] = 2;
    num_enzyme_termini["C-term unspecific"] = 8;
    num_enzyme_termini["N-term unspecific"] = 9;

    os << "search_enzyme_number = " << enzyme_number << "\n";                // choose from list at end of this params file
    os << "num_enzyme_termini = " << num_enzyme_termini[getStringOption_("num_enzyme_termini")] << "\n";                  // 1 (semi-digested), 2 (fully digested, default), 8 C-term unspecific , 9 N-term unspecific
    os << "allowed_missed_cleavage = " << getIntOption_("allowed_missed_cleavages") << "\n";             // maximum value is 5; for enzyme search

    // Up to 9 variable modifications are supported
    // format:  <mass> <residues> <0=variable/else binary> <max_mods_per_peptide> <term_distance> <n/c-term> <required>
    //     e.g. 79.966331 STY 0 3 -1 0 0
    vector<String> variable_modifications_names = getStringList_("variable_modifications");
    vector<ResidueModification> variable_modifications = getModifications_(variable_modifications_names);
    if (variable_modifications.size() > 9)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Comet only supports 9 variable modifications. " + String(variable_modifications.size()) + " provided.");
    }

    Size var_mod_index = 0;

    // write out user specified modifications
    for (; var_mod_index < variable_modifications.size(); ++var_mod_index)
    {
      const ResidueModification mod = variable_modifications[var_mod_index];
      double mass = mod.getDiffMonoMass();
      String residues = mod.getOrigin();  // TODO: check if origin contains C-term string or similar. Should not be passed to comet as residue string
      String variable = "0";
      String max_mods_per_peptide = "3";
      String term_distance = "-1";
      String nc_term = "0";

      if (residues=="N-term")
      {
          residues="n";
      }
      if (residues=="C-term")
      {
          residues="c";
      }
      if (mod.getTermSpecificity() == ResidueModification::C_TERM)
      {
        term_distance = 0;
        nc_term = "3";
      }
      else if (mod.getTermSpecificity() == ResidueModification::N_TERM)
      {
        term_distance = 0;
        nc_term = "2";
      }
      else if (mod.getTermSpecificity() == ResidueModification::PROTEIN_N_TERM) // not yet available
      {
        term_distance = 0;
        nc_term = "0";
      }
      else if (mod.getTermSpecificity() == ResidueModification::PROTEIN_C_TERM) // not yet available
      {
        term_distance = 0;
        nc_term = "1";
      }
      String required = "0";

      os << "variable_mod0" << var_mod_index+1 << " = " << mass << " " << residues << " " << variable << " " << max_mods_per_peptide << " " << term_distance << " " << nc_term << " " << required << "\n";
    }

    // fill remaining modification slots (if any) in Comet with "no modification"
    for (; var_mod_index < 9; ++var_mod_index)
    {
      os << "variable_mod0" << var_mod_index+1 << " = " << "0.0 X 0 3 -1 0 0" << "\n";
    }

    os << "max_variable_mods_in_peptide = " << getIntOption_("max_variable_mods_in_peptide") << "\n";
    os << "require_variable_mod = " << 0 << "\n";

    // fragment ions
    // ion trap ms/ms:  1.0005 tolerance, 0.4 offset (mono masses), theoretical_fragment_ions = 1
    // high res ms/ms:    0.02 tolerance, 0.0 offset (mono masses), theoretical_fragment_ions = 0

    map<string,int> instrument;
    instrument["high_res"] = 0;
    instrument["low_res"] = 1;

    os << "fragment_bin_tol = " << getDoubleOption_("fragment_bin_tolerance") << "\n";               // binning to use on fragment ions
    os << "fragment_bin_offset = " << getDoubleOption_("fragment_bin_offset")  << "\n";              // offset position to start the binning (0.0 to 1.0)
    os << "theoretical_fragment_ions = " << instrument[getStringOption_("instrument")] << "\n";           // 0=use flanking peaks, 1=M peak only
    os << "use_A_ions = " << (int)(getStringOption_("use_A_ions")=="true") << "\n";
    os << "use_B_ions = " << (int)(getStringOption_("use_B_ions")=="true") << "\n";
    os << "use_C_ions = " << (int)(getStringOption_("use_C_ions")=="true") << "\n";
    os << "use_X_ions = " << (int)(getStringOption_("use_X_ions")=="true") << "\n";
    os << "use_Y_ions = " << (int)(getStringOption_("use_Y_ions")=="true") << "\n";
    os << "use_Z_ions = " << (int)(getStringOption_("use_Z_ions")=="true") << "\n";
    os << "use_NL_ions = " << (int)(getStringOption_("use_NL_ions")=="true") << "\n";                         // 0=no, 1=yes to consider NH3/H2O neutral loss peaks

    // output
    os << "output_sqtstream = " << 0 << "\n";                    // 0=no, 1=yes  write sqt to standard output
    os << "output_sqtfile = " << 0 << "\n";                      // 0=no, 1=yes  write sqt file
    os << "output_txtfile = " << 0 << "\n";                     // 0=no, 1=yes  write tab-delimited txt file
    os << "output_pepxmlfile = " << 1 << "\n";                   // 0=no, 1=yes  write pep.xml file

    os << "output_percolatorfile = " << !getStringOption_("pin_out").empty() << "\n";              // 0=no, 1=yes  write Percolator tab-delimited input file
    os << "output_outfiles = " <<  0 << "\n";                    // 0=no, 1=yes  write .out files
    os << "print_expect_score = " << 1 << "\n";                  // 0=no, 1=yes to replace Sp with expect in out & sqt
    os << "num_output_lines = " << getIntOption_("num_hits") << "\n";                    // num peptide results to show
    os << "show_fragment_ions = " << 0 << "\n";                  // 0=no, 1=yes for out files only
    os << "sample_enzyme_number = " << enzyme_number << "\n";                // Sample enzyme which is possibly different than the one applied to the search.

    // mzXML parameters
    map<string,int> override_charge;
    override_charge["keep any known"] = 0;
    override_charge["ignore known"] = 1;
    override_charge["ignore outside range"] = 2;
    override_charge["keep known search unknown"] = 3;

    int precursor_charge_min, precursor_charge_max;
    parseRange_(getStringOption_("precursor_charge"), precursor_charge_min, precursor_charge_max);

    os << "scan_range = " << "0 0" << "\n";                        // start and scan scan range to search; 0 as 1st entry ignores parameter
    os << "precursor_charge = " << precursor_charge_min << " " << precursor_charge_max << "\n";                  // precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter
    os << "override_charge = " << override_charge[getStringOption_("override_charge")] << "\n";                     // 0=no, 1=override precursor charge states, 2=ignore precursor charges outside precursor_charge range, 3=see online
    os << "ms_level = " << getIntOption_("ms_level") << "\n";                            // MS level to analyze, valid are levels 2 (default) or 3
    os << "activation_method = " << getStringOption_("activation_method") << "\n";                 // activation method; used if activation method set; allowed ALL, CID, ECD, ETD, PQD, HCD, IRMPD

    // misc parameters
    int digest_mass_range_min, digest_mass_range_max;
    parseRange_(getStringOption_("digest_mass_range"), digest_mass_range_min, digest_mass_range_max);

    os << "digest_mass_range = " << digest_mass_range_min << " " << digest_mass_range_max << "\n";        // MH+ peptide mass range to analyze
    os << "num_results = " << 100 << "\n";                       // number of search hits to store internally
    os << "skip_researching = " << 1 << "\n";                    // for '.out' file output only, 0=search everything again (default), 1=don't search if .out exists
    os << "max_fragment_charge = " << getIntOption_("max_fragment_charge") << "\n";                 // set maximum fragment charge state to analyze (allowed max 5)
    os << "max_precursor_charge = " << getStringOption_("max_precursor_charge") << "\n";                // set maximum precursor charge state to analyze (allowed max 9)
    os << "nucleotide_reading_frame = " << 0 << "\n";            // 0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six
    os << "clip_nterm_methionine = " << (int)(getStringOption_("clip_nterm_methionine")=="true") << "\n";              // 0=leave sequences as-is; 1=also consider sequence w/o N-term methionine
    os << "spectrum_batch_size = " << getIntOption_("spectrum_batch_size") << "\n";                 // max. // of spectra to search at a time; 0 to search the entire scan range in one loop
    os << "decoy_prefix = " << "rev_" << "\n";                 // decoy entries are denoted by this string which is pre-pended to each protein accession
    os << "output_suffix = " << "" << "\n";                      // add a suffix to output base names i.e. suffix "-C" generates base-C.pep.xml from base.mzXML input
    os << "mass_offsets = " << getDoubleOption_("mass_offsets") << "\n";                       // one or more mass offsets to search (values subtracted from deconvoluted precursor mass)

    // spectral processing
    map<string,int> remove_precursor_peak;
    remove_precursor_peak["no"] = 0;
    remove_precursor_peak["yes"] = 1;
    remove_precursor_peak["all"] = 2;

    double clear_mz_range_min, clear_mz_range_max;
    parseRange_(getStringOption_("clear_mz_range"), clear_mz_range_min, clear_mz_range_max);

    os << "minimum_peaks = " << getIntOption_("minimum_peaks") << "\n";                      // required minimum number of peaks in spectrum to search (default 10)
    os << "minimum_intensity = " << getIntOption_("minimum_intensity") << "\n";                   // minimum intensity value to read in
    os << "remove_precursor_peak = " << remove_precursor_peak[getStringOption_("remove_precursor_peak")] << "\n";               // 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)
    os << "remove_precursor_tolerance = " << getIntOption_("remove_precursor_tolerance") << "\n";        // +- Da tolerance for precursor removal
    os << "clear_mz_range = " << clear_mz_range_min << " " << clear_mz_range_max << "\n";                // for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range


    // write fixed modifications - if not specified residue parameter is zero
    // Aminoacid:
    //      add_AA.OneletterCode_AA.ThreeLetterCode = xxx
    // Terminus:
    //      add_N/Cterm_peptide = xxx       protein not available yet
    vector<String> fixed_modifications_names = getStringList_("fixed_modifications");
    vector<ResidueModification> fixed_modifications = getModifications_(fixed_modifications_names);
    for (vector<ResidueModification>::const_iterator it = fixed_modifications.begin(); it != fixed_modifications.end(); ++it)
    {
      String AA = it->getOrigin();
      if ((AA!="N-term") && (AA!="C-term"))
      {
      const Residue* r = ResidueDB::getInstance()->getResidue(AA);
      String name = r->getName();
      os << "add_" << r->getOneLetterCode() << "_" << name.toLower() << " = " << it->getDiffMonoMass() << endl;
      }
      else
      {
      os << "add_" << AA.erase(1,1) << "_peptide = " << it->getDiffMonoMass() << endl;
      }
    }

    // COMET_ENZYME_INFO _must_ be at the end of this parameters file
    os << "[COMET_ENZYME_INFO]" << "\n";
    os << "0.  No_enzyme              0      -           -" << "\n";
    os << "1.  Trypsin                1      KR          P" << "\n";
    os << "2.  Trypsin/P              1      KR          -" << "\n";
    os << "3.  Lys_C                  1      K           P" << "\n";
    os << "4.  Lys_N                  0      K           -" << "\n";
    os << "5.  Arg_C                  1      R           P" << "\n";
    os << "6.  Asp_N                  0      D           -" << "\n";
    os << "7.  CNBr                   1      M           -" << "\n";
    os << "8.  Glu_C                  1      DE          P" << "\n";
    os << "9.  PepsinA                1      FL          P" << "\n";
    os << "10. Chymotrypsin           1      FWYL        P" << "\n";
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String inputfile_name = getStringOption_("in");
    writeDebug_(String("Input file: ") + inputfile_name, 1);
    if (inputfile_name.empty())
    {
      writeLog_("No input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    String out = getStringOption_("out");
    writeDebug_(String("Output file___real one: ") + out, 1);
    if (out.empty())
    {
      writeLog_("No output file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    String db_name(getStringOption_("database"));
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

    //tmp_dir
    //const String tmp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/").toQString());
    const String tmp_dir = makeTempDirectory_(); //OpenMS::File::getTempDirectory() + "/";
    writeDebug_("Creating temporary directory '" + tmp_dir + "'", 1);
    //QDir d;
    //d.mkpath(tmp_dir.toQString());
    String tmp_pepxml = tmp_dir + "result.pep.xml";
    String tmp_pin = tmp_dir + "result.pin";
    String default_params = getStringOption_("default_params_file");
    String tmp_file;

    //default params given or to be written
    if (default_params.empty())
    {
        tmp_file = tmp_dir + "param.txt";
        ofstream os(tmp_file.c_str());
        createParamFile_(os);
        os.close();
    }
    else
    {
        tmp_file = default_params;
    }

    PeakMap exp;
    MzMLFile mzml_file;
    mzml_file.getOptions().addMSLevel(2); // only load msLevel 2 //TO DO: setMSLevels or clearMSLevels
    mzml_file.setLogType(log_type_);
    mzml_file.load(inputfile_name, exp);

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

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    String paramP = "-P" + tmp_file;
    String paramN = "-N" + File::removeExtension(File::removeExtension(tmp_pepxml));
    QStringList process_params;
    process_params << paramP.toQString() << paramN.toQString() << inputfile_name.toQString();
    qDebug() << process_params;

    String comet_executable = getStringOption_("comet_executable");
    int status = QProcess::execute(comet_executable.toQString(),process_params); // does automatic escaping etc...
    if (status != 0)
    {
      writeLog_("Comet problem. Aborting! Calling command was: '" + comet_executable + " \"" + inputfile_name + "\"'.\nDoes the Comet executable exist?");
      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // writing IdXML output
    //-------------------------------------------------------------

    // read the pep.xml put of Comet and write it to idXML

    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;

    writeDebug_("load PepXMLFile", 1);
    PepXMLFile().load(tmp_pepxml, protein_identifications, peptide_identifications);
    writeDebug_("write idXMLFile", 1);
    writeDebug_(out, 1);
    IdXMLFile().store(out, protein_identifications, peptide_identifications);

    //-------------------------------------------------------------
    // create (move) optional pin output
    //-------------------------------------------------------------

    String pin_out = getStringOption_("pin_out");
    if (!pin_out.empty())
    {
      // existing file? Qt won't overwrite, so try to remove it:
      if (QFile::exists(pin_out.toQString()) && !QFile::remove(pin_out.toQString()))
      {
        writeLog_("Fatal error: Could not overwrite existing file '" + pin_out + "'");
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      // move the temporary file to the actual destination:
      if (!QFile::rename(tmp_pin.toQString(), pin_out.toQString()))
      {
        writeLog_("Fatal error: Could not move temporary mzid file to '" + pin_out + "'");
        return CANNOT_WRITE_OUTPUT_FILE;
      }
    }

    // remove tempdir
    if (this->debug_level_ == 0)
    {
        removeTempDir_(tmp_dir);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << tmp_dir << "'" << std::endl;
    }
    else
    {
      LOG_WARN << "Keeping the temporary files at '" << tmp_pepxml << "'. Set debug level to 0 to remove them." << std::endl;
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPCometAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
