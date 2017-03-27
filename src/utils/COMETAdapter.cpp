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
// $Authors: Leon Bichmann, Timo Sachsenberg, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
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
#include <iostream>     // std::cout, std::ostream, std::ios
#include <fstream>

using namespace OpenMS;         
using namespace std;             

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_COMETAdapter COMETAdapter

    @brief Identifies peptides in MS/MS spectra via COMET.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ COMETAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em COMET must be installed before this wrapper can be used. This wrapper
    has been successfully tested with several versions of COMET.

    This adapter supports relative database filenames, which (when not found in the current working directory) is looked up in
    the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

    COMET settings not exposed by this adapter can be directly adjusted using an XML configuration file.
    By default, all (!) parameters available explicitly via this wrapper take precedence over the XML configuration file.
    The parameter "default_config_file" can be used to specify such a custom configuration.
    An example of a configuration file (named "default_input.xml") is contained in the "bin" folder of the
    @em COMET installation and the OpenMS installation under OpenMS/share/CHEMISTRY/COMET_default_input.xml.
    The latter is loaded by default.
    If you want to use the XML configuration file and @em ignore most of the parameters set via this adapter, use the '-ignore_adapter_param'
    flag. Then, the config given in '-default_config_file' is used exclusively and only '-in', '-out', '-database' and '-COMET_executable' are
    taken from this adapter.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_COMETAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_COMETAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPCOMETAdapter :
  public TOPPBase
{
public:
  TOPPCOMETAdapter() :
    TOPPBase("COMETAdapter", "Annotates MS/MS spectra using COMET.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {

    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerOutputFile_("pin_out", "<file>", "", "Output file - for percolator input");
    setValidFormats_("pin_out", ListUtils::create<String>("pin"), false);
    registerInputFile_("database", "<file>", "", "FASTA file or pro file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));
    registerInputFile_("comet_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      // COMET compiles as tandem on OSX and tandem.exe on any other platform
      "comet.exe",                     
      "Comet executable of the installation e.g. 'comet.exe'", true, false, ListUtils::create<String>("skipexists"));

    addEmptyLine_();
    //
    // Optional parameters (if '-ignore_adapter_param' is set)
    //

    registerDoubleOption_("peptide_mass_tolerance", "<tolerance>", 10.0, "peptide_mass_tolerance", false);
    registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.3, "Fragment mass error", false);

    registerIntOption_("peptide_mass_units", "<search_enzyme_number>", 2, "0=amu, 1=mmu, 2=ppm", false);
    registerStringOption_("precursor_error_units", "<unit>", "ppm", "Parent monoisotopic mass error units", false);
    registerStringOption_("fragment_error_units", "<unit>", "Da", "Fragment monoisotopic mass error units", false);
    vector<String> valid_strings = ListUtils::create<String>("ppm,Da");

    registerIntOption_("mass_type_parent", "<num>", 1, "0=average masses, 1=monoisotopic masses", false);
    registerIntOption_("mass_type_fragment", "<num>", 1, "0=average masses, 1=monoisotopic masses", false);
    registerIntOption_("precursor_tolerance_type", "<num>", 0, "0=average masses, 1=monoisotopic masses", false);
    registerIntOption_("isotope_error", "<num>", 0, "0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)", false);

    registerIntOption_("num_enzyme_termini", "<num>", 2, "1 (semi-digested), 2 (fully digested, default), 8 C-term unspecific , 9 N-term unspecific", false);
    registerIntOption_("allowed_missed_cleavages", "<num>", 1, "Number of possible cleavage sites missed by the enzyme, maximum value is 5; for enzyme search", false);

    registerStringOption_("decoy_prefix", "<unit>", "rev_", "decoy entries are denoted by this string which is pre-pended to each protein accession", false);

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
     registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    addEmptyLine_();
    
    vector<String> all_enzymes;
    EnzymesDB::getInstance()->getAllXTandemNames(all_enzymes);
    registerStringOption_("cleavage_site", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("cleavage_site", all_enzymes);
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
    if (tmp_dir.empty()) return; // no temp. dir. created

    if (debug_level_ >= 2)
    {
      writeDebug_("Keeping temporary files in directory '" + tmp_dir + "'. Set debug level to 1 or lower to remove them.", 2);
    }
    else
    {
      if (debug_level_ == 1) writeDebug_("Deleting temporary directory '" + tmp_dir + "'. Set debug level to 2 or higher to keep it.", 1);
      File::removeDirRecursively(tmp_dir);
    }
  }

  void createParamFile_(ostream& os)
  {
    os << "# comet_version 2016.01 rev. 2\n";               //required as first line in the param file
    os << "# Comet MS/MS search engine parameters file.\n";
    os << "# Everything following the '#' symbol is treated as a comment.\n";
    os << "database_name = " << getStringOption_("database") << "\n";
    os << "decoy_search = " << 0 << "\n"; // 0=no (default), 1=concatenated search, 2=separate search
    os << "num_threads = " << getIntOption_("threads") << "\n";  // 0=poll CPU to set num threads; else specify num threads directly (max 64)

    // masses
    os << "peptide_mass_tolerance = " << getDoubleOption_("peptide_mass_tolerance") << "\n";
    os << "peptide_mass_units = " << getIntOption_("peptide_mass_units") << "\n";                  // 0=amu, 1=mmu, 2=ppm
    os << "mass_type_parent = " << getIntOption_("mass_type_parent") << "\n";                    // 0=average masses, 1=monoisotopic masses
    os << "mass_type_fragment = " << getIntOption_("mass_type_fragment") << "\n";                  // 0=average masses, 1=monoisotopic masses
    os << "precursor_tolerance_type = " << getIntOption_("precursor_tolerance_type") << "\n";            // 0=MH+ (default), 1=precursor m/z; only valid for amu/mmu tolerances
    os << "isotope_error = " << getIntOption_("isotope_error") << "\n";                      // 0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)

    // search enzyme

    // TODO: complete
    map<String, Size> map_oms2comet;
    map_oms2comet["Trypsin"] = 1;
    map_oms2comet["Arg-C"] = 5;
    map_oms2comet["Asp-N"] = 6;
    map_oms2comet["Chymotrypsin"] = 10;
    map_oms2comet["CNBr"] = 7;
    map_oms2comet["Lys-C"] = 3;
    map_oms2comet["PepsinA"] = 9;
    map_oms2comet["Trypsin/P"] = 2;
    map_oms2comet["no cleavage"] = 0;   

    String enzyme_name = getStringOption_("cleavage_site");
    Size enzyme_number = 1;
    if (map_oms2comet.find(enzyme_name) != map_oms2comet.end())
    {
      enzyme_number = map_oms2comet.at(enzyme_name);
    }
    else
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Enzyme not supported. " + enzyme_name);
    }

    os << "search_enzyme_number = " << enzyme_number << "\n";                // choose from list at end of this params file
    os << "num_enzyme_termini = " << getIntOption_("num_enzyme_termini") << "\n";                  // 1 (semi-digested), 2 (fully digested, default), 8 C-term unspecific , 9 N-term unspecific
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

    Size var_mod_index = 1;

    // write out user specified modifications
    for (; var_mod_index <= variable_modifications.size(); ++var_mod_index)
    {
      const ResidueModification mod = variable_modifications[var_mod_index];
      double mass = mod.getDiffMonoMass();
      String residues = mod.getOrigin();  // TODO: check if origin contains C-term string or similar. Should not be passed to commet as residue string
      String variable = "0";
      String max_mods_per_peptide = "3";
      String term_distance = "-1";
      String nc_term = "0";

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
      else if (mod.getTermSpecificity() == ResidueModification::PROTEIN_N_TERM)
      {
        term_distance = 0;
        nc_term = "0";
      }
      else if (mod.getTermSpecificity() == ResidueModification::PROTEIN_C_TERM)
      {
        term_distance = 0;
        nc_term = "1";
      } 
      String required = "0";

      os << "variable_mod0" << var_mod_index << " = " << mass << " " << residues << " " << variable << " " << max_mods_per_peptide << " " << term_distance << " " << nc_term << " " << required << "\n"; 
    }

    // fill remaining modification slots (if any) in Comet with "no modification"
    for (; var_mod_index <= 9; ++var_mod_index)
    {
      os << "variable_mod0" << var_mod_index << " = " << "0.0 X 0 3 -1 0 0" << "\n";
    }

    os << "max_variable_mods_in_peptide = " << 5 << "\n";
    os << "require_variable_mod = " << 0 << "\n";

    // fragment ions
    // ion trap ms/ms:  1.0005 tolerance, 0.4 offset (mono masses), theoretical_fragment_ions = 1
    // high res ms/ms:    0.02 tolerance, 0.0 offset (mono masses), theoretical_fragment_ions = 0
    os << "fragment_bin_tol = " << 1.0005 << "\n";               // binning to use on fragment ions
    os << "fragment_bin_offset = " << 0.4  << "\n";              // offset position to start the binning (0.0 to 1.0)
    os << "theoretical_fragment_ions = " << 1 << "\n";           // 0=use flanking peaks, 1=M peak only
    os << "use_A_ions = " << 0 << "\n";
    os << "use_B_ions = " << 1 << "\n";
    os << "use_C_ions = " << 0 << "\n";
    os << "use_X_ions = " << 0 << "\n";
    os << "use_Y_ions = " << 1 << "\n";
    os << "use_Z_ions = " << 0 << "\n";
    os << "use_NL_ions = " << 0 << "\n";                         // 0=no, 1=yes to consider NH3/H2O neutral loss peaks

    // output
    os << "output_sqtstream = " << 0 << "\n";                    // 0=no, 1=yes  write sqt to standard output
    os << "output_sqtfile = " << 0 << "\n";                      // 0=no, 1=yes  write sqt file
    os << "output_txtfile = " << 0 << "\n";                     // 0=no, 1=yes  write tab-delimited txt file
    os << "output_pepxmlfile = " << 1 << "\n";                   // 0=no, 1=yes  write pep.xml file

    os << "output_percolatorfile = " << !getStringOption_("pin_out").empty() << "\n";              // 0=no, 1=yes  write Percolator tab-delimited input file
    os << "output_outfiles = " <<  0 << "\n";                    // 0=no, 1=yes  write .out files
    os << "print_expect_score = " << 1 << "\n";                  // 0=no, 1=yes to replace Sp with expect in out & sqt
    os << "num_output_lines = " << 5 << "\n";                    // num peptide results to show
    os << "show_fragment_ions = " << 0 << "\n";                  // 0=no, 1=yes for out files only
    os << "sample_enzyme_number = " << 0 << "\n";                // Sample enzyme which is possibly different than the one applied to the search.

    // mzXML parameters
    os << "scan_range = " << "0 0" << "\n";                        // start and scan scan range to search; 0 as 1st entry ignores parameter
    os << "precursor_charge = " << "0 0" << "\n";                  // precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter
    os << "override_charge = " << 0 << "\n";                     // 0=no, 1=override precursor charge states, 2=ignore precursor charges outside precursor_charge range, 3=see online
    os << "ms_level = " << 2 << "\n";                            // MS level to analyze, valid are levels 2 (default) or 3
    os << "activation_method = " << "ALL" << "\n";                 // activation method; used if activation method set; allowed ALL, CID, ECD, ETD, PQD, HCD, IRMPD

    // misc parameters
    os << "digest_mass_range = " << "600.0 5000.0" << "\n";        // MH+ peptide mass range to analyze
    os << "num_results = " << 100 << "\n";                       // number of search hits to store internally
    os << "skip_researching = " << 1 << "\n";                    // for '.out' file output only, 0=search everything again (default), 1=dont search if .out exists
    os << "max_fragment_charge = " << 3 << "\n";                 // set maximum fragment charge state to analyze (allowed max 5)
    os << "max_precursor_charge = " << 4 << "\n";                // set maximum precursor charge state to analyze (allowed max 9)
    os << "nucleotide_reading_frame = " << 0 << "\n";            // 0=proteinDB, 1-6, 7=forward three, 8=reverse three, 9=all six
    os << "clip_nterm_methionine = " << 0 << "\n";              // 0=leave sequences as-is; 1=also consider sequence w/o N-term methionine
    os << "spectrum_batch_size = " << 0 << "\n";                 // max. // of spectra to search at a time; 0 to search the entire scan range in one loop
    os << "decoy_prefix = " << getStringOption_("decoy_prefix") << "\n";                 // decoy entries are denoted by this string which is pre-pended to each protein accession
    os << "output_suffix = " << "" << "\n";                      // add a suffix to output base names i.e. suffix "-C" generates base-C.pep.xml from base.mzXML input
    os << "mass_offsets = " << "" << "\n";                       // one or more mass offsets to search (values substracted from deconvoluted precursor mass)

    // spectral processing
    os << "minimum_peaks = " << 10 << "\n";                      // required minimum number of peaks in spectrum to search (default 10)
    os << "minimum_intensity = " << 0 << "\n";                   // minimum intensity value to read in
    os << "remove_precursor_peak = " << 0 << "\n";               // 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD)
    os << "remove_precursor_tolerance = " << 1.5 << "\n";        // +- Da tolerance for precursor removal
    os << "clear_mz_range = " << "0.0 0.0" << "\n";                // for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range

    // additional modifications
    os << "add_Cterm_peptide = " << 0.0 << "\n";
    os << "add_Nterm_peptide = " <<  0.0 << "\n";
    os << "add_Cterm_protein = " << 0.0 << "\n";
    os << "add_Nterm_protein = " << 0.0 << "\n";
    os << "add_G_glycine = " << 0.0000 << "\n";                  // added to G - avg.  57.0513, mono.  57.02146
    os << "add_A_alanine = " << 0.0000 << "\n";                  // added to A - avg.  71.0779, mono.  71.03711
    os << "add_S_serine = " << 0.0000 << "\n";                  // added to S - avg.  87.0773, mono.  87.03203
    os << "add_P_proline = " << 0.0000 << "\n";                  // added to P - avg.  97.1152, mono.  97.05276
    os << "add_V_valine = " << 0.0000 << "\n";                   // added to V - avg.  99.1311, mono.  99.06841
    os << "add_T_threonine = " << 0.0000 << "\n";                // added to T - avg. 101.1038, mono. 101.04768
    os << "add_C_cysteine = " << 0.0000 << "\n";             // added to C - avg. 103.1429, mono. 103.00918
    os << "add_L_leucine = " << 0.0000 << "\n";                  // added to L - avg. 113.1576, mono. 113.08406
    os << "add_I_isoleucine = " << 0.0000 << "\n";               // added to I - avg. 113.1576, mono. 113.08406
    os << "add_N_asparagine = " << 0.0000 << "\n";               // added to N - avg. 114.1026, mono. 114.04293
    os << "add_D_aspartic_acid = " << 0.0000 << "\n";            // added to D - avg. 115.0874, mono. 115.02694
    os << "add_Q_glutamine = " << 0.0000 << "\n";                // added to Q - avg. 128.1292, mono. 128.05858
    os << "add_K_lysine = " << 0.0000 << "\n";                   // added to K - avg. 128.1723, mono. 128.09496
    os << "add_E_glutamic_acid =" <<  0.0000 << "\n";            // added to E - avg. 129.1140, mono. 129.04259
    os << "add_M_methionine =" <<  0.0000 << "\n";               // added to M - avg. 131.1961, mono. 131.04048
    os << "add_O_ornithine = " << 0.0000 << "\n";                // added to O - avg. 132.1610, mono  132.08988
    os << "add_H_histidine = " << 0.0000 << "\n";                // added to H - avg. 137.1393, mono. 137.05891
    os << "add_F_phenylalanine = " << 0.0000 << "\n";            // added to F - avg. 147.1739, mono. 147.06841
    os << "add_U_selenocysteine = " << 0.0000 << "\n";           // added to U - avg. 150.3079, mono. 150.95363
    os << "add_R_arginine = " << 0.0000 << "\n";                 // added to R - avg. 156.1857, mono. 156.10111
    os << "add_Y_tyrosine = " << 0.0000 << "\n";                 // added to Y - avg. 163.0633, mono. 163.06333
    os << "add_W_tryptophan = " << 0.0000 << "\n";               // added to W - avg. 186.0793, mono. 186.07931
    os << "add_B_user_amino_acid = " << 0.0000 << "\n";          // added to B - avg.   0.0000, mono.   0.00000
    os << "add_J_user_amino_acid = " << 0.0000 << "\n";          // added to J - avg.   0.0000, mono.   0.00000
    os << "add_X_user_amino_acid = " << 0.0000 << "\n";          // added to X - avg.   0.0000, mono.   0.00000
    os << "add_Z_user_amino_acid = " << 0.0000 << "\n";          // added to Z - avg.   0.0000, mono.   0.00000

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
    if (inputfile_name == "")
    {
      writeLog_("No input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    String outputfile_name = getStringOption_("out");
    writeDebug_(String("Output file: ") + outputfile_name, 1);
    if (outputfile_name == "")
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
    const String tmp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/").toQString());
    writeDebug_("Creating temporary directory '" + tmp_dir + "'", 1);
    QDir d;
    d.mkpath(tmp_dir.toQString());

    String tmp_file = tmp_dir + "param.txt";
    String tmp_pepxml = File::removeExtension(inputfile_name) + ".pep.xml";
    String tmp_pin = File::removeExtension(inputfile_name) + ".pin";

    ofstream os(tmp_file);
    createParamFile_(os);
    os.close();

    PeakMap exp;
    MzMLFile mzml_file;
    mzml_file.getOptions().addMSLevel(2); // only load msLevel 2
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

    String param = "-P" + tmp_file;
    QStringList process_params;
    process_params << param.toQString() << inputfile_name.toQString();
    //qDebug() << process_params;

    String comet_executable = getStringOption_("comet_executable");
    //int status = QProcess::execute(comet_executable.toQString(), QStringList(inputfile_name.toQString())); // does automatic escaping etc...
    int status = QProcess::execute(comet_executable.toQString(),process_params); // does automatic escaping etc...
    if (status != 0)
    {
      writeLog_("Comet problem. Aborting! Calling command was: '" + comet_executable + " \"" + inputfile_name + "\"'.\nDoes the Comet executable exist?");
      // clean temporary files
      if (this->debug_level_ < 2)
      {
        removeTempDir_(tmp_dir);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << tmp_dir << "'" << std::endl;
      }
      else
      {
        LOG_WARN << "Keeping the temporary files at '" << tmp_dir << "'. Set debug level to <2 to remove them." << std::endl;
      }
      //return EXTERNAL_PROGRAM_ERROR;
    }

    // read the pep.xml output of COMET and write it to idXML

    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;

    writeDebug_("write PepXMLFile", 1);
    PepXMLFile().load(tmp_pepxml, protein_identifications, peptide_identifications, inputfile_name);

    if (this->debug_level_ == 0)
    {
      File::remove(tmp_pepxml);
      LOG_WARN << "Set debug level to >0 to keep the temporary pep.xml and pin files at '" << tmp_pepxml << "'" << std::endl;
    }
    else
    {
      LOG_WARN << "Keeping the temporary files at '" << tmp_pepxml << "'. Set debug level to 0 to remove them." << std::endl;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    IdXMLFile().store(outputfile_name, protein_identifications, peptide_identifications);

  }

};


int main(int argc, const char** argv)
{
  TOPPCOMETAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
