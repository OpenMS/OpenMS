// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <boost/regex.hpp>

#include <fstream>

using namespace std;

namespace OpenMS
{

  XTandemInfile::XTandemInfile() :
    Internal::XMLFile(),
    fragment_mass_tolerance_(0.3),
    precursor_mass_tolerance_plus_(2.0),
    precursor_mass_tolerance_minus_(2.0),
    fragment_mass_error_unit_(XTandemInfile::DALTONS),
    precursor_mass_error_unit_(XTandemInfile::DALTONS),
    fragment_mass_type_(XTandemInfile::MONOISOTOPIC),
    precursor_mass_type_(XTandemInfile::MONOISOTOPIC),
    max_precursor_charge_(4),
    precursor_lower_mz_(500.0),
    fragment_lower_mz_(200.0),
    number_of_threads_(1),
    modifications_(),
    input_filename_(""),
    output_filename_(""),
    cleavage_site_("[KR]|{P}"),
    semi_cleavage_(false),
    allow_isotope_error_(false),
    number_of_missed_cleavages_(1),
    default_parameters_file_(""),
    output_results_("valid"),
    max_valid_evalue_(0.01),
    force_default_mods_(false)
  {
  }

  XTandemInfile::~XTandemInfile() = default;

  void XTandemInfile::write(const String& filename, bool ignore_member_parameters, bool force_default_mods)
  {
    if (!File::writable(filename))
    {
      throw (Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename));
    }
    force_default_mods_ = force_default_mods;
    ofstream os(filename.c_str());
    writeTo_(os, ignore_member_parameters);
    return;
  }

  String XTandemInfile::convertModificationSet_(const set<ModificationDefinition>& mods, map<String, double>& affected_origins) const
  {
    // check if both "Glu->pyro-Glu (N-term E)" and "Gln->pyro-Glu (N-term Q)"
    // are specified:
    bool has_pyroglu_e = false, has_pyroglu_q = false;
    for (set<ModificationDefinition>::const_iterator it = mods.begin();
         it != mods.end(); ++it)
    {
      if (it->getModificationName() == "Glu->pyro-Glu (N-term E)")
      {
        has_pyroglu_e = true;
      }
      else if (it->getModificationName() == "Gln->pyro-Glu (N-term Q)")
      {
        has_pyroglu_q = true;
      }
      if (has_pyroglu_e && has_pyroglu_q)
      {
        break;
      }
    }

    map<String, double> origin_set;
    StringList xtandem_mods;
    for (set<ModificationDefinition>::const_iterator it = mods.begin();
         it != mods.end(); ++it)
    {
      if (!force_default_mods_ &&
          // @TODO: change Acetyl spec. to "protein N-term" once it's supported
          ((it->getModificationName() == "Acetyl (N-term)") ||
           // for the pyro-Glus, only skip if both are present:
           ((it->getModificationName() == "Gln->pyro-Glu (N-term Q)") &&
            has_pyroglu_e) ||
           ((it->getModificationName() == "Glu->pyro-Glu (N-term E)") &&
            has_pyroglu_q)))
      {
        continue;
      }

      double mod_mass = it->getModification().getDiffMonoMass();

      String orig = it->getModification().getOrigin();
      ResidueModification::TermSpecificity ts = it->getModification().getTermSpecificity();
      if ((ts != ResidueModification::ANYWHERE) && !orig.empty())
      {
        OPENMS_LOG_WARN << "Warning: X! Tandem doesn't support modifications with both residue and terminal specificity. Using only terminal specificity for modification '" << it->getModificationName() << "'." << endl;
      }

      if (ts == ResidueModification::C_TERM)
      {
        orig = "]";
      }
      else if (ts == ResidueModification::N_TERM)
      {
        orig = "[";
      }
      // check double usage
      if (origin_set.find(orig) != origin_set.end())
      {
        OPENMS_LOG_WARN << "X! Tandem config file: Duplicate modification assignment to origin '" << orig << "'. "
                 << "X! Tandem will ignore the first modification '" << origin_set.find(orig)->second << "'!\n";
      }
      // check if already used before (i.e. we are currently looking at variable mods)
      if (affected_origins.find(orig) != affected_origins.end())
      {
        OPENMS_LOG_INFO << "X! Tandem config file: Fixed modification and variable modification to origin '" << orig << "' detected. "
                 << "Using corrected mass of " << mod_mass - affected_origins.find(orig)->second << " instead of " << mod_mass << ".\n";
        mod_mass -= affected_origins.find(orig)->second;
      }
      // insert the (corrected) value
      origin_set.insert(make_pair(orig, mod_mass));

      String mod_string;
      if (mod_mass >= 0)
      {
        mod_string = String("+") + String(mod_mass); // prepend a "+"
      }
      else
      {
        mod_string = String(mod_mass); // the '-' is implicit
      }
      mod_string += "@" + orig;
      xtandem_mods.push_back(mod_string);
    }

    // copy now; above we need an independent set, in case 'affected_origins' was non-empty
    affected_origins = origin_set;

    return ListUtils::concatenate(xtandem_mods, ",");
  }

  void XTandemInfile::writeTo_(ostream& os, bool ignore_member_parameters)
  {
    os << "<?xml version=\"1.0\"?>" << "\n"
       << R"(<?xml-stylesheet type="text/xsl" href="tandem-input-style.xsl"?>)" << "\n"
       << "<bioml>" << "\n";

    writeNote_(os, "spectrum, path", input_filename_);
    writeNote_(os, "output, path", output_filename_);
    writeNote_(os, "list path, taxonomy information", taxonomy_file_); // contains the FASTA database
    if (!default_parameters_file_.empty())
    {
      writeNote_(os, "list path, default parameters", default_parameters_file_);
    }
    // these are needed for finding and parsing the results:
    writeNote_(os, "output, path hashing", false);
    writeNote_(os, "output, proteins", true);
    writeNote_(os, "output, spectra", true);
    writeNote_(os, "output, sort results by", "spectrum");
    // required by Percolator to recognize output file
    // (see https://github.com/percolator/percolator/issues/180):
    writeNote_(os, "output, xsl path", "tandem-style.xsl");
    // to help diagnose problems:
    writeNote_(os, "output, parameters", true);

    if (!ignore_member_parameters)
    {
      //////////////// spectrum parameters
      //<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
      writeNote_(os, "spectrum, fragment monoisotopic mass error", String(fragment_mass_tolerance_));
      //<note type="input" label="spectrum, parent monoisotopic mass error plus">100</note>
      writeNote_(os, "spectrum, parent monoisotopic mass error plus", String(precursor_mass_tolerance_plus_));
      //<note type="input" label="spectrum, parent monoisotopic mass error minus">100</note>
      writeNote_(os, "spectrum, parent monoisotopic mass error minus", String(precursor_mass_tolerance_minus_));
      //<note type="input" label="spectrum, parent monoisotopic mass isotope error">yes</note>
      String allow = allow_isotope_error_ ? "yes" : "no";
      writeNote_(os, "spectrum, parent monoisotopic mass isotope error", allow);
      //<note type="input" label="spectrum, fragment monoisotopic mass error units">Daltons</note>
      //<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
      if (fragment_mass_error_unit_ == XTandemInfile::DALTONS)
      {
        writeNote_(os, "spectrum, fragment monoisotopic mass error units", "Daltons");
      }
      else
      {
        writeNote_(os, "spectrum, fragment monoisotopic mass error units", "ppm");
      }

      //<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>
      //<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
      if (precursor_mass_error_unit_ == XTandemInfile::PPM)
      {
        writeNote_(os, "spectrum, parent monoisotopic mass error units", "ppm");
      }
      else
      {
        writeNote_(os, "spectrum, parent monoisotopic mass error units", "Daltons");
      }

      //<note type="input" label="spectrum, fragment mass type">monoisotopic</note>
      //<note>values are monoisotopic|average </note>
      if (fragment_mass_type_ == XTandemInfile::MONOISOTOPIC)
      {
        writeNote_(os, "spectrum, fragment mass type", "monoisotopic"); // default
      }
      else
      {
        writeNote_(os, "spectrum, fragment mass type", "average");
      }
      ////////////////////////////////////////////////////////////////////////////////


      //////////////// spectrum conditioning parameters
      //<note type="input" label="spectrum, dynamic range">100.0</note>
      //<note>The peaks read in are normalized so that the most intense peak
      //is set to the dynamic range value. All peaks with values of less that
      //1, using this normalization, are not used. This normalization has the
      //overall effect of setting a threshold value for peak intensities.</note>
      //writeNote_(os, "spectrum, dynamic range", String(dynamic_range_);

      //<note type="input" label="spectrum, total peaks">50</note>
      //<note>If this value is 0, it is ignored. If it is greater than zero (lets say 50),
      //then the number of peaks in the spectrum with be limited to the 50 most intense
      //peaks in the spectrum. X! tandem does not do any peak finding: it only
      //limits the peaks used by this parameter, and the dynamic range parameter.</note>
      //writeNote_(os, "spectrum, total peaks", String(total_number_peaks_);

      //<note type="input" label="spectrum, maximum parent charge">4</note>
      writeNote_(os, "spectrum, maximum parent charge", String(max_precursor_charge_));

      // <note type="input" label="spectrum, use noise suppression">yes</note>
      //writeNote_(os, "spectrum, use noise suppression", noise_suppression_);

      //<note type="input" label="spectrum, minimum parent m+h">500.0</note>
      //writeNote_(os, "spectrum, minimum parent m+h", String(precursor_lower_mz_));

      //<note type="input" label="spectrum, minimum fragment mz">150.0</note>
      //writeNote_(os, "spectrum, minimum fragment mz", String(fragment_lower_mz_));

      //<note type="input" label="spectrum, minimum peaks">15</note>
      //writeNote_(os, "spectrum, minimum peaks", String(min_number_peaks_));

      //<note type="input" label="spectrum, threads">1</note>
      writeNote_(os, "spectrum, threads", String(number_of_threads_));

      //<note type="input" label="spectrum, sequence batch size">1000</note>
      //writeNote_(os, "spectrum, sequence batch size", String(batch_size_));
      ////////////////////////////////////////////////////////////////////////////////


      //////////////// protein parameters
      //<note type="input" label="protein, taxon">other mammals</note>
      //<note>This value is interpreted using the information in taxonomy.xml.</note>
      writeNote_(os, "protein, taxon", taxon_);

      //<note type="input" label="protein, cleavage site">[RK]|{P}</note>
      //<note>this setting corresponds to the enzyme trypsin. The first characters
      //in brackets represent residues N-terminal to the bond - the '|' pipe -
      //and the second set of characters represent residues C-terminal to the
      //bond. The characters must be in square brackets (denoting that only
      //these residues are allowed for a cleavage) or french brackets (denoting
      //that these residues cannot be in that position). Use UPPERCASE characters.
      //To denote cleavage at any residue, use [X]|[X] and reset the
      //scoring, maximum missed cleavage site parameter (see below) to something like 50.
      //</note>
      writeNote_(os, "protein, cleavage site", cleavage_site_);
      
      //////////////// semi cleavage parameter
      //<note type="input" label="protein, cleavage semi">yes</note>
      writeNote_(os, "protein, cleavage semi", semi_cleavage_);

      //<note type="input" label="protein, modified residue mass file"></note>
      //writeNote_(os, "protein, modified residue mass file", modified_residue_mass_file_);

      //<note type="input" label="protein, cleavage C-terminal mass change">+17.002735</note>
      //writeNote_(os, "protein, cleavage C-terminal mass change", String(cleavage_c_term_mass_change_));

      //<note type="input" label="protein, cleavage N-terminal mass change">+1.007825</note>
      //writeNote_(os, "protein, cleavage N-terminal mass change", String(cleavage_n_term_mass_change_));

      //<note type="input" label="protein, N-terminal residue modification mass">0.0</note>
      //writeNote_(os, "protein, N-terminal residue modification mass", String(protein_n_term_mod_mass_));

      //<note type="input" label="protein, C-terminal residue modification mass">0.0</note>
      //writeNote_(os, "protein, C-terminal residue modification mass", String(protein_c_term_mod_mass_));

      //<note type="input" label="protein, homolog management">no</note>
      //<note>if yes, an upper limit is set on the number of homologues kept for a particular spectrum</note>
      //writeNote_(os, "protein, homolog management", protein_homolog_management_);

      // special cases for default (N-terminal) modifications:
      set<String> var_mods = modifications_.getVariableModificationNames();
      // Ron Beavis: "If a variable modification is set for the peptide N-terminus, the 'quick acetyl' and 'quick pyrolidone' are turned off so that they don't interfere with the specified variable modification." -> check for that
      boost::regex re(" \\(N-term( .)?\\)$");
      for (set<String>::iterator vm_it = var_mods.begin();
           vm_it != var_mods.end(); ++vm_it)
      {
        if (boost::regex_search(*vm_it, re) && (*vm_it != "Acetyl (N-term)") &&
            (*vm_it != "Gln->pyro-Glu (N-term Q)") &&
            (*vm_it != "Glu->pyro-Glu (N-term E)"))
        {
          force_default_mods_ = true;
        }
      }

      if (!force_default_mods_ &&
          (var_mods.find("Gln->pyro-Glu (N-term Q)") != var_mods.end()) &&
          (var_mods.find("Glu->pyro-Glu (N-term E)") != var_mods.end()))
      {
        writeNote_(os, "protein, quick pyrolidone", true);
        OPENMS_LOG_INFO << "Modifications 'Gln->pyro-Glu (N-term Q)' and 'Glu->pyro-Glu (N-term E)' are handled implicitly by the X! Tandem option 'protein, quick pyrolidone'. Set the 'force' flag in XTandemAdapter to force explicit inclusion of these modifications." << endl;
      }

      // special case for "Acetyl (N-term)" modification:
      if (!force_default_mods_ &&
          (var_mods.find("Acetyl (N-term)") != var_mods.end()))
      {
        writeNote_(os, "protein, quick acetyl", true);
        OPENMS_LOG_INFO << "Modification 'Acetyl (N-term)' is handled implicitly by the X! Tandem option 'protein, quick acetyl'. Set the 'force' flag in XTandemAdapter to force explicit inclusion of this modification." << endl;
      }

      ////////////////////////////////////////////////////////////////////////////////


      //////////////// residue modification parameters
      //<note type="input" label="residue, modification mass">57.022@C</note>
      //<note>The format of this parameter is m@X, where m is the modification
      //mass in Daltons and X is the appropriate residue to modify. Lists of
      //modifications are separated by commas. For example, to modify M and C
      //with the addition of 16.0 Daltons, the parameter line would be
      //+16.0@M,+16.0@C
      //Positive and negative values are allowed.
      //</note>

      map<String, double> affected_origins;
      writeNote_(os, "residue, modification mass", convertModificationSet_(modifications_.getFixedModifications(), affected_origins));

      //<note type="input" label="residue, potential modification mass"></note>
      //<note>The format of this parameter is the same as the format
      //for residue, modification mass (see above).</note>
      writeNote_(os, "residue, potential modification mass", convertModificationSet_(modifications_.getVariableModifications(), affected_origins));

      //<note type="input" label="residue, potential modification motif"></note>
      //<note>The format of this parameter is similar to residue, modification mass,
      //with the addition of a modified PROSITE notation sequence motif specification.
      //For example, a value of 80@[ST!]PX[KR] indicates a modification
      //of either S or T when followed by P, and residue and the a K or an R.
      //A value of 204@N!{P}[ST]{P} indicates a modification of N by 204, if it
      //is NOT followed by a P, then either an S or a T, NOT followed by a P.
      //Positive and negative values are allowed.
      //</note>
      //    writeNote_(os, "residue, potential modification motif", variable_modification_motif_);
      ////////////////////////////////////////////////////////////////////////////////


      //////////////// model refinement parameters
      //<note type="input" label="refine">yes</note>
      //writeNote_(os, "refine", refine_);
      //<note type="input" label="refine, modification mass"></note>
      //writeNote_(os, "refine, modification mass", String(refine_mod_mass_));
      //<note type="input" label="refine, sequence path"></note>
      //writeNote_(os, "refine, sequence path", refine_sequence_path_);
      //<note type="input" label="refine, tic percent">20</note>
      //writeNote_(os, "refine, tic percent", String(refine_tic_percent_));
      //<note type="input" label="refine, spectrum synthesis">yes</note>
      //writeNote_(os, "refine, spectrum synthesis", refine_spectrum_synthesis_);
      //<note type="input" label="refine, maximum valid expectation value">0.1</note>
      //writeNote_(os, "refine, maximum valid expectation value", String(refine_max_valid_evalue_));
      //<note type="input" label="refine, potential N-terminus modifications">+42.010565@[</note>
      //writeNote_(os, "refine, potential N-terminus modifications", refine_variable_n_term_mods_);
      //<note type="input" label="refine, potential C-terminus modifications"></note>
      //writeNote_(os, "refine, potential C-terminus modifications", refine_variable_c_term_mods_);
      //<note type="input" label="refine, unanticipated cleavage">yes</note>
      //writeNote_(os, "refine, unanticipated cleavage", refine_unanticipated_cleavage_);
      //<note type="input" label="refine, potential modification mass"></note>
      //writeNote_(os, "refine, potential modification mass", String(variable_mod_mass_));
      //<note type="input" label="refine, point mutations">no</note>
      //writeNote_(os, "refine, point mutations", refine_point_mutations_);
      //<note type="input" label="refine, use potential modifications for full refinement">no</note>
      //writeNote_(os, "refine, use potential modifications for full refinement", use_var_mod_for_full_refinement_);
      //<note type="input" label="refine, potential modification motif"></note>
      //<note>The format of this parameter is similar to residue, modification mass,
      //with the addition of a modified PROSITE notation sequence motif specification.
      //For example, a value of 80@[ST!]PX[KR] indicates a modification
      //of either S or T when followed by P, and residue and the a K or an R.
      //A value of 204@N!{P}[ST]{P} indicates a modification of N by 204, if it
      //is NOT followed by a P, then either an S or a T, NOT followed by a P.
      //Positive and negative values are allowed.
      //</note>

      //writeNote_(os, "refine, potential modification motif", refine_var_mod_motif_);
      ////////////////////////////////////////////////////////////////////////////////


      //////////////// scoring parameters
      //<note type="input" label="scoring, minimum ion count">4</note>
      //writeNote_(os, "scoring, minimum ion count", String(scoring_min_ion_count_));
      //<note type="input" label="scoring, maximum missed cleavage sites">1</note>
      writeNote_(os, "scoring, maximum missed cleavage sites", String(number_of_missed_cleavages_));
      //<note type="input" label="scoring, x ions">no</note>
      //writeNote_(os, "scoring, x ions", score_x_ions_);
      //<note type="input" label="scoring, y ions">yes</note>
      //writeNote_(os, "scoring, y ions", score_y_ions_);
      //<note type="input" label="scoring, z ions">no</note>
      //writeNote_(os, "scoring, z ions", score_z_ions_);
      //<note type="input" label="scoring, a ions">no</note>
      //writeNote_(os, "scoring, a ions", score_a_ions_);
      //<note type="input" label="scoring, b ions">yes</note>
      //writeNote_(os, "scoring, b ions", score_b_ions_);
      //<note type="input" label="scoring, c ions">no</note>
      //writeNote_(os, "scoring, c ions", score_c_ions_);
      //<note type="input" label="scoring, cyclic permutation">no</note>
      //<note>if yes, cyclic peptide sequence permutation is used to pad the scoring histograms</note>
      //writeNote_(os, "scoring, cyclic permutation", scoring_cyclic_permutation_);
      //<note type="input" label="scoring, include reverse">no</note>
      //<note>if yes, then reversed sequences are searched at the same time as forward sequences</note>
      //writeNote_(os, "scoring, include reverse", scoring_include_reverse_);
      ////////////////////////////////////////////////////////////////////////////////


      //////////////// output parameters
      //<note type="input" label="output, log path"></note>
      //<note type="input" label="output, message">...</note>
      //writeNote_(os, "output, message", String("..."));
      //<note type="input" label="output, one sequence copy">no</note>
      //<note type="input" label="output, sequence path"></note>
      //<note type="input" label="output, path">output.xml</note>
      //writeNote_(os, "output, path", output_filename_);
      //<note type="input" label="output, sort results by">protein</note>
      //<note>values = protein|spectrum (spectrum is the default)</note>
      //<note type="input" label="output, path hashing">yes</note>
      //<note>values = yes|no</note>
      //<note type="input" label="output, xsl path">tandem-style.xsl</note>
      //<note type="input" label="output, parameters">yes</note>
      //<note>values = yes|no</note>
      //<note type="input" label="output, performance">yes</note>
      //<note>values = yes|no</note>
      //<note type="input" label="output, spectra">yes</note>
      //<note>values = yes|no</note>
      //<note type="input" label="output, histograms">yes</note>
      //<note>values = yes|no</note>
      //<note type="input" label="output, proteins">yes</note>
      //<note>values = yes|no</note>
      //<note type="input" label="output, sequences">yes</note>
      //<note>values = yes|no</note>
      //<note type="input" label="output, one sequence copy">no</note>
      //<note>values = yes|no, set to yes to produce only one copy of each protein sequence in the output xml</note>
      //<note type="input" label="output, results">valid</note>
      //<note>values = all|valid|stochastic</note>
      writeNote_(os, "output, results", output_results_);
 
      //<note type="input" label="output, maximum valid expectation value">0.1</note>
      writeNote_(os, "output, maximum valid expectation value", String(max_valid_evalue_));

      //<note>value is used in the valid|stochastic setting of output, results</note>
      //<note type="input" label="output, histogram column width">30</note>
      //<note>values any integer greater than 0. Setting this to '1' makes cutting and pasting histograms
      //into spread sheet programs easier.</note>
      //<note type="description">ADDITIONAL EXPLANATIONS</note>
      //<note type="description">Each one of the parameters for X! tandem is entered as a labeled note
      //  node. In the current version of X!, keep those note nodes
      //  on a single line.
      //</note>
      //<note type="description">The presence of the type 'input' is necessary if a note is to be considered
      //  an input parameter.
      //</note>
      //<note type="description">Any of the parameters that are paths to files may require alteration for a
      //  particular installation. Full path names usually cause the least trouble,
      //  but there is no reason not to use relative path names, if that is the
      //  most convenient.
      //</note>
      //<note type="description">Any parameter values set in the 'list path, default parameters' file are
      //  reset by entries in the normal input file, if they are present. Otherwise,
      //  the default set is used.
      //</note>
      //<note type="description">The 'list path, taxonomy information' file must exist.
      //</note>
      //<note type="description">The directory containing the 'output, path' file must exist: it will not be created.
      //</note>
      //<note type="description">The 'output, xsl path' is optional: it is only of use if a good XSLT style sheet exists.
      //</note>
      ////////////////////////////////////////////////////////////////////////////////
    }

    os << "</bioml>\n";
  }

  void XTandemInfile::writeNote_(ostream& os, const String& label, const String& value)
  {
    os << "\t<note type=\"input\" label=\"" << label << "\">" << value << "</note>\n";
  }

  void XTandemInfile::writeNote_(ostream& os, const String& label, const char* value)
  {
    String val(value);
    writeNote_(os, label, val);
  }

  void XTandemInfile::writeNote_(ostream& os, const String& label, bool value)
  {
    String val = value ? "yes" : "no";
    writeNote_(os, label, val);
  }

  void XTandemInfile::setOutputFilename(const String& filename)
  {
    output_filename_ = filename;
  }

  const String& XTandemInfile::getOutputFilename() const
  {
    return output_filename_;
  }

  void XTandemInfile::setInputFilename(const String& filename)
  {
    input_filename_ = filename;
  }

  const String& XTandemInfile::getInputFilename() const
  {
    return input_filename_;
  }

  void XTandemInfile::setTaxonomyFilename(const String& filename)
  {
    taxonomy_file_ = filename;
  }

  const String& XTandemInfile::getTaxonomyFilename() const
  {
    return taxonomy_file_;
  }

  void XTandemInfile::setDefaultParametersFilename(const String& filename)
  {
    default_parameters_file_ = filename;
  }

  const String& XTandemInfile::getDefaultParametersFilename() const
  {
    return default_parameters_file_;
  }

  void XTandemInfile::setModifications(const ModificationDefinitionsSet& mods)
  {
    modifications_ = mods;
  }

  const ModificationDefinitionsSet& XTandemInfile::getModifications() const
  {
    return modifications_;
  }

  void XTandemInfile::setTaxon(const String& taxon)
  {
    taxon_ = taxon;
  }

  const String& XTandemInfile::getTaxon() const
  {
    return taxon_;
  }

  void XTandemInfile::setPrecursorMassTolerancePlus(double tolerance)
  {
    precursor_mass_tolerance_plus_ = tolerance;
  }

  double XTandemInfile::getPrecursorMassTolerancePlus() const
  {
    return precursor_mass_tolerance_plus_;
  }

  void XTandemInfile::setPrecursorMassToleranceMinus(double tolerance)
  {
    precursor_mass_tolerance_minus_ = tolerance;
  }

  double XTandemInfile::getPrecursorMassToleranceMinus() const
  {
    return precursor_mass_tolerance_minus_;
  }

  void XTandemInfile::setPrecursorMassErrorUnit(ErrorUnit unit)
  {
    precursor_mass_error_unit_ = unit;
  }

  XTandemInfile::ErrorUnit XTandemInfile::getPrecursorMassErrorUnit() const
  {
    return precursor_mass_error_unit_;
  }

  void XTandemInfile::setFragmentMassErrorUnit(ErrorUnit unit)
  {
    fragment_mass_error_unit_ = unit;
  }

  XTandemInfile::ErrorUnit XTandemInfile::getFragmentMassErrorUnit() const
  {
    return fragment_mass_error_unit_;
  }

  void XTandemInfile::setMaxPrecursorCharge(Int max_charge)
  {
    max_precursor_charge_ = max_charge;
  }

  Int XTandemInfile::getMaxPrecursorCharge() const
  {
    return max_precursor_charge_;
  }

  void XTandemInfile::setFragmentMassTolerance(double tolerance)
  {
    fragment_mass_tolerance_ = tolerance;
  }

  double XTandemInfile::getFragmentMassTolerance() const
  {
    return fragment_mass_tolerance_;
  }

  void XTandemInfile::setNumberOfThreads(UInt num_threads)
  {
    number_of_threads_ = num_threads;
  }

  UInt XTandemInfile::getNumberOfThreads() const
  {
    return number_of_threads_;
  }

  XTandemInfile::MassType XTandemInfile::getPrecursorErrorType() const
  {
    return precursor_mass_type_;
  }

  void XTandemInfile::setPrecursorErrorType(const MassType mass_type)
  {
    precursor_mass_type_ = mass_type;
  }

  void XTandemInfile::setMaxValidEValue(double value)
  {
    max_valid_evalue_ = value;
  }

  double XTandemInfile::getMaxValidEValue() const
  {
    return max_valid_evalue_;
  }

  void XTandemInfile::setNumberOfMissedCleavages(UInt missed_cleavages)
  {
    number_of_missed_cleavages_ = missed_cleavages;
  }

  UInt XTandemInfile::getNumberOfMissedCleavages() const
  {
    return number_of_missed_cleavages_;
  }

  void XTandemInfile::setOutputResults(const String& result)
  {
    if (result == "valid" || result == "all" || result == "stochastic")
    {
      output_results_ = result;
    }
    else
    {
      throw OpenMS::Exception::FailedAPICall(__FILE__, __LINE__, __FUNCTION__, "Invalid result type provided (must be either all, valid or stochastic).: '" + result + "'");
    }
  }

  String XTandemInfile::getOutputResults() const
  {
    return output_results_;
  }

  void XTandemInfile::setSemiCleavage(const bool semi_cleavage)
  {
    semi_cleavage_ = semi_cleavage;
  }

  void XTandemInfile::setAllowIsotopeError(const bool allow_isotope_error)
  {
    allow_isotope_error_ = allow_isotope_error;
  }
  
  void XTandemInfile::setCleavageSite(const String& cleavage_site)
  {
    cleavage_site_ = cleavage_site;
  }

  const String& XTandemInfile::getCleavageSite() const
  {
    return cleavage_site_;
  }

} // namespace OpenMS
