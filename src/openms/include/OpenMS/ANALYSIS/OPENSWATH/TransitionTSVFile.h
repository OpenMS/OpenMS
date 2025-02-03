// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <fstream>

namespace OpenMS
{

  /**
      @brief This class supports reading and writing of OpenSWATH transition
      lists.

      The transition lists can be either comma- or tab-separated plain
      text files (CSV or TSV format).  Modifications should be provided in
      UniMod format<sup>1</sup>, but can also be provided in TPP format. 
      For another file format that stores transitions, see also TransitionPQPFile.

      The following columns are required:

      <table>
        <tr> <td BGCOLOR="#EBEBEB">PrecursorMz*</td> <td>float</td> <td>This describes the precursor ion \a m/z</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ProductMz*</td> <td>float; synonyms: FragmentMz</td> <td>This specifies the product ion \a m/z</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">LibraryIntensity*</td> <td>float; synonyms: RelativeFragmentIntensity</td> <td>This specifies the relative intensity of the fragment ion</td></tr>
        <tr> <td BGCOLOR="#EBEBEB">NormalizedRetentionTime*</td> <td>float; synonyms: RetentionTime, Tr_recalibrated, iRT, RetentionTimeCalculatorScore</td> <td>This specifies the expected retention time (normalized retention time) </td></tr>
      </table>

  For targeted proteomics files, the following additional columns should be provided:
        <table>
          <tr> <td BGCOLOR="#EBEBEB">GeneName**</td> <td>free text; </td><td> Gene name (unique gene identifier)</td></tr>
          <tr> <td BGCOLOR="#EBEBEB">ProteinId**</td> <td>free text; synonyms: ProteinName</td><td> Protein identifier</td></tr>
          <tr> <td BGCOLOR="#EBEBEB">PeptideSequence**</td> <td>free text</td> <td> sequence only (no modifications); synonyms: Sequence, StrippedSequence</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">ModifiedPeptideSequence**</td> <td>free text</td> <td> should contain modifications<sup>1</sup>; synonyms: FullUniModPeptideName, FullPeptideName, ModifiedSequence</td>  </tr>
          <tr> <td BGCOLOR="#EBEBEB">PrecursorCharge**</td> <td>integer; synonyms: Charge</td> <td> contains the charge of the precursor ion</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">ProductCharge**</td> <td>integer; synonyms: FragmentCharge</td> <td> contains the fragment charge</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">FragmentType</td> <td>free text</td> <td> contains the type of the fragment, e.g. "b" or "y"</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">FragmentSeriesNumber</td> <td>integer; synonyms: FragmentNumber</td> <td> e.g. for y7 use "7" here</td> </tr>
        </table>

  OpenSWATH uses grouped transitions to detect candidate analyte signals. These groups are by default generated based on the input, but can also be manually specified:

        <table>
          <tr> <td BGCOLOR="#EBEBEB">TransitionGroupId**</td> <td>free text; synomys: TransitionGroupName, transition_group_id</td><td> designates the transition group [e.g. peptide] to which this transition belongs</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">TransitionId**</td> <td>free text; synonyms: TransitionName, transition_name </td> <td> needs to be unique for each transition [in this file]</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">Decoy</td> <td>1: decoy, 0: target; synomys: decoy, IsDecoy </td> <td> determines whether the transition is a decoy transition or not</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">PeptideGroupLabel</td> <td>free text </td> <td> designates to which peptide label group (as defined in MS:1000893) the peptide belongs to<sup>2</sup></td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">DetectingTransition</td> <td> 0 or 1; synonyms: detecting_transition</td> <td>1: use transition to detect peak group, 0: don't use transition for detection</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">IdentifyingTransition</td> <td> 0 or 1; synonyms: identifying_transition</td> <td>1: use transition for peptidoform inference in the <a href="http://openswath.org/en/latest/docs/ipf.html">IPF Workflow</a>, 0: don't use transition for identification</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">QuantifyingTransition</td> <td> 0 or 1; synonyms: quantifying_transition</td> <td>1: use transition to quantify peak group, 0: don't use transition for quantification</td> </tr>
        </table>

  Optionally, the following columns can be specified but they are not actively used by OpenSWATH:
        <table>
          <tr>  <td BGCOLOR="#EBEBEB">CollisionEnergy </td><td>float; synonyms: CE</td><td>Collision energy </td></tr>
          <tr>  <td BGCOLOR="#EBEBEB">Annotation </td><td>free text</td><td>Transition-level annotation, e.g. y7</td> </tr>
          <tr>  <td BGCOLOR="#EBEBEB">UniprotId </td><td>free text; synonyms: UniprotID</td> <td>A Uniprot identifier </td></tr>
          <tr>  <td BGCOLOR="#EBEBEB">LabelType </td><td>free text</td><td>optional description of which label was used, e.g. heavy or light</td> </tr>
        </table>

  For targeted metabolomics files, the following fields are also supported:

        <table>
          <tr> <td BGCOLOR="#EBEBEB">CompoundName**</td> <td>free text; synonyms: CompoundId</td><td>Should be unique for the analyte, if present the file will be interpreted as a metabolomics file </td></tr>
          <tr> <td BGCOLOR="#EBEBEB">SMILES</td><td>free text</td><td>SMILES identifier of the compound</td></tr>
          <tr> <td BGCOLOR="#EBEBEB">SumFormula</td><td>free text</td><td>molecular formula of the compound (e.g. H2O)</td></tr>
        </table>

  Fields indicated with * are strictly required to create a TraML file. Fields
  indicated with ** are recommended, but only required for a specific
  application (such as using the transition list for an analysis tool such as
  \ref TOPP_OpenSwathWorkflow) or in a specific context (proteomics or metabolomics).

<p>
Remarks:
</p>
<ul>
  <li>
    1. modifications should be supplied inside the sequence using UniMod
      identifiers or freetext identifiers that are understood by %OpenMS. See also @ref OpenMS::AASequence for more information. For example:
      <ul>
      <li> PEPT(Phosphorylation)IDE(UniMod:27)A ) </li>
      </ul>
  </li>
  <li>
    2. peptide label groups designate groups of peptides that are isotopically
    modified forms of the same peptide species. For example, the heavy and
    light forms of the same peptide will both be assigned the same peptide
    group label. For example:
      <ul>
      <li> PEPTIDEAK -> gets label "PEPTIDEAK_gr1"  </li>
      <li> PEPTIDEAK[+8] -> gets label "PEPTIDEAK_gr1"  </li>
      <li> PEPT(Phosphorylation)IDEAK -> gets label "PEPTIDEAK_gr2"  </li>
      <li> PEPT(Phosphorylation)IDEAK[+8] -> gets label "PEPTIDEAK_gr2"  </li>
      </ul>
  </li>
</ul>
</p>



  @htmlinclude OpenMS_TransitionTSVFile.parameters

  */
  class OPENMS_DLLAPI TransitionTSVFile :
    public ProgressLogger,
    public DefaultParamHandler
  {

protected:
    /**
      @brief Internal structure to represent a transition

      Internal structure to represent a single line from a transition input
      file (one transition).
    */
    struct TSVTransition
    {
      double precursor = -1; ///< Precursor m/z
      double product = -1; ///< Product m/z (fragment ion m/z)
      double rt_calibrated = -1; ///< Normalized RT
      String transition_name = ""; ///< Unique transition name
      double CE = -1; ///< Collision Energy
      double library_intensity = -1; ///< Library intensity of fragment ion (relative)
      String group_id = ""; ///< Transition group identifier (grouping transitions of the same analyte)
      bool decoy = false; ///< Whether the transition is a decoy transition
      String PeptideSequence; ///< Peptide sequence (only AA sequence)
      std::vector<String> ProteinName; ///< List of protein identifiers
      String GeneName; ///< Gene identifier
      String Annotation; ///< Fragment ion annotation
      String FullPeptideName; ///< Full peptide sequence with UniMod modifications
      String CompoundName; ///< Compound name (for metabolomics)
      String SMILES; ///< SMILES identifier (for metabolomics)
      String SumFormula; ///< Molecular formula (for metabolomics)
      String Adducts; ///< Adducts (for metabolomics)
      String precursor_charge; ///< Precursor charge state
      String peptide_group_label; ///< Peptide group identifier (grouping isotopically labelled peptides)
      String label_type; ///< Type of label that was used (e.g. "heavy" or "light")
      String fragment_charge = "NA"; ///< Fragment ion charge state
      int fragment_nr = -1; ///< Fragment number (e.g. "7" for a y7 ion)
      double fragment_mzdelta = -1; ///< Fragment m/z delta to theoretical ion
      double drift_time = -1; ///< Ion mobility drift time
      int fragment_modification = 0; ///< Fragment modification
      String fragment_type; ///< Fragment type (e.g. "y" for a y7 ion)
      std::vector<String> uniprot_id; ///< List of UniProt identifiers of associated proteins
      bool detecting_transition = true; ///< Whether to use transition to detect peak group,
      bool identifying_transition = false; ///< Whether to use transition for peptidoform inference using IPF
      bool quantifying_transition = true; ///< Whether to use transition to quantify peak group
      std::vector<String> peptidoforms; ///< List of peptidoforms

      /// Whether the transition represents a peptide (by convention, if the
      /// (metabolic) compound name field is empty, it is a peptide.)
      bool isPeptide() const
      {
        return CompoundName.empty() || CompoundName == "NA";
      }
    };

    /** @name  Conversion functions from TSVTransition objects to OpenMS datastructures
     *
     * These functions convert the relevant data from a TSVTransition to the
     * datastructures used by the TraML handler or by the internal LightTargetedExperiment.
     *
   */
    //@{

    /** @brief Convert a list of TSVTransition to a TargetedExperiment
     *
     * Converts the list (read from csv/mrm) file into a object model using the
     * TargetedExperiment with proper hierarchical structure from Transition to
     * Peptide to Protein.
     *
    */
    void TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenMS::TargetedExperiment& exp);

    /** @brief Convert a list of TSVTransition to a LightTargetedExperiment
     *
     * Converts the list (read from csv/mrm) file into a object model using the
     * LightTargetedExperiment with proper hierarchical structure from
     * Transition to Peptide to Protein.
     *
    */
    void TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenSwath::LightTargetedExperiment& exp);

    /// Convert an OpenMS transition to a TSVTransition for output writing
    TransitionTSVFile::TSVTransition convertTransition_(const ReactionMonitoringTransition* it, OpenMS::TargetedExperiment& targeted_exp);
    //@}

    /// Synchronize members with param class
    void updateMembers_() override;

private:

    // Members
    String retentionTimeInterpretation_;
    bool override_group_label_check_;
    bool force_invalid_mods_;

    // Typedefs
    typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
    typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
    typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

    static const char* strarray_[];

    static const std::vector<std::string> header_names_;

    /** @name Reader helper functions
     *
    */
    //@{

    /** @brief Determine separator in a CSV file and check for correct headers
     *
     * @param line The header to be parsed
     * @param delimiter The delimiter which will be determined from the input
     * @param header_dict The map which maps the fields in the header to their position
     *
    */
    void getTSVHeader_(const std::string& line, char& delimiter, std::map<std::string, int>& header_dict) const;

    /** @brief Read tab or comma separated input with columns defined by their column headers only
     *
     * @param filename The input file
     * @param filetype The type of file ("mrm" or "tsv")
     * @param transition_list The output list of transitions
     *
    */
    void readUnstructuredTSVInput_(const char* filename, FileTypes::Type filetype, std::vector<TSVTransition>& transition_list);

    /// Extract retention time from a SpectraST comment string
    void spectrastRTExtract(const String& str_inp, double & value, bool & spectrast_legacy);

    /// Extract annotation from a SpectraST comment string
    bool spectrastAnnotationExtract(const String& str_inp, TSVTransition & mytransition);

    /** @brief Cleanup of the read fields (removing quotes etc.)
    */
    void cleanupTransitions_(TSVTransition& mytransition);
    //@}

    /** @name Conversion helper functions
     *
    */
    //@{

    /** @brief Resolve cases where the same peptide label group has different sequences.
     *
     * Since members in a peptide label group (MS:1000893) should only be
     * isotopically modified forms of the same peptide, having different
     * peptide sequences (different AA sequences) within the same group most likely
     * constitutes an error. This function will fix the error by erasing the
     * provided "peptide group label" for a peptide and replace it with the
     * peptide identifier (transition group id).
     *
     * @param transition_list The list of transitions to be fixed.
     *
     */
    void resolveMixedSequenceGroups_(std::vector<TSVTransition>& transition_list) const;

    /// Populate a new ReactionMonitoringTransition object from a row in the csv
    void createTransition_(std::vector<TSVTransition>::iterator& tr_it,
                           OpenMS::ReactionMonitoringTransition& rm_trans);

    /// Populate a new TargetedExperiment::Protein object from a row in the csv
    void createProtein_(String protein_name, const String& uniprot_id,
                        OpenMS::TargetedExperiment::Protein& protein);

    /// Helper function to assign retention times to compounds and peptides
    void interpretRetentionTime_(std::vector<TargetedExperiment::RetentionTime>& retention_times,
                                 const OpenMS::DataValue& rt_value);

    /// Populate a new TargetedExperiment::Peptide object from a row in the csv
    void createPeptide_(std::vector<TSVTransition>::const_iterator tr_it,
                        OpenMS::TargetedExperiment::Peptide& peptide);

    /// Populate a new TargetedExperiment::Compound object (a metabolite) from a row in the csv
    void createCompound_(std::vector<TSVTransition>::const_iterator tr_it,
                         OpenMS::TargetedExperiment::Compound& compound);

    /// Add a modification at the specified location
    void addModification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
                          int location,
                          const ResidueModification& rmod);
    //@}

    /** @brief Write a TargetedExperiment to a file
     *
     * @param filename Name of the output file
     * @param targeted_exp The data structure to be written to the file
    */
    void writeTSVOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

public:

    //@{
    /// Constructor
    TransitionTSVFile();

    /// Destructor
    ~TransitionTSVFile() override;
    //@}

    /** @brief Write out a targeted experiment (TraML structure) into a tsv file
     *
     * @param filename The output file
     * @param targeted_exp The targeted experiment
     *
    */
    void convertTargetedExperimentToTSV(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

    /** @brief Read in a tsv/mrm file and construct a targeted experiment (TraML structure)
     *
     * @param filename The input file
     * @param filetype The type of file ("mrm" or "tsv")
     * @param targeted_exp The output targeted experiment
     *
    */
    void convertTSVToTargetedExperiment(const char* filename, FileTypes::Type filetype, OpenMS::TargetedExperiment& targeted_exp);

    /** @brief Read in a tsv file and construct a targeted experiment (Light transition structure)
     *
     * @param filename The input file
     * @param filetype The type of file ("mrm" or "tsv")
     * @param targeted_exp The output targeted experiment
     *
    */
    void convertTSVToTargetedExperiment(const char* filename, FileTypes::Type filetype, OpenSwath::LightTargetedExperiment& targeted_exp);

    /// Validate a TargetedExperiment (check that all ids are unique)
    void validateTargetedExperiment(const OpenMS::TargetedExperiment& targeted_exp);

  };
}

