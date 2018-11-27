// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
      
      The following columns are required:

      <table>
        <tr> <td BGCOLOR="#EBEBEB">PrecursorMz*</td> <td>float</td> <td>This describes the precursor ion \a m/z</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ProductMz*</td> <td>float; synonyms: FragmentMz</td> <td>This specifies the product ion \a m/z</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">LibraryIntensity*</td> <td>float; synonyms: RelativeFragmentIntensity</td> <td>This specifies the relative intensity of the fragment ion</td></tr>
        <tr> <td BGCOLOR="#EBEBEB">NormalizedRetentionTime*</td> <td>float; synonyms: RetentionTime, Tr_recalibrated, iRT, RetentionTimeCalculatorScore</td> <td>This specifies the expected retention time (normalized retention time) </td></tr>
      </table>

  For targeted proteomics files, the following additional columns should be provided:
        <table>
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
          <tr> <td BGCOLOR="#EBEBEB">Decoy</td> <td>1: decoy, 0: target; synomys: decoy, isDecoy </td> <td> determines whether the transition is a decoy transition or not</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">PeptideGroupLabel</td> <td>free text </td> <td> designates to which peptide label group (as defined in MS:1000893) the peptide belongs to<sup>2</sup></td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">DetectingTransition</td> <td> 0 or 1; synonyms: detecting_transition</td> <td>1: use transition to detect peak group, 0: don't use transition for detection</td> </tr>
          <tr> <td BGCOLOR="#EBEBEB">IdentifyingTransition</td> <td> 0 or 1; synonyms: identifying_transition</td> <td>1: use transition for peptidoform inference using IPF, 0: don't use transition for identification</td> </tr>
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
          <tr> <td BGCOLOR="#EBEBEB">SMILES</td><td>free text</td><td>SMILES identifier</td></tr>
          <tr> <td BGCOLOR="#EBEBEB">SumFormula</td><td>free text</td><td>molecular formula </td></tr>
        </table>

  Fields indicated with * are strictly required to create a TraML file. Fields
  indicated with ** are recommended, but only required for a specific
  application (such as using the transition list for an analysis tool such as
  OpenSwath) or in a specific context (proteomics or metabolomics).

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
      double precursor;
      double product;
      double rt_calibrated;
      String transition_name;
      double CE;
      double library_intensity;
      String group_id;
      bool decoy;
      String PeptideSequence;
      String ProteinName;
      String Annotation;
      String FullPeptideName;
      String CompoundName;
      String SMILES;
      String SumFormula;
      String precursor_charge;
      String peptide_group_label;
      String label_type;
      String fragment_charge;
      int fragment_nr;
      double fragment_mzdelta;
      double drift_time;
      int fragment_modification;
      String fragment_type;
      String uniprot_id;
      bool detecting_transition;
      bool identifying_transition;
      bool quantifying_transition;
      std::vector<String> peptidoforms;

      TSVTransition() :
        precursor(-1),
        product(-1),
        rt_calibrated(-1),
        CE(-1),
        library_intensity(-1),
        decoy(false),
        fragment_charge("NA"),
        fragment_nr(-1),
        fragment_mzdelta(-1),
        drift_time(-1),
        fragment_modification(0),
        detecting_transition(true),
        identifying_transition(false),
        quantifying_transition(true)
      {}

      // By convention, if there is no (metabolic) compound name, it is a peptide 
      bool isPeptide() const
      {
        return CompoundName.empty();
      }
    };

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

    /** @name  Conversion functions from TSVTransition objects to TraML datastructures
     *
     * These functions convert the relevant data from a TSVTransition to the
     * datastructures used by the TraML handler, namely
     * ReactionMonitoringTransition, TargetedExperiment::Protein and
     * TargetedExperiment::Peptide.
     *
   */
    //@{

    TransitionTSVFile::TSVTransition convertTransition_(const ReactionMonitoringTransition* it, OpenMS::TargetedExperiment& targeted_exp);

    /// Synchronize members with param class
    void updateMembers_() override;

private:
    /// Members
    String retentionTimeInterpretation_;
    bool override_group_label_check_;
    bool force_invalid_mods_;

    /// Typedefs
    typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
    typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
    typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

    static const char* strarray_[];

    static const std::vector<std::string> header_names_;

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

    void spectrastRTExtract(const String str_inp, double & value, bool & spectrast_legacy);

    bool spectrastAnnotationExtract(const String str_inp, TSVTransition & mytransition);

    /** @brief Cleanup of the read fields (removing quotes etc.)
    */
    void cleanupTransitions_(TSVTransition& mytransition);

    /** Resolve cases where the same peptide label group has different sequences.
     *
     * Since members in a peptide label group (MS:1000893) should only be
     * isotopically modified forms of the same peptide, having different
     * peptide sequences (different AA order) within the same group most likely
     * constitutes an error. This function will fix the error by erasing the
     * provided "peptide group label" for a peptide and replace it with the
     * peptide id (transition group id).
     *
     * @param transition_list The list of read transitions to be fixed.
     *
     */
    void resolveMixedSequenceGroups_(std::vector<TSVTransition>& transition_list);

    /// Populate a new ReactionMonitoringTransition object from a row in the csv
    void createTransition_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::ReactionMonitoringTransition& rm_trans);

    /// Populate a new TargetedExperiment::Protein object from a row in the csv
    void createProtein_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Protein& protein);

    /// Helper function to assign retention times to compounds and peptides
    void interpretRetentionTime_(std::vector<TargetedExperiment::RetentionTime>& retention_times, const OpenMS::DataValue rt_value);

    /// Populate a new TargetedExperiment::Peptide object from a row in the csv
    void createPeptide_(std::vector<TSVTransition>::const_iterator tr_it, OpenMS::TargetedExperiment::Peptide& peptide);

    /// Populate a new TargetedExperiment::Compound object (a metabolite) from a row in the csv
    void createCompound_(std::vector<TSVTransition>::const_iterator tr_it, OpenMS::TargetedExperiment::Compound& compound);

    void addModification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
                          int location, const ResidueModification& rmod);
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

