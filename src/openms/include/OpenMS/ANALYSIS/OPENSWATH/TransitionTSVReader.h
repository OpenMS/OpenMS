// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_TRANSITIONTSVREADER_H
#define OPENMS_ANALYSIS_OPENSWATH_TRANSITIONTSVREADER_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <fstream>

namespace OpenMS
{

  /**
      @brief This class can convert TraML and TSV files into each other

      The TSV are tab-separated and need to have the following columns:

      PrecursorMz (float)
      ProductMz (float)
      Tr_calibrated (float)
      transition_name (free text, needs to be unique for each transition [in this file])
      CE (float)
      LibraryIntensity (float)
      transition_group_id (free text, designates the transition group [e.g. peptide] to which this transition belongs)
      decoy (1==decoy, 0== no decoy; determines whether the transition is a decoy transition or not)
      PeptideSequence  (free text, sequence only (no modifications) )
      ProteinName  (free text)
      Annotation  (free text, e.g. y7)
      FullUniModPeptideName  (free text, should contain modifications*)
      MissedCleavages
      Replicates
      NrModifications
      PrecursorCharge (integer)
      PeptideGroupLabel (free text, designates to which peptide label group (as defined in MS:1000893) the peptide belongs to)
      LabelType (free text, optional description of which label was used, e.g. heavy or light)
      detecting_transition (bool, should this transition be used for peak-picking (detection) of the peak group?)
      identifying_transition (bool, should this transition be used for UIS identification together with the detecting transitions?)
      quantifying_transition (bool, should this transition be used for quantification?)

  @htmlinclude OpenMS_TransitionTSVReader.parameters

  */
  class OPENMS_DLLAPI TransitionTSVReader :
    public ProgressLogger,
    public DefaultParamHandler
  {

private:
    /// Members
    String retentionTimeInterpretation_;
    bool override_group_label_check_;

    /// Typedefs
    typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
    typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
    typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

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
      int decoy;
      String PeptideSequence;
      String ProteinName;
      String Annotation;
      String FullPeptideName;
      int precursor_charge;
      String peptide_group_label;
      String label_type;
      int product_charge;
      int fragment_nr;
      double fragment_mzdelta;
      int fragment_modification;
      String fragment_type;
      String uniprot_id;
      bool detecting_transition;
      bool identifying_transition;
      bool quantifying_transition;
    };

    static const char* strarray_[];

    static const std::vector<std::string> header_names_;

    /** @brief Determine separator in a CSV file and check for correct headers
     *
     * @param line The header to be parsed
     * @param delimiter The delimiter which will be determined from the input
     * @param header The fields of the header
     * @param header_dict The map which maps the fields in the header to their position
     *
    */
    void getTSVHeader_(const std::string& line, char& delimiter, std::vector<std::string> header, std::map<std::string, int>& header_dict);

    /** @brief Read tab or comma separated input with columns defined by their column headers only
     *
     * @param filename The input file
     * @param filetype The type of file ("mrm" or "tsv")
     * @param transition_list The output list of transitions
     *
    */
    void readUnstructuredTSVInput_(const char* filename, FileTypes::Type filetype, std::vector<TSVTransition>& transition_list);

    /** @brief Cleanup of the read fields (removing quotes etc.)
    */
    void cleanupTransitions_(TSVTransition& mytransition);

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

    void createTransition_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::ReactionMonitoringTransition& rm_trans);

    void createProtein_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Protein& protein);

    void createPeptide_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Peptide& peptide);

    void addModification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
                          int location, ResidueModification& rmod, const String& name);
    //@}


    /** @brief Write a TargetedExperiment to a file
     *
     * @param filename Name of the output file
     * @param targeted_exp The data structure to be written to the file
    */
    void writeTSVOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

protected:

    /// Synchronize members with param class
    void updateMembers_();

public:

    //@{
    /// Constructor
    TransitionTSVReader();

    /// Destructor
    ~TransitionTSVReader();
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
    void validateTargetedExperiment(OpenMS::TargetedExperiment& targeted_exp);

  };
}

#endif
