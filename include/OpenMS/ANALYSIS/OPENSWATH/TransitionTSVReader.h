// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
      Labelgroup (free text, e.g. heavy or light)

  */
  class OPENMS_DLLAPI TransitionTSVReader :
    public ProgressLogger
  {

private:

    typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
    typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
    typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

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
      String group_label;
      int fragment_charge;
      int fragment_nr;
      String fragment_type;
      String uniprot_id;
    };

    static const char* strarray_[];

    static const std::vector<std::string> header_names_;

    /// read TSV input with columns in defined order
    void readTSVInput_(const char* filename, std::vector<TSVTransition>& transition_list);

    /// determine separater in a CSV file and check for correct headers
    void getTSVHeader_(std::string & line, char & delimiter, std::vector<std::string> header, std::map<std::string, int> & header_dict);

    /// read tab or comma separated input with columns defined by their column headers only
    void readUnstructuredTSVInput_(const char* filename, std::vector<TSVTransition>& transition_list);

    /// do post-processing on read input data (removing quotes etc)
    void cleanupTransitions_(TSVTransition & mytransition);

    /// store a list of TSVTransition objects properly in a TargetedExperiment
    void TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenMS::TargetedExperiment& exp);

    void TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenSwath::LightTargetedExperiment& exp);

    void createTransition_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::ReactionMonitoringTransition& rm_trans);

    void createProtein_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Protein& protein);

    void createPeptide_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Peptide& peptide);

    void addModification_(std::vector<TargetedExperiment::Peptide::Modification> & mods,
          int location, ResidueModification & rmod, const String & name);


    /// write a TargetedExperiment to a file
    void writeTSVOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

public:
    /// Write out a targeted experiment (TraML structure) into a tsv file
    void convertTargetedExperimentToTSV(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

    /// Read in a tsv file and construct a targeted experiment (TraML structure)
    void convertTSVToTargetedExperiment(const char* filename, OpenMS::TargetedExperiment& targeted_exp);

    /// Read in a tsv file and construct a targeted experiment (Light transition structure)
    void convertTSVToTargetedExperiment(const char* filename, OpenSwath::LightTargetedExperiment& targeted_exp);

    /// Validate a TargetedExperiment (check that all ids are unique)
    void validateTargetedExperiment(OpenMS::TargetedExperiment& targeted_exp);

  };
}

#endif

