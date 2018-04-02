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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_ANALYSIS_TARGETED_PRECURSORIONSELECTION_H
#define OPENMS_ANALYSIS_TARGETED_PRECURSORIONSELECTION_H

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/TARGETED/PSLPFormulation.h>
#include <set>

namespace OpenMS
{
  class PrecursorIonSelectionPreprocessing;
  class PSProteinInference;
  /**
    @brief This class implements different precursor ion selection strategies.

    @htmlinclude OpenMS_PrecursorIonSelection.parameters
  */
  class OPENMS_DLLAPI PrecursorIonSelection :
    public DefaultParamHandler
  {
public:

    /**
      @brief Precursor ion selection type (iterative, static, upshift, downshift, dynamic exclusion).

      The iterative strategy changes the ranking of possible precursors based on
      identification results from previous iterations.

      The upshift strategy assigns a higher priority to precursors  whose masses are matching
      peptide masses of potential protein identifications to enable a safe identification in the next iterations.

      The downshift strategy assigns a lower priority to precursors  whose masses are matching
      peptide masses of safe protein identifications.

      The dynamic exclusion excludes precursors whose masses are matching peptide masses of safe protein identifications.

      The static selection uses precomputed scores to rank the precursor, the order of precursors isn't
      changed throughout the run.
    */
    enum Type
    {
      IPS,
      ILP_IPS,
      SPS,
      UPSHIFT,
      DOWNSHIFT,
      DEX
    };

    PrecursorIonSelection();
    PrecursorIonSelection(const PrecursorIonSelection & source);
    ~PrecursorIonSelection() override;

    const double & getMaxScore() const;
    void setMaxScore(const double & max_score);


    /// Compare by score
    struct TotalScoreMore :
      std::binary_function<Feature, Feature, bool>
    {
      inline bool operator()(Feature const & left, Feature const & right) const
      {
        return (double)left.getMetaValue("msms_score") > (double)right.getMetaValue("msms_score");
      }

    };

    /// Compare by score
    struct SeqTotalScoreMore :
      std::binary_function<Feature, Feature, bool>
    {
      inline bool operator()(Feature const & left, Feature const & right) const
      {
        if (left.getRT() < right.getRT()) return true;
        else if (left.getRT() > right.getRT()) return false;
        else return (double)left.getMetaValue("msms_score") > (double)right.getMetaValue("msms_score");
      }

    };


    /**
      @brief Sort features by total score.
    */
    void sortByTotalScore(FeatureMap & features)
    {
      FeatureMap::Iterator beg = features.begin();
      FeatureMap::Iterator end  = features.end();
      std::sort(beg, end, TotalScoreMore());
    }

    /**
      @brief Returns features with highest score for MS/MS

      @param features FeatureMap with all possible precursors
      @param next_features FeatureMap with next precursors
      @param number Number of features to be reported
    */
    void getNextPrecursors(FeatureMap & features, FeatureMap & next_features, UInt number);
    void getNextPrecursorsSeq(FeatureMap & features, FeatureMap & next_features, UInt number, double & rt);
    void getNextPrecursors(std::vector<Int> & solution_indices, std::vector<PSLPFormulation::IndexTriple> & variable_indices, std::set<Int> & measured_variables,
                           FeatureMap & features, FeatureMap & new_features, UInt step_size, PSLPFormulation & ilp);

//      /**
//        @brief Change scoring of features using peptide identifications only from spectra of the last
//        iteration
//
//        @param features FeatureMap with all possible precursors
//        @param new_pep_ids Peptide identifications
//        @param preprocessed_db Information from preprocessed database
//
//     */
//     void rescoreIncremental(FeatureMap& features,std::vector<PeptideIdentification>& new_pep_ids,
//                                                      std::vector<ProteinIdentification>& prot_ids,
//                                                      PrecursorIonSelectionPreprocessing& preprocessed_db);


    /**
      @brief Change scoring of features using peptide identifications from all spectra.

      @param features FeatureMap with all possible precursors
      @param new_pep_ids Peptide identifications
      @param prot_ids Protein identifications
      @param preprocessed_db Information from preprocessed database
      @param check_meta_values True if the FeatureMap should be checked for the presence of required meta values
    */
    void rescore(FeatureMap & features, std::vector<PeptideIdentification> & new_pep_ids,
                 std::vector<ProteinIdentification> & prot_ids,
                 PrecursorIonSelectionPreprocessing & preprocessed_db, bool check_meta_values = true);


    /**
      @brief Simulate the iterative precursor ion selection.

      @param features FeatureMap with all possible precursors
      @param pep_ids Peptide identifications
      @param prot_ids Protein identifications
      @param preprocessed_db Information from preprocessed database
      @param step_size Number of MS/MS spectra considered per iteration
      @param path Path to output file
    */
    void simulateRun(FeatureMap & features, std::vector<PeptideIdentification> & pep_ids,
                     std::vector<ProteinIdentification> & prot_ids,
                     PrecursorIonSelectionPreprocessing & preprocessed_db,
                     String path, PeakMap & experiment, String precursor_path = "");

    void setLPSolver(LPWrapper::SOLVER solver)
    {
      solver_ = solver;
      std::cout << " LPSolver set to " << solver_ << std::endl;
    }

    LPWrapper::SOLVER getLPSolver()
    {
      return solver_;
    }

    void reset();

    const std::map<String, std::set<String> > & getPeptideProteinCounter()
    {
      return prot_id_counter_;
    }

private:
    void simulateILPBasedIPSRun_(FeatureMap & features, PeakMap & experiment,
                                 std::vector<PeptideIdentification> & pep_ids,
                                 std::vector<ProteinIdentification> & prot_ids,
                                 PrecursorIonSelectionPreprocessing & preprocessed_db,
                                 String output_path, String precursor_path = "");

    void simulateRun_(FeatureMap & features, std::vector<PeptideIdentification> & pep_ids,
                      std::vector<ProteinIdentification> & prot_ids,
                      PrecursorIonSelectionPreprocessing & preprocessed_db, String path, String precursor_path = "");

    void shiftDown_(FeatureMap & features, PrecursorIonSelectionPreprocessing & preprocessed_db, String protein_acc);

    void shiftUp_(FeatureMap & features, PrecursorIonSelectionPreprocessing & preprocessed_db, String protein_acc);

    /// update members method from DefaultParamHandler to update the members
    void updateMembers_() override;

    void rescore_(FeatureMap & features, std::vector<PeptideIdentification> & new_pep_ids,
                  PrecursorIonSelectionPreprocessing & preprocessed_db, PSProteinInference & protein_inference);

    /**
      @brief Adds user params, required for the use of IPS, to a feature map using default values.

      @param features FeatureMap with all possible precursors
    */
    void checkForRequiredUserParams_(FeatureMap & features);

    /**
      @brief Groups protein identifications that cannot be distinguished by their peptide identifications.

      @param prot_ids All protein identifications.
    */
    UInt filterProtIds_(std::vector<ProteinIdentification> & prot_ids);

    std::vector<PeptideIdentification> filterPeptideIds_(std::vector<PeptideIdentification> & pep_ids);

    void convertPeptideIdScores_(std::vector<PeptideIdentification> & pep_ids);

    /// minimal number of peptides identified for a protein to be declared identified
    UInt min_pep_ids_;
    /// maximal score in the FeatureMap
    double max_score_;
    /// precursor ion selection strategy
    Type type_;
    /// stores the peptide sequences for all protein identifications
    std::map<String, std::set<String> > prot_id_counter_;
    /// stores the number of selected precursors per fraction
    std::vector<Size> fraction_counter_;
    /// precursor ion error tolerance
    double mz_tolerance_;
    /// precursor ion error tolerance unit (ppm or Da)
    String mz_tolerance_unit_;
    /// maximal number of iterations
    UInt max_iteration_;
    Size x_variable_number_;

    LPWrapper::SOLVER solver_;

  };

}

#endif // #ifndef OPENMS_ANALYSIS_ID_PRECURSORIONSELECTION_H
