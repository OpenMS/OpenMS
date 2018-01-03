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

#ifndef OPENMS_ANALYSIS_TARGETED_PSLPFORMULATION_H
#define OPENMS_ANALYSIS_TARGETED_PSLPFORMULATION_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

//#define DEBUG_OPS
namespace OpenMS
{
  class PrecursorIonSelectionPreprocessing;
  class PSProteinInference;

  /**
    @brief Implements ILP formulation of precursor selection problems

    @htmlinclude OpenMS_PSLPFormulation.parameters
  */
  class OPENMS_DLLAPI PSLPFormulation :
    public DefaultParamHandler
  {



public:

    PSLPFormulation();

    ~PSLPFormulation() override;

    /**
      @brief Struct that holds the indices of the precursors in the feature map and the ilp formulation.
    */
    struct IndexTriple
    {
      Size feature;
      Int scan;
      Size variable;
      double rt_probability;
      double signal_weight;
      String prot_acc;
    };


    /**
      @brief Encode ILP formulation for a given LC-MS map, but unknown protein sample.

      @param features FeatureMap with all possible precursors
      @param experiment Input raw data
      @param variable_indices Assignment of feature indices and ILP variables
      @param mass_ranges Feature borders as indices in the raw data
      @param charges_set Allowed charge states
      @param ms2_spectra_per_rt_bin Allowed number of precursors per rt bin
      @param solution_indices Indices of ILP variables that are in the optimal solution
    */
    void createAndSolveILPForKnownLCMSMapFeatureBased(const FeatureMap& features,
                                                      const PeakMap& experiment,
                                                      std::vector<IndexTriple>& variable_indices,
                                                      std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                                                      std::set<Int>& charges_set, UInt ms2_spectra_per_rt_bin,
                                                      std::vector<int>& solution_indices);

    /**
      @brief Find a set of precursors, so that the protein coverage is maximal
      and that the number of precursors per bin is not exceeded
    */
    void createAndSolveILPForInclusionListCreation(PrecursorIonSelectionPreprocessing& preprocessing,
                                                   UInt ms2_spectra_per_rt_bin, UInt max_list_size,
                                                   FeatureMap& precursors,
                                                   bool solve_ILP = true);

    void createAndSolveCombinedLPForKnownLCMSMapFeatureBased(const FeatureMap& features,
                                                             const PeakMap& experiment,
                                                             std::vector<IndexTriple>& variable_indices,
                                                             std::vector<int>& solution_indices,
                                                             std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                                                             std::set<Int>& charges_set, UInt ms2_spectra_per_rt_bin,
                                                             Size step_size = 0, bool sequential_order = false);

    void updateStepSizeConstraint(Size iteration, UInt step_size);
    void updateFeatureILPVariables(FeatureMap& new_features, std::vector<IndexTriple>& variable_indices, std::map<Size, std::vector<String> >& feature_constraints_map);
    void updateRTConstraintsForSequentialILP(Size& rt_index, UInt ms2_spectra_per_rt_bin, Size max_rt_index);
    void updateCombinedILP(FeatureMap& features, PrecursorIonSelectionPreprocessing& preprocessed_db, std::vector<IndexTriple>& variable_indices,
                           std::vector<String>& new_protein_accs, std::vector<String>& protein_accs, PSProteinInference& prot_inference, Size& variable_counter,
                           std::map<String, std::vector<Size> >& protein_feature_map, Feature& new_feature, std::map<String, Size>& protein_variable_index_map,
                           std::map<String, std::set<String> >& prot_id_counter);


    /**
       @brief Solve the ILP.
    */
    void solveILP(std::vector<int>& solution_indices);

    void setLPSolver(LPWrapper::SOLVER solver)
    {
      solver_ = solver;
    }

    LPWrapper::SOLVER getLPSolver()
    {
      return solver_;
    }

    struct IndexLess :
      std::binary_function<IndexTriple, IndexTriple, bool>
    {
      inline bool operator()(IndexTriple  const& left,
                             IndexTriple const& right) const
      {
        return left.feature < right.feature;
      }

    };


    struct ScanLess :
      std::binary_function<IndexTriple, IndexTriple, bool>
    {
      inline bool operator()(IndexTriple  const& left,
                             IndexTriple  const& right) const
      {
        return left.scan < right.scan;
      }

    };

    struct VariableIndexLess :
      std::binary_function<IndexTriple, IndexTriple, bool>
    {
      inline bool operator()(IndexTriple  const& left,
                             IndexTriple  const& right) const
      {
        return left.variable < right.variable;
      }

    };

protected:

    void getXIC_(const std::vector<std::pair<Size, Size> >& end_points,
                 std::vector<double>& weights,
                 const PeakMap& experiment,
                 const bool normalize);

    /**
      @brief Calculates the XICs for all features.
    */
    void calculateXICs_(std::vector<std::vector<double> >& xics,
                        const FeatureMap& features,
                        const PeakMap& experiment,
                        const std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                        const bool normalize);

    /**
      @brief Creates and solves the ILP.
    */
    void createAndSolveILP_(const FeatureMap& features, std::vector<std::vector<double> >& intensity_weights,
                            std::set<Int>& charges_set, std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                            std::vector<IndexTriple>& variable_indices, std::vector<int>& solution_indices,
                            UInt ms2_spectra_per_rt_bin, Size number_of_scans);

    void createAndSolveCombinedLPFeatureBased_(const FeatureMap& features, std::vector<std::vector<double> >& intensity_weights,
                                               std::set<Int>& charges_set, std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                                               std::vector<IndexTriple>& variable_indices, std::vector<Int>& solution_indices,
                                               UInt ms2_spectra_per_rt_bin, Size number_of_scans, Size step_size = 0, bool sequential_order = false);

    void addProteinToILP_(PrecursorIonSelectionPreprocessing& preprocessing,
                          std::map<String, std::vector<double> >::const_iterator map_iter,
                          Size& counter, Size& pep_counter, Size& feature_counter,
                          std::vector<IndexTriple>& variable_indices,
                          std::map<String, Size>& protein_penalty_index_map, FeatureMap& precursors);

    void addPrecursorAcquisitionNumberConstraint_(std::vector<IndexTriple>& variable_indices, Size number_of_features, UInt number_of_msms_per_precursor);

    void addMaxInclusionListSizeConstraints_(std::vector<IndexTriple>& variable_indices, /*Size number_of_features,*/ UInt max_list_size);

    void addRTBinCapacityConstraint_(std::vector<IndexTriple>& variable_indices,
                                     Size max_rt_index, UInt ms2_spectra_per_rt_bin, bool sequential_order = false);

    void addProteinCoverageConstraint_(std::vector<IndexTriple>& variable_indices,
                                       PrecursorIonSelectionPreprocessing& preprocessing,
                                       std::map<String, Size> protein_variable_index_map);

    void addStepSizeConstraint_(std::vector<IndexTriple>& variable_indices, UInt step_size);


    void assembleInclusionListForProteinBasedLP_(std::vector<IndexTriple>& variable_indices, FeatureMap& precursors, std::vector<int>& solution_indices, PrecursorIonSelectionPreprocessing& preprocessing);

    void updateObjFunction_(String acc, FeatureMap& features, PrecursorIonSelectionPreprocessing& preprocessed_db, std::vector<IndexTriple>& variable_indices);


    Int getNumberOfPrecsInSpectrum_(Int constr_idx);

    LPWrapper* model_;
    LPWrapper::SOLVER solver_;
  };

  inline OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const PSLPFormulation::IndexTriple& triple)
  {
    os << "feature: " << triple.feature << " scan: " << triple.scan << " variable: " << triple.variable << " prot_acc: " << triple.prot_acc;
    return os;
  }

} // namespace

#endif // OPENMS_ANALYSIS_ID_PSLPFORMULATION_H
