// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_TARGETED_PSLPFORMULATION_H
#define OPENMS_ANALYSIS_TARGETED_PSLPFORMULATION_H
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/LPWrapper.h>

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

    virtual ~PSLPFormulation();

    /**
      @brief Struct that holds the indices of the precursors in the feature map and the ilp formulation.
    */
    struct IndexTriple
    {
      Size feature;
      Int scan;
      Size variable;
      DoubleReal rt_probability;
      DoubleReal signal_weight;
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
    template <typename InputPeakType>
    void createAndSolveILPForKnownLCMSMapFeatureBased(const FeatureMap<> & features,
                                                      const MSExperiment<InputPeakType> & experiment,
                                                      std::vector<IndexTriple> & variable_indices,
                                                      std::vector<std::vector<std::pair<Size, Size> > > & mass_ranges,
                                                      std::set<Int> & charges_set, UInt ms2_spectra_per_rt_bin,
                                                      std::vector<int> & solution_indices);

    /**
      @brief Find a set of precursors, so that the protein coverage is maximal
      and that the number of precursors per bin is not exceeded
    */
    void createAndSolveILPForInclusionListCreation(PrecursorIonSelectionPreprocessing & preprocessing,
                                                   UInt ms2_spectra_per_rt_bin, UInt max_list_size,
                                                   FeatureMap<> & precursors,
                                                   bool solve_ILP = true);

    template <typename InputPeakType>
    void createAndSolveCombinedLPForKnownLCMSMapFeatureBased(const FeatureMap<> & features,
                                                             const MSExperiment<InputPeakType> & experiment,
                                                             std::vector<IndexTriple> & variable_indices,
                                                             std::vector<int> & solution_indices,
                                                             std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
                                                             std::set<Int> & charges_set, UInt ms2_spectra_per_rt_bin,
                                                             Size step_size = 0, bool sequential_order = false);

    void updateStepSizeConstraint(Size iteration, UInt step_size);
    void updateFeatureILPVariables(FeatureMap<> & new_features, std::vector<IndexTriple> & variable_indices, std::map<Size,std::vector<String> > & feature_constraints_map);
    void updateRTConstraintsForSequentialILP(Size & rt_index, UInt ms2_spectra_per_rt_bin, Size max_rt_index);
    void updateCombinedILP(FeatureMap<> & features, PrecursorIonSelectionPreprocessing & preprocessed_db, std::vector<IndexTriple> & variable_indices,
                           std::vector<String> & new_protein_accs, std::vector<String> & protein_accs, PSProteinInference & prot_inference, Size & variable_counter,
                           std::map<String,std::vector<Size> > & protein_feature_map, Feature& new_feature, std::map<String,Size> & protein_variable_index_map,
                           std::map<String,std::set<String> > & prot_id_counter);

    
    /**
       @brief Solve the ILP.
    */
    void solveILP(std::vector<int> & solution_indices);
    
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
      inline bool operator()(IndexTriple  const & left,
                             IndexTriple const & right) const
      {
        return left.feature < right.feature;
      }

    };


    struct ScanLess :
      std::binary_function<IndexTriple, IndexTriple, bool>
    {
      inline bool operator()(IndexTriple  const & left,
                             IndexTriple  const & right) const
      {
        return left.scan < right.scan;
      }

    };

    struct VariableIndexLess :
      std::binary_function<IndexTriple, IndexTriple, bool>
    {
      inline bool operator()(IndexTriple  const & left,
                             IndexTriple  const & right) const
      {
        return left.variable < right.variable;
      }

    };

protected:

    template <typename InputPeakType>
    void getXIC_(const std::vector<std::pair<Size, Size> > & end_points,
                 std::vector<DoubleReal> & weights,
                 const MSExperiment<InputPeakType> & experiment,
                 const bool normalize);

    /**
      @brief Calculates the XICs for all features.
    */
    template <typename InputPeakType>
    void calculateXICs_(std::vector<std::vector<DoubleReal> > & xics,
                        const FeatureMap<> & features,
                        const MSExperiment<InputPeakType> & experiment,
                        const std::vector<std::vector<std::pair<Size, Size> > > & mass_ranges,
                        const bool normalize);

    /**
      @brief Creates and solves the ILP.
    */
    void createAndSolveILP_(const FeatureMap<> & features, std::vector<std::vector<DoubleReal> > & intensity_weights,
                            std::set<Int> & charges_set, std::vector<std::vector<std::pair<Size, Size> > > & mass_ranges,
                            std::vector<IndexTriple> & variable_indices, std::vector<int> & solution_indices,
                            UInt ms2_spectra_per_rt_bin, Size number_of_scans);

    void createAndSolveCombinedLPFeatureBased_(const FeatureMap<> & features, std::vector<std::vector<DoubleReal> > & intensity_weights,
                                               std::set<Int> & charges_set, std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
                                               std::vector<IndexTriple> & variable_indices, std::vector<Int> & solution_indices,
                                               UInt ms2_spectra_per_rt_bin, Size number_of_scans, Size step_size = 0, bool sequential_order = false);

    void addProteinToILP_(PrecursorIonSelectionPreprocessing & preprocessing,
                          std::map<String, std::vector<DoubleReal> >::const_iterator map_iter,
                          Size & counter, Size & pep_counter, Size & feature_counter,
                          std::vector<IndexTriple> & variable_indices,
                          std::map<String, Size> & protein_penalty_index_map, FeatureMap<> & precursors);

    void addPrecursorAcquisitionNumberConstraint_(std::vector<IndexTriple> & variable_indices, Size number_of_features, UInt number_of_msms_per_precursor);
    
    void addMaxInclusionListSizeConstraints_(std::vector<IndexTriple> & variable_indices, /*Size number_of_features,*/ UInt max_list_size);

    void addRTBinCapacityConstraint_(std::vector<IndexTriple> & variable_indices,
                                     Size max_rt_index, UInt ms2_spectra_per_rt_bin, bool sequential_order = false);

    void addProteinCoverageConstraint_(std::vector<IndexTriple> & variable_indices,
                                       PrecursorIonSelectionPreprocessing & preprocessing,
                                       std::map<String, Size> protein_variable_index_map);

    void addStepSizeConstraint_(std::vector<IndexTriple> & variable_indices, UInt step_size);


    void assembleInclusionListForProteinBasedLP_(std::vector<IndexTriple> & variable_indices, FeatureMap<> & precursors, std::vector<int> & solution_indices, PrecursorIonSelectionPreprocessing & preprocessing);

    void updateObjFunction_(String acc, FeatureMap<> & features, PrecursorIonSelectionPreprocessing & preprocessed_db, std::vector<IndexTriple>& variable_indices);


   Int getNumberOfPrecsInSpectrum_(Int constr_idx);
   
    LPWrapper * model_;
    LPWrapper::SOLVER solver_;
  };

  template <typename InputPeakType>
  void PSLPFormulation::getXIC_(const std::vector<std::pair<Size, Size> > & end_points,
                                std::vector<DoubleReal> & weights,
                                const MSExperiment<InputPeakType> & experiment,
                                const bool normalize)
  {
    DoubleReal max_weight = 0.;
    weights.clear();
    for (Size i = 0; i < end_points.size(); i += 2)
    {
      DoubleReal weight = 0.;
      for (Size j = end_points[i].second; j <= end_points[i + 1].second; ++j)
      {
        weight += experiment[end_points[i].first][j].getIntensity();
        // std::cout << " add "<<experiment[end_points[i].first][j].getIntensity()<<std::endl;
      }
      if (weight > max_weight)
        max_weight = weight;

      weights.push_back(weight);
    }

    if (normalize)
    {
      // normalize weights
      for (Size i = 0; i < weights.size(); ++i)
      {
#ifdef DEBUG_OPS
        if (end_points.size() >= i)
        {
          std::cout << "scan " << end_points[i].first << " " << weights[i] << " " << max_weight
                    << " " << weights[i] / max_weight << std::endl;
        }
#endif
        weights[i] /= max_weight;
      }
    }
  }

  template <typename InputPeakType>
  void PSLPFormulation::calculateXICs_(std::vector<std::vector<DoubleReal> > & xics,
                                       const FeatureMap<> & features,
                                       const MSExperiment<InputPeakType> & experiment,
                                       const std::vector<std::vector<std::pair<Size, Size> > > & mass_ranges,
                                       const bool normalize)
  {
    xics.clear();
    xics.resize(features.size());
    for (Size i = 0; i < features.size(); ++i)
    {
      getXIC_(mass_ranges[i], xics[i], experiment, normalize);
    }
  }

  template <typename InputPeakType>
  void PSLPFormulation::createAndSolveILPForKnownLCMSMapFeatureBased(const FeatureMap<> & features,
                                                                     const MSExperiment<InputPeakType> & experiment,
                                                                     std::vector<IndexTriple> & variable_indices,
                                                                     std::vector<std::vector<std::pair<Size, Size> > > & mass_ranges,
                                                                     std::set<Int> & charges_set, UInt ms2_spectra_per_rt_bin,
                                                                     std::vector<int> & solution_indices)
  {

    std::vector<std::vector<DoubleReal> > intensity_weights;
    calculateXICs_(intensity_weights, features, experiment, mass_ranges, true);
#ifdef DEBUG_OPS
    std::cout << "got xics" << std::endl;
#endif

    createAndSolveILP_(features, intensity_weights, charges_set, mass_ranges, variable_indices, solution_indices,
                       ms2_spectra_per_rt_bin, experiment.size());
  }

  inline OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const PSLPFormulation::IndexTriple & triple)
  {
    os << "feature: " << triple.feature << " scan: " << triple.scan << " variable: " << triple.variable << " prot_acc: " << triple.prot_acc;
    return os;
  }

  template <typename InputPeakType>
	void PSLPFormulation::createAndSolveCombinedLPForKnownLCMSMapFeatureBased(const FeatureMap<> & features,
                                                                            const MSExperiment<InputPeakType> & experiment,
                                                                            std::vector<IndexTriple> & variable_indices,
                                                                            std::vector<Int> & solution_indices,
                                                                            std::vector<std::vector<std::pair<Size,Size> > > & mass_ranges,
                                                                            std::set<Int> & charges_set, UInt ms2_spectra_per_rt_bin,
                                                                            Size step_size, bool sequential_order)
	{
    
		std::vector<std::vector<DoubleReal> > intensity_weights;
		calculateXICs_(intensity_weights, features, experiment, mass_ranges, true);
#ifdef DEBUG_OPS
		std::cout << "got xics"<<std::endl;
#endif
		
		createAndSolveCombinedLPFeatureBased_(features, intensity_weights, charges_set, mass_ranges, variable_indices, solution_indices, ms2_spectra_per_rt_bin,
                                          experiment.size(), step_size, sequential_order);
	}
	

 

  
} // namespace

#endif // OPENMS_ANALYSIS_ID_PSLPFORMULATION_H
