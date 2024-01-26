// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantHelper.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
#include <algorithm>
#include <fstream>

using namespace std;
namespace OpenMS
{
  class OPENMS_DLLAPI FLASHQuantAlgorithm :
      public ProgressLogger,
      public DefaultParamHandler
{
  public:
    typedef FLASHQuantHelper::FeatureSeed FeatureSeed;
    typedef FLASHQuantHelper::FeatureGroup FeatureGroup;
    typedef FLASHQuantHelper::Feature Feature;

    /// Default constructor
    FLASHQuantAlgorithm();

    /// Default destructor
    ~FLASHQuantAlgorithm() = default;

    /// copy constructor
    FLASHQuantAlgorithm(const FLASHQuantAlgorithm& ) = delete;

    /// move constructor
    FLASHQuantAlgorithm(FLASHQuantAlgorithm&& other) = default;

    /// assignment operator
    FLASHQuantAlgorithm& operator=(const FLASHQuantAlgorithm& fd);

    /// main method of FeatureFindingMetabo
    void run(std::vector<MassTrace> &input_mtraces, std::vector<FeatureGroup> &out_fgroups);

    // test purpose
    String output_file_path_;

  protected:
    void updateMembers_() override;

  private:
    Param getFLASHDeconvParams_() const;

    // equivalent to FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_
    void getFeatureFromSpectrum_(std::vector<FeatureSeed*> &local_traces, std::vector<FeatureGroup> &local_fgroup, const double &rt);

    void buildMassTraceGroups_(std::vector<FeatureSeed> &in_seeds, std::vector<FeatureGroup> &features);

    bool scoreAndFilterFeatureGroup_(FeatureGroup& fg, double min_iso_score = -1) const;

    void refineFeatureGroups_(std::vector<FeatureGroup>& features);

    bool rescoreFeatureGroup_(FeatureGroup& fg) const;

    bool doFWHMbordersOverlap_(const std::pair<double, double>& border1, const std::pair<double, double>& border2) const;

    bool doMassTraceIndicesOverlap(const FeatureGroup& fg1, const FeatureGroup& fg2, const double overlap_percentage_threshold = 0.5, const bool charge_specific = true) const;

    void clusterFeatureGroups_(std::vector<FeatureGroup>& fgroups,
                               std::vector<MassTrace>& input_mtraces);

    void resolveConflictInCluster_(std::vector<FeatureGroup>& feature_groups,
                                   std::vector<MassTrace> &input_masstraces,
                                   std::vector<std::vector<Size> >& shared_m_traces_indices,
                                   std::set<Size>& hypo_indices,
                                   std::vector<FeatureGroup>& out_features);

    void filterOutIneligibleFeatureGroupsInCluster(std::vector<FeatureGroup>& feature_groups,
                                                   std::vector<std::vector<Size> >& shared_m_traces_indices,
                                                   std::set<Size>& hypo_indices) const;

    void setOptionalDetailedOutput_();

    void writeMassTracesOfFeatureGroup_(const FeatureGroup &fgroup,
                                       const Size &fgroup_idx,
                                       const std::vector<std::vector<Size>> &shared_m_traces_indices,
                                       const bool &is_before_resolution);

    void writeTheoreticalShapeForConflictResolution_(const Size &fgroup_idx,
                                                     const FeatureSeed &shared_mt,
                                                     const std::vector<double> &theo_intensities,
                                                     const double &calculated_ratio);

    void resolveConflictRegion_(std::vector<Feature> &conflicting_features,
                                const std::vector<MassTrace> &conflicting_mts,
                                const std::vector<Size> &conflicting_mt_indices);

    MassTrace updateMassTrace_(const MassTrace& ref_trace, const double &ratio) const;

    void fitTraceModelFromUniqueTraces_(Feature const& tmp_feat, EGHTraceFitter* fitted_model) const;

    void runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces, EGHTraceFitter* fitter) const;

    void updateFeatureWithFitModel(std::vector<Feature>& conflicting_features, Size mt_index,
                                   const MassTrace& obs_masstrace, const Size& org_index_of_obs_mt,
                                   Matrix<int>& pointer_to_components, vector<std::vector<double>>& components);

    void getMostAbundantMassTraceFromFeatureGroup_(const FeatureGroup &fgroup,
                                                  const int &ignore_this_charge,
                                                  FeatureSeed* &most_abundant_mt_ptr,
                                                  const std::vector<std::vector<Size>>& shared_m_traces) const;

    void getFLASHDeconvConsensusResult();

    bool isThisMassOneOfTargets_(const double &candi_mass, const double &candi_rt) const;

    void makeMSSpectrum_(std::vector<FeatureSeed *> &local_traces, MSSpectrum &spec, const double &rt) const;

    void setFeatureGroupMembersForResultWriting_(std::vector<FeatureGroup> &f_groups) const;

    bool isEligibleFeatureForConflictResolution_(Feature &new_feature, std::vector<std::vector<Size>> &shared_m_traces_indices, FeatureGroup &feat_group) const;

    /// default parameters
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    double min_mass_;
    double max_mass_;
    double mz_tolerance_; // ppm
    bool resolving_shared_signal_;

    double mass_tolerance_da_; // Da, for feature mass collection
    double iso_da_distance_ = Constants::ISOTOPE_MASSDIFF_55K_U;

    // advanced parameter?
    Size min_nr_mtraces_ = 3; //  minimum number of mass traces to support feature group
    Size min_nr_charges_ = 2; //  minimum number of consecutive bridges to support feature group
    Size min_nr_peaks_in_mtraces_ = 4; // at least 4 is needed for EGHTraceFitter
    bool use_smoothed_intensities_;

    /// variables for internal use (not for user input)
    FLASHDeconvHelperStructs::PrecalculatedAveragine iso_model_;
    Size max_nr_traces_; // calculated from iso_model_ (setAveragineModel())

    /// cosine threshold between observed and theoretical isotope patterns for MS1
    double min_isotope_cosine_;

    /// SpectralDeconvolution class for deconvolution
    SpectralDeconvolution deconv_;

    // loop up table
    std::vector<std::pair<double, double>> target_masses_; // mass and rt
    bool with_target_masses_ = false;

    /// for detailed outputs (mostly test purpose)
    bool shared_output_requested_;
    std::fstream shared_out_stream_ = std::fstream();
};
}