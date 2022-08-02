// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <algorithm>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include "boost/dynamic_bitset.hpp"
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <fstream>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvQuantHelper.h>

using namespace std;
namespace OpenMS
{
  class OPENMS_DLLAPI FLASHDeconvQuantAlgorithm :
      public ProgressLogger,
      public DefaultParamHandler
{
  public:
    typedef FLASHDeconvQuantHelper::FeatureSeed FeatureSeed;
    typedef FLASHDeconvQuantHelper::FeatureGroup FeatureGroup;
    typedef FLASHDeconvQuantHelper::FeatureElement FeatureElement;

    /// Default constructor
    FLASHDeconvQuantAlgorithm();

    /// Default destructor
    ~FLASHDeconvQuantAlgorithm() override;

    /// main method of FeatureFindingMetabo
    void run(std::vector<MassTrace> &input_mtraces, FeatureMap &output_featmap);

    String inputfile_path;
    String outfile_path;
    bool outFeatureXML = false;

  protected:
    void updateMembers_() override;

  private:
    Param getFLASHDeconvParams_();

    // equivalent to FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_
    void getFeatureFromSpectrum_(std::vector<FeatureSeed*> &local_traces, std::vector<FeatureGroup> &local_fgroup, const double &rt);

    void buildMassTraceGroups_(std::vector<FeatureSeed> &in_seeds, std::vector<FeatureGroup> &features);

    bool scoreAndFilterFeatureGroup_(FeatureGroup& fg) const;

    void calculatePerChargeIsotopeIntensity_(std::vector<double> &per_isotope_intensity,
                                                        std::vector<double> &per_charge_intensity,
                                                        FeatureGroup &fg) const;

    void refineFeatureGroups_(std::vector<FeatureGroup>& features);

    bool rescoreFeatureGroup_(FeatureGroup& fg, bool score_always = false) const;

    void setFeatureGroupScore_(FeatureGroup &fg) const;

    double scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const;

    double scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const;

    double computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const;

    bool doFWHMbordersOverlap(const std::pair<double, double>& border1, const std::pair<double, double>& border2) const;

    bool doMassTraceIndicesOverlap(const FeatureGroup& fg1, const FeatureGroup& fg2) const;

    void clusterFeatureGroups_(std::vector<FeatureGroup>& fgroups,
                               std::vector<MassTrace>& input_mtraces) const;

    void resolveConflictInCluster_(std::vector<FeatureGroup>& feature_groups,
                                   const std::vector<MassTrace> & input_masstraces,
                                   const std::vector<std::vector<Size> >& shared_m_traces_indices,
                                   const std::set<Size>& hypo_indices,
                                   std::vector<FeatureGroup>& out_features) const;

    void writeMassTracesOfFeatureGroup(const std::vector<FeatureGroup>& featgroups,
                                       const std::vector<std::vector<Size> >& shared_m_traces_indices) const;

    void writeFeatureGroupsInTsvFile(std::vector<FeatureGroup>& featgroups) const;

    void writeOutputInFeatureXML_(const std::vector<FeatureGroup> &feature_groups,
                                  const std::vector<std::vector<Size>> &shared_m_traces_indices) const;

    void storeFeatureGroupInOpenMSFeature(std::vector<FeatureGroup> &feature_groups,
                                          FeatureMap &out_featmap) const;

    void resolveConflictRegion_(std::vector<FeatureElement> &feat,
                                const std::vector<Size> &feature_idx_in_current_conflict_region,
                                const std::vector<const MassTrace*> &conflicting_mts) const;

    void runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces, EGHTraceFitter* fitter) const;

    void getMostAbundantMassTraceFromFeatureGroup(const FeatureGroup &fgroup,
                                                  const int &ignore_this_charge,
                                                  FeatureSeed* &most_abundant_mt_ptr,
                                                  const std::vector<std::vector<Size>>& shared_m_traces) const;

    void getFLASHDeconvConsensusResult();

    bool isThisMassOneOfTargets(const double &candi_mass, const double &candi_rt) const;

    void makeMSSpectrum_(std::vector<FeatureSeed *> &local_traces, MSSpectrum &spec, const double &rt) const;

    /// parameter stuff
//    double local_rt_range_;
    double local_mz_range_;
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    int charge_range_;
    double min_mass_;
    double max_mass_;
    double mz_tolerance_; // ppm

    double total_intensity_;

    const double mass_tolerance_da_ = 3; // Da, for feature mass collection
//    const double mass_tolerance_ppm_ = 20;

    // advanced parameter?
    Size min_nr_mtraces_ = 3; // minimum number of consecutive bridges among mass traces to support feature
    bool use_smoothed_intensities_;
    double rt_window_ = 1; // TODO : remove?

    /// variables for internal use (not for user input)
    FLASHDeconvHelperStructs::PrecalculatedAveragine iso_model_;
    Size max_nr_traces_; // calculated from iso_model_ (setAveragineModel())

    /// cosine threshold between observed and theoretical isotope patterns for MS1
    double min_isotope_cosine_ = 0.90;

    /// FLASHDeconvAlgorithm class for deconvolution
    FLASHDeconvAlgorithm fd_;

    // loop up table
    std::vector<std::pair<double, double>> target_masses_; // mass and rt
    bool with_target_masses_ = false;

  };
}