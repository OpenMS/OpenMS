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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <map> //insert

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h> //SiriusMSFile

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  /**
  @brief This class provides methods for the extraction of targeted assays for metabolomics.

  */

  class OPENMS_DLLAPI MetaboTargetedAssay
  {
    public:

    /**
    @brief MetaboTargetedAssay is able to store a precursor and its metadata as well as a reference to a compound.

    */

    double precursor_int;
    double transition_quality_score;
    String compound_name;
    String compound_adduct;
    TargetedExperiment::Compound potential_cmp;
    std::vector<ReactionMonitoringTransition> potential_rmts;

    /**
    @brief CompoundSpectrumPair stores a pair of CompoundInfo and MSSpectrum

    */
    struct CompoundSpectrumPair
    {
      std::pair <SiriusMSFile::CompoundInfo, MSSpectrum> compoundspectrumpair;
    };


    /**
    @brief Extract a vector of MetaboTargetedAssays without using fragment annotation

    @return Vector of MetaboTargetedAssay



    @param spectra: Input of MSExperiment with spectra information
    @param feature_ms2_spectra_map: FeatureMapping class with associated MS2 spectra
    @param precursor_rt_tol: Retention time tolerance of the precursor
    @param precursor_mz_distance: Max m/z distance of the precursor entries of two spectra to be merged
    @param cosine_sim_threshold: Cosine similarty threshold for the usage of SpectraMerger
    @param transition_threshold: Intensity threshold for MS2 peak used in MetaboTargetedAssay
    @param method_consensus_spectrum: Boolean to use consensus spectrum method
    @param exclude_ms2_precursor: Boolean to exclude MS2 precursor from MetaboTargetedAssay
    @param file_counter: Count if multiple files are used.

    */

    static std::vector<MetaboTargetedAssay> extractMetaboTargetedAssay(const MSExperiment& spectra,
                                                                       const FeatureMapping::FeatureToMs2Indices& feature_ms2_index,
                                                                       const double& precursor_rt_tol,
                                                                       const double& precursor_mz_distance,
                                                                       const double& cosine_sim_threshold,
                                                                       const double& transition_threshold,
                                                                       const bool& method_consensus_spectrum,
                                                                       const bool& exclude_ms2_precursor,
                                                                       const unsigned int& file_counter);

    /**
    @brief Extract a vector of MetaboTargetedAssays using fragment annotation

    @return Vector of MetaboTargetedAssay


    @param v_cmp_spec: Vector of CompoundInfo with associated fragment annotated MSspectrum
    @param transition_threshold: Intensity threshold for MS2 peak used in MetaboTargetedAssay
    @param use_exact_mass: Boolean if exact mass should be used as peak mass for annotated fragments
    @param exclude_ms2_precursor: Boolean to exclude MS2 precursor from MetaboTargetedAssay
    @param file_counter: Count if multiple files are used.

    */

    static std::vector<MetaboTargetedAssay> extractMetaboTargetedAssayFragmentAnnotation(const std::vector< CompoundSpectrumPair >& v_cmp_spec,
                                                                                         const double& transition_threshold,
                                                                                         const bool& use_exact_mass,
                                                                                         const bool& exclude_ms2_precursor,
                                                                                         const unsigned int& file_counter);

    protected:

    /**
    @brief Compare two peaks based on their intensity
    */

    static bool intensityLess_(Peak1D a, Peak1D b);

};

} // namespace OpenMS
