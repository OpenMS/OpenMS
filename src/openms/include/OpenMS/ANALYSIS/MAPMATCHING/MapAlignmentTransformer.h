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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTTRANSFORMER_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTTRANSFORMER_H

#include <vector>
#include <OpenMS/config.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  class TransformationDescription;
  class ConsensusMap;
  class PeptideIdentification;
  class ConsensusFeature;

  /**
   * @brief This class collects functions for applying retention time transformations to data structures.
   */
  class OPENMS_DLLAPI MapAlignmentTransformer
  {

  public:
    /// Applies the given transformation to a peak map
    static void transformRetentionTimes(MSExperiment<>& msexp,
                                        const TransformationDescription& trafo,
                                        bool store_original_rt = false);

    /// Applies the given transformation to a feature map
    static void transformRetentionTimes(
      FeatureMap& fmap, const TransformationDescription& trafo,
      bool store_original_rt = false);

    /// Applies the given transformation to a consensus map
    static void transformRetentionTimes(
      ConsensusMap& cmap, const TransformationDescription& trafo,
      bool store_original_rt = false);

    /// Applies the given transformation to peptide identifications
    static void transformRetentionTimes(
      std::vector<PeptideIdentification>& pep_ids,
      const TransformationDescription& trafo, bool store_original_rt = false);

  private:
    /// Applies a transformation to a feature
    static void applyToFeature_(Feature& feature,
                                const TransformationDescription& trafo,
                                bool store_original_rt = false);

    /// Applies a transformation to a basic feature
    static void applyToBaseFeature_(BaseFeature& feature,
                                    const TransformationDescription& trafo,
                                    bool store_original_rt = false);

    /// Applies a transformation to a consensus feature
    static void applyToConsensusFeature_(
      ConsensusFeature& feature, const TransformationDescription& trafo,
      bool store_original_rt = false);

    /**
       @brief Stores the original RT in a meta value

       The original RT is written to a meta value "original_RT", but only if that meta value doesn't already exist.

       @returns True if the meta value was written, false if it already exists.
    */
    static bool storeOriginalRT_(MetaInfoInterface& interface,
                                 double original_rt);

  };
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTTRANSFORMER_H
