// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  class OPENMS_DLLAPI MapConversion
  {
public:

    /**
      @brief Similar to @p convert for FeatureMaps.

      Only the @p n most intense elements are copied.

      Currently PeakMap does not have a unique id but ConsensusMap has
      one, so we assign a new one here.

      @param input_map_index The index of the input map.
      @param input_map The input map to be converted.
      @param output_map The resulting ConsensusMap.
      @param n The maximum number of elements to be copied.
    */
    static void convert(UInt64 const input_map_index,
                        PeakMap& input_map,
                        ConsensusMap& output_map,
                        Size n = -1);

    /**
      @brief Convert a ConsensusMap to a FeatureMap (of any feature type).

      The previous content of output_map is cleared. UID's of the elements and
      the container is copied if the @p keep_uids flag is set.

      @param input_map The container to be converted.
      @param keep_uids Shall the UID's of the elements and the container be kept or created anew
      @param output_map The resulting ConsensusMap.
    */
    static void convert(ConsensusMap const& input_map,
                        const bool keep_uids,
                        FeatureMap& output_map);

    /**
      @brief Convert a FeatureMap (of any feature type) to a ConsensusMap.

      Each ConsensusFeature contains a map index, so this has to be given as
      well. The previous content of @p output_map is cleared. An arguable
      design decision is that the unique id of the FeatureMap is copied (!) to
      the ConsensusMap, because that is the way it is meant to be used in the
      algorithms.

      Only the first (!) @p n elements are copied. (This parameter exists
      mainly for compatibility with @p convert for MSExperiments. To use it in
      a meaningful way, apply one of the sorting methods to @p input_map
      beforehand.)

      @param input_map_index The index of the input map.
      @param input_map The container to be converted.
      @param output_map The resulting ConsensusMap.
      @param n The maximum number of elements to be copied.
    */
    static void convert(UInt64 const input_map_index,
                        FeatureMap const& input_map,
                        ConsensusMap& output_map,
                        Size n = -1);
  };
} // namespace OpenMS

