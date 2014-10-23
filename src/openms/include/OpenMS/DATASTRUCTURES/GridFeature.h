// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_GRIDFEATURE_H
#define OPENMS_DATASTRUCTURES_GRIDFEATURE_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <set>

namespace OpenMS
{
  class BaseFeature;
  class AASequence;

  /**
   * @brief Representation of a feature in a hash grid.
   *
   * A GridFeature can be stored in a HashGrid and points to a BaseFeature (Feature or ConsensusFeature). Used for QT feature grouping (see QTClusterFinder).
   */
  class OPENMS_DLLAPI GridFeature
  {
private:
    /// Reference to the contained feature
    const BaseFeature & feature_;

    /// Index of the feature map or consensus map
    Size map_index_;

    /// Index of the feature in the map
    Size feature_index_;

    /// Set of peptide sequences annotated to the feature
    std::set<AASequence> annotations_;

public:
    /**
     * @brief Detailed constructor
     * @param feature Reference to the contained feature
     * @param map_index Index of the feature map or consensus map
     * @param feature_index Index of the feature in the map
     */
    GridFeature(const BaseFeature & feature, Size map_index, Size feature_index);

    /// Returns the feature
    const BaseFeature & getFeature() const;

    /// Destructor
    virtual ~GridFeature();

    /// Returns the map index
    Size getMapIndex() const;

    /// Returns the feature index
    Size getFeatureIndex() const;

    /// Returns the ID of the GridFeature (same as the feature index)
    Int getID() const;

    /// Returns the set of peptide sequences annotated to the cluster center
    const std::set<AASequence> & getAnnotations() const;

    /// Returns the feature RT
    double getRT() const;

    /// Returns the feature m/z
    double getMZ() const;
  };
}

#endif // OPENMS_DATASTRUCTURES_GRIDFEATURE_H
