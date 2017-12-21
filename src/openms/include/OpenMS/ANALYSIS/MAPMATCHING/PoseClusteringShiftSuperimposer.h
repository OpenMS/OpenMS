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
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGSHIFTSUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGSHIFTSUPERIMPOSER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>

namespace OpenMS
{

  /**
   @brief  A superimposer that uses a voting scheme, also known as pose clustering,
   to find a good shift transformation.

   This algorithm works on two consensus maps.  It computes an shift
   transformation that maps the elements of second map as near as possible
   to the elements in the first map.

   The voting scheme hashes shift transformations between features in
   map one and features in map two.  Each such pair defines a
   (potential) "pose" of the second map relative to the first.
   Then it finds a cluster in the parameter space of these poses.
   The shift transformation is then computed from this
   cluster of potential poses, hence the name pose clustering.

   @sa PoseClusteringAffineSuperimposer

   @htmlinclude OpenMS_PoseClusteringShiftSuperimposer.parameters

   @ingroup MapAlignment
   */
  class OPENMS_DLLAPI PoseClusteringShiftSuperimposer :
    public BaseSuperimposer
  {
public:

    /// Default ctor
    PoseClusteringShiftSuperimposer();

    /// Destructor
    ~PoseClusteringShiftSuperimposer() override
    {}

    /**
      @brief Estimates the transformation and fills the given mapping function. (Has a precondition!)

      @note Exactly two input maps must be given.

      @pre For performance reasons, we trust that (the equivalent of:)
      <code>
        maps[0].updateRanges();
        maps[1].updateRanges();
      </code>
      has been done <i>before</i> calling this.  You have been warned!

      @exception IllegalArgument is thrown if the input maps are invalid.
    */
    void run(const ConsensusMap & map_model, const ConsensusMap & map_scene, TransformationDescription & transformation) override;

    /// Returns an instance of this class
    static BaseSuperimposer * create()
    {
      return new PoseClusteringShiftSuperimposer();
    }

    /// Returns the name of this module
    static const String getProductName()
    {
      return "poseclustering_shift";
    }

  };
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGSHIFTSUPERIMPOSER_H
