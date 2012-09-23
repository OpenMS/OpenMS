// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMIDENTIFICATION_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMIDENTIFICATION_H

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

namespace OpenMS
{
  /**
     @deprecated Deprecated in OpenMS 1.7.

   @brief A map feature grouping algorithm for identified features.

   It takes many maps and searches for corresponding features.
   The corresponding features must be aligned, but may have small position deviations.

   @htmlinclude OpenMS_FeatureGroupingAlgorithmIdentification.parameters

   @ingroup FeatureGrouping
   */
  class OPENMS_DLLAPI FeatureGroupingAlgorithmIdentification :
    public FeatureGroupingAlgorithm
  {
public:
    /// Default constructor
    FeatureGroupingAlgorithmIdentification();

    /// Destructor
    virtual
    ~FeatureGroupingAlgorithmIdentification();

    /**
     @brief Applies the algorithm

     @exception IllegalArgument is thrown if less than two input maps are given.
     */
    virtual void
    group(const std::vector<FeatureMap<> > & maps, ConsensusMap & out);

    /// Creates a new instance of this class (for Factory)
    static FeatureGroupingAlgorithm *
    create()
    {
      return new FeatureGroupingAlgorithmIdentification();
    }

    /// Returns the product name (for the Factory)
    static String
    getProductName()
    {
      return "identification";
    }

private:

    /// Copy constructor intentionally not implemented -> private
    FeatureGroupingAlgorithmIdentification(const FeatureGroupingAlgorithmIdentification &);
    /// Assignment operator intentionally not implemented -> private
    FeatureGroupingAlgorithmIdentification &
    operator=(const FeatureGroupingAlgorithmIdentification &);

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMIDENTIFICATION_H
