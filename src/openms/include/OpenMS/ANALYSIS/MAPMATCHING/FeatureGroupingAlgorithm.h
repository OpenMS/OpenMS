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
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHM_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

  /**
      @brief Base class for all feature grouping algorithms

      These algorithms group corresponding features in one map or across maps.
  */
  class OPENMS_DLLAPI FeatureGroupingAlgorithm :
    public DefaultParamHandler
  {
public:
    /// Default constructor
    FeatureGroupingAlgorithm();

    /// Destructor
    ~FeatureGroupingAlgorithm() override;

    ///Applies the algorithm. The features in the input @p maps are grouped and the output is written to the consensus map @p out
    virtual void group(const std::vector<FeatureMap > & maps, ConsensusMap & out) = 0;

    ///Applies the algorithm. The consensus features in the input @p maps are grouped and the output is written to the consensus map @p out
    /// Algorithms not supporting ConsensusMap input should simply not override this method,
    /// as the base implementation will forward the data to the FeatureMap version of group()
    virtual void group(const std::vector<ConsensusMap> & maps, ConsensusMap & out);

    /// Transfers subelements (grouped features) from input consensus maps to the result consensus map
    void transferSubelements(const std::vector<ConsensusMap> & maps, ConsensusMap & out) const;

    /// Register all derived classes in this method
    static void registerChildren();

private:
    ///Copy constructor is not implemented -> private
    FeatureGroupingAlgorithm(const FeatureGroupingAlgorithm &);
    ///Assignment operator is not implemented -> private
    FeatureGroupingAlgorithm & operator=(const FeatureGroupingAlgorithm &);

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHM_H
