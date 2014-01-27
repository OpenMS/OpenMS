// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  /**@brief The main feature finder class.

      - Stores the flags for (indices of) data points ("used", "unused")
      - The algorithm itself is a factory product (derived from FeatureFinderAlgorithm)
      - The main method is run(), which is a template so that we can deal with different types of input and output
      - The run() method takes five arguments: algorithm_name, input_map, output, parameters, seeds
      .

      @ingroup FeatureFinder
  */
  class OPENMS_DLLAPI FeatureFinder :
    public ProgressLogger,
    public FeatureFinderDefs
  {

public:
    /// Default constructor.
    FeatureFinder();

    /// Destructor
    virtual ~FeatureFinder();

    /**
        @brief Executes the FeatureFinder using the given algorithm

        There are several constraints for the @p input_map.  They are tested before
        the algorithm starts.  It must only contain MS 1 level scans and you
        have to call updateRanges() before passing it to this method.
The input map is sorted by RT & m/z if that's not the case.
Furthermore we throw an Exception if the data contains negative m/z values,
as this will disturb most algorithms.

        @param algorithm_name Name of the feature finding algorithm to use
        @param input_map Input peak map
        @param features Output feature map
        @param param Algorithm parameters
        @param seeds List of seeds to use

        Implemented in FeatureFinder_impl.h
    */
    template <class PeakType, class FeatureType>
    void run(const String & algorithm_name, MSExperiment<PeakType> & input_map, FeatureMap<FeatureType> & features, const Param & param, const FeatureMap<FeatureType> & seeds);

    /// Returns a non-mutable reference to a peak flag
    const Flag & getPeakFlag(const IndexPair & index) const
    {
      return flags_[index.first][index.second];
    }

    /// Returns mutable reference to a peak flag
    Flag & getPeakFlag(const IndexPair & index)
    {
      return flags_[index.first][index.second];
    }

    /// Returns the default parameters for the algorithm with name @p algorithm_name
    Param getParameters(const String & algorithm_name) const;

protected:

    /// Container for flags attached to input data
    std::vector<std::vector<Flag> > flags_;

  };   // class FeatureFinder

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
