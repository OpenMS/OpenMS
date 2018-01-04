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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_LABELEDPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_LABELEDPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

#include <boost/math/tr1.hpp>
#include <cmath>

namespace OpenMS
{
  /**
      @brief The LabeledPairFinder allows the matching of labeled features (features with a fixed distance).

      Finds feature pairs that have a defined distance in RT and m/z in the same map.

      @htmlinclude OpenMS_LabeledPairFinder.parameters

      @todo Implement support for labeled MRM experiments, Q1 m/z value and charges. (Andreas)
      @todo Implement support for more than one mass delta, e.g. from missed cleavages and so on (Andreas)

      @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI LabeledPairFinder :
    public BaseGroupFinder
  {

public:

    /// Default constructor
    LabeledPairFinder();

    /// Destructor
    inline ~LabeledPairFinder() override
    {
    }

    /// Returns an instance of this class
    static BaseGroupFinder * create()
    {
      return new LabeledPairFinder();
    }

    /// Returns the name of this module
    static const String getProductName()
    {
      return "labeled_pair_finder";
    }

    /**
      @brief Run the algorithm

      @note Exactly one @em input map has to be provided.
      @note The @em output map has to have two file descriptions, containing
      the same file name. The file descriptions have to be labeled 'heavy'
      and 'light'.

      @exception Exception::IllegalArgument is thrown if the input data is not valid.
    */
    void run(const std::vector<ConsensusMap> & input_maps, ConsensusMap & result_map) override;

protected:

    /// return the p-value at position x for the bi-Gaussian distribution with mean @p m and standard deviation @p sig1 (left) and @p sig2 (right)
    inline double PValue_(double x, double m, double sig1, double sig2)
    {
      if (m < x)
      {
        return 1 - boost::math::tr1::erf((x - m) / sig2 / 0.707106781);
      }
      else
      {
        return 1 - boost::math::tr1::erf((m - x) / sig1 / 0.707106781);
      }
    }

private:

    /// Copy constructor not implemented => private
    LabeledPairFinder(const LabeledPairFinder & source);

    /// Assignment operator not implemented => private
    LabeledPairFinder & operator=(const LabeledPairFinder & source);

  };   // end of class LabeledPairFinder

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_LABELEDPAIRFINDER_H
