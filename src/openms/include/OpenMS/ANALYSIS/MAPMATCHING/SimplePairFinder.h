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
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>

#define V_SimplePairFinder(bla) // std::cout << bla << std::endl;

namespace OpenMS
{

  /** @brief This class implements a simple point pair finding algorithm.

    It offers a method to find element pairs across two element maps.

    The similarity value should express our confidence that one element might
    possibly be matched to the other.  Larger quality values are better, the
    maximal similarity is one.  Let \f$\Delta_\textit{RT}\f$ and
    \f$\Delta_\textit{MZ}\f$ be the absolute values of the RT and MZ differences
    in the data.  Then the similarity value is
    \f[
    \frac{1}{
    \big( 1 + \Delta_\textit{RT} \cdot \textit{diff\_intercept\_RT} \big)^\textit{diff\_exponent\_RT}
    \cdot
    \big( 1 + \Delta_\textit{MZ} \cdot \textit{diff\_intercept\_MZ} \big)^\textit{diff\_exponent\_MZ}
    }
    \f]

    Choosing diff_exponent: This parameter controls the growth rate of the
    penalty for differences.  It is for example possible to search for pairs
    using the absolute distance in RT (which should not be very susceptible to
    outliers) and the squared distance in MZ (where small difference occur
    frequently, but large differences indicate a mismatch).

    Choosing diff_intercept: Since we are taking the reciprocal value ("1/..."),
    we include an offset to avoid division by zero in case \f$\Delta=0\f$.  To
    set this parameter, ask yourself: How much worse is a difference of 1
    compared to no difference?

    The following image illustrates the influence of these parameters:
    \image html SimplePairFinder.png "Influence of the parameters intercept and exponent"

    @htmlinclude OpenMS_SimplePairFinder.parameters

    @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI SimplePairFinder :
    public BaseGroupFinder
  {
public:
    ///Base class
    typedef BaseGroupFinder Base;

    /// Constructor
    SimplePairFinder();
    /// Destructor
    ~SimplePairFinder() override
    {
    }

    /// returns an instance of this class
    static BaseGroupFinder * create()
    {
      return new SimplePairFinder();
    }

    /// returns the name of this module
    static const String getProductName()
    {
      return "simple";
    }

    /**
      @brief Run the algorithm

      @note Exactly two @em input maps must be provided.
      @note All two @em input maps must be provided.

      @exception Exception::IllegalArgument is thrown if the input data is not valid.
    */
    void run(const std::vector<ConsensusMap> & input_maps, ConsensusMap & result_map) override;

protected:

    //docu in base class
    void updateMembers_() override;

    /// A parameter for similarity_().
    double diff_exponent_[2];

    /// A parameter for similarity_().
    double diff_intercept_[2];

    /// Minimal pair quality
    double pair_min_quality_;

    /**@brief Compute the similarity for a pair of elements.
    */
    double similarity_(ConsensusFeature const & left, ConsensusFeature const & right) const;

  }; // SimplePairFinder

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_SimplePairFinder_H

/*

gnuplot history - how the plot was created - please do not delete this receipt

f(x,intercept,exponent)=1/(1+(abs(x)*intercept)**exponent)
set terminal postscript enhanced color
set output "choosingsimplepairfinderparams.ps"
set size ratio .3
plot [-3:3] [0:1] f(x,1,1), f(x,2,1), f(x,1,2), f(x,2,2)

*/
