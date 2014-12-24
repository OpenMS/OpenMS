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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DECHARGING_ILPDCWRAPPER_H
#define OPENMS_ANALYSIS_DECHARGING_ILPDCWRAPPER_H

#include <vector>
#include <set>

#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{

  class MassExplainer;
  class FeatureMap;

  class OPENMS_DLLAPI ILPDCWrapper
  {

public:
    typedef std::vector<ChargePair> PairsType;
    typedef PairsType::size_type PairsIndex;

    ///Constructor
    ILPDCWrapper();

    ///Destructor
    virtual ~ILPDCWrapper();

    /// Compute optimal solution and return value of objective function
    /// If the input feature map is empty, a warning is issued and -1 is returned.
    /// @return value of objective function
    /// and @p pairs will have all realized edges set to "active"
    double compute(const FeatureMap fm, PairsType & pairs, Size verbose_level) const;

private:

    /// slicing the problem into subproblems
    double computeSlice_(const FeatureMap fm,
                             PairsType & pairs,
                             const PairsIndex margin_left,
                             const PairsIndex margin_right,
                             const Size verbose_level) const;

    /// slicing the problem into subproblems
    double computeSliceOld_(const FeatureMap fm,
                                PairsType & pairs,
                                const PairsIndex margin_left,
                                const PairsIndex margin_right,
                                const Size verbose_level) const;

    /// calculate a score for the i_th edge
    double getLogScore_(const PairsType::value_type & pair, const FeatureMap & fm) const;

    typedef Map<String, std::set<Size> > FeatureType_;

    // add another charge annotation variant for a feature
    void updateFeatureVariant_(FeatureType_ & f_set, const String & rota_l, const Size & v) const;



  };   // !class

} // !namespace

#endif // OPENMS_ANALYSIS_DECHARGING_ILPDCWRAPPER_H
