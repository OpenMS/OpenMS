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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ELUTIONMODELFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ELUTIONMODELFITTER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

namespace OpenMS
{
  /**
     @brief Helper class for fitting elution models to features

     @htmlinclude OpenMS_ElutionModelFitter.parameters
  */
  class OPENMS_DLLAPI ElutionModelFitter :
    public DefaultParamHandler
  {

  public:
    /// Default constructor
    ElutionModelFitter();

    /// Destructor
    ~ElutionModelFitter() override;

    /**
       @brief Fit models of elution profiles to all features (and validate them)

       Assumptions (not checked!):
       - all features have meta values "left-"/"rightWidth" giving RT start/end
       - all features have subordinates (for the mass traces/transitions)
       - each subordinate has an appropriate meta value "isotope_probability"
       - each subordinate has one convex hull
       - all convex hulls in one feature contain the same number (> 0) of points
       - the y coordinates of the hull points store the intensities
    */
    void fitElutionModels(FeatureMap& features);

  protected:
    typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
    typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

    /// Calculate quality of model fit (mean relative error)
    double calculateFitQuality_(const TraceFitter* fitter, 
                                const MassTraces& traces);
  };
}

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ELUTIONMODELFITTER_H
