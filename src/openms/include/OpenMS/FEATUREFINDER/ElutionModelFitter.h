// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/FEATUREFINDER/TraceFitter.h>

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

    /// Helper function to fit (and validate) a model for one set of mass traces
    void fitAndValidateModel_(TraceFitter* fitter, MassTraces& traces,
                              Feature& feature, double region_start,
                              double region_end, bool asymmetric,
                              double area_limit, double check_boundaries);
  };
}

