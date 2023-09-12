// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

namespace OpenMS
{
  /**
    @brief Implements the isotope wavelet feature finder.

    The FeatureFinderAlgorithmIsotopeWavelet class has been designed for finding features in 1D
    or 2D MS data sets using the isotope wavelet. In the case of two dimensional data, the class
    provides additionally the sweep line algorithm. Please note that the algorithm implemented
    here is only marginally related to the algorithm presented in
    Schulz-Trieglaff et al. (2007, 2008), as no fitting procedure is applied anymore after the
    wavelet-based seeding step. The wavelet has been designed to extract even very lowly-abundant
    features (see Hussong et al. (2007, 2009)), usually featuring a very low
    signal-to-noise ratio. The wavelet in its current implementation is not able to resolve
    overlapping patterns (see also Hussong et al. (2009)) and slightly shifts masses to the right
    due to the construction of the wavelet.

    @htmlinclude OpenMS_FeatureFinderAlgorithmIsotopeWavelet.parameters

    @ingroup FeatureFinder
  */
  class OPENMS_DLLAPI FeatureFinderAlgorithmIsotopeWavelet :
    public FeatureFinderAlgorithm
  {
public:

    typedef FeatureFinderAlgorithm Base;
    typedef Peak1D PeakType;

    /** @brief Default Constructor */
    FeatureFinderAlgorithmIsotopeWavelet();

    /** @brief Destructor. */
    ~FeatureFinderAlgorithmIsotopeWavelet() override;

    MSSpectrum* createHRData(const UInt i);

    /** @brief The working horse of this class. */
    void run() override;

    static const String getProductName();

    static FeatureFinderAlgorithm* create();

protected:

    /** @brief Internally used data structure for the sweep line algorithm. */
    struct BoxElement
    {
      double mz;
      UInt c; ///<Note, this is not the charge (it is charge-1!!!)
      double score;
      double intens;
      double RT; ///<The elution time (not the scan index)
    };

    typedef std::map<UInt, BoxElement> Box; ///<Key: RT (index), value: BoxElement

    UInt max_charge_; ///<The maximal absolute charge state we will consider
    double intensity_threshold_; ///<The only parameter of the isotope wavelet
    UInt RT_votes_cutoff_, real_RT_votes_cutoff_, RT_interleave_; ///<The number of subsequent scans a pattern must cover in order to be considered as signal
    String use_gpus_, intensity_type_;
    bool check_PPMs_, hr_data_;
    std::vector<UInt> gpu_ids_; ///< A list of all GPU devices that can be used

    Int progress_counter_;

    void updateMembers_() override;

  };

} //namespace

