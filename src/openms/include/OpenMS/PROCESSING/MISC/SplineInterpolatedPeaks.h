// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/PROCESSING/MISC/SplinePackage.h>

#include <vector>

namespace OpenMS
{
/**
 * @brief Data structure for spline interpolation of MS1 spectra and chromatograms
 *
 * The data structure consists of a set of splines, each interpolating the MS1 spectrum (or chromatogram) in a
 * certain m/z (or RT) range. Between these splines no raw data points exist and the intensity is identical to zero.
 * 
 * A spline on non-equi-distant input data is not well supported in regions without data points. Hence, a spline tends to
 * swing wildly in these regions and cannot be used for reliable interpolation. We assume that in m/z (or RT) regions
 * without data points, the spectrum (or chromatogram) is identical to zero.
 *
 * @see SplinePackage
 * @see MSSpectrum
 * @see MSChromatogram
 */
class OPENMS_DLLAPI SplineInterpolatedPeaks
{
  public:
    /**
     * @brief constructor taking two vectors
     * (and an optional scaling factor for the m/z (or RT) step width)
     * 
     * @note Vectors are assumed to be sorted by m/z (or RT)!
     */
    SplineInterpolatedPeaks(const std::vector<double>& pos, const std::vector<double>& intensity);

    /**
     * @brief constructor taking an MSSpectrum
     * (and an optional scaling factor for the m/z step width)
     */
    SplineInterpolatedPeaks(const MSSpectrum& raw_spectrum);

    /**
     * @brief constructor taking an MSChromatogram
     * (and an optional scaling factor for the RT step width)
     */
    SplineInterpolatedPeaks(const MSChromatogram& raw_chromatogram);

    /**
     * @brief destructor
     */
    ~SplineInterpolatedPeaks();

    /**
     * @brief returns the minimum m/z (or RT) of the spectrum
     */
    double getPosMin() const;

    /**
     * @brief returns the maximum m/z (or RT) of the spectrum
     */
    double getPosMax() const;

    /** 
     * @brief Get number of spline packages found during initialization
     *
     * Note that this function should be called right after the C'tor to ensure the spectrum
     * has some usable data to work on.
     * In case there are no packages, a subsequent call to getNavigator() will throw an exception.
     */
    size_t size() const;

    /**
    * @brief iterator class for access of spline packages
    */
    class OPENMS_DLLAPI Navigator
    {
      public:
        /**
        * @brief constructor of iterator
        * 
        * @param packages Spline packages to be accessed
        * @param pos_max Maximum in m/z (or RT) of the spectrum (or chromatogram)
        * @param scaling    The step width can be scaled by this factor. Often it is advantageous to iterate
        * in slightly smaller steps over the spectrum (or chromatogram).
        */
        Navigator(const std::vector<SplinePackage>* packages, double pos_max, double scaling);

        /**
        * @brief constructor (for pyOpenMS)
        */
        Navigator();

        /**
        * @brief destructor
        */
        ~Navigator();

        /**
        * @brief returns spline interpolated intensity at this position
        * (fast access since we can start search from lastPackage)
        */
        double eval(double pos);

        /**
        * @brief returns the next sensible m/z (or RT) position for scanning through a spectrum (or chromatogram)
        * (fast access since we can start search from lastPackage)
        * 
        * In the middle of a package, we increase the position by the average spacing of the input data (times a scaling factor).
        * At the end of a package, we jump straight to the beginning of the next package.
        */
        double getNextPos(double pos);

      private:
        
        /**
        * @brief list of spline packages to be accessed
        */
        const std::vector<SplinePackage> * packages_;

        /**
        * @brief index of spline package last accessed
        */
        size_t last_package_;

        /**
        * @brief m/z (or RT) limits of the spectrum (or chromatogram)
        */
        double pos_max_;
        
        /**
        * @brief scaling of the step width
        * 
        * Each package stores its own step width, which is the average spacing of the input data points.
        * This step width can be adjusted by the scaling factor. Often it is advantageous to use a step width
        * which is somewhat smaller than the average raw data spacing.
        * 
        * @see getNextPos() 
        */
        double pos_step_width_scaling_;
    };

    /**
    * @brief returns an iterator for access of spline packages
    *
    * Will throw an exception if no packages were found during construction.
    * Check using getSplineCount().
    *
    * Make sure that the underlying SplineInterpolatedPeaks does not run out-of-scope since the
    * Navigator relies on its data.
    * 
    * @param scaling    step width scaling parameter
    *
    * @throw Exception::InvalidSize if packages is empty
    */
    SplineInterpolatedPeaks::Navigator getNavigator(double scaling = 0.7);

  private:

    /// hide default C'tor
    SplineInterpolatedPeaks();
    
    /**
     * @brief m/z (or RT) limits of the spectrum
     */
    double pos_min_;
    double pos_max_;

    /**
     * @brief set of spline packages each interpolating in a certain m/z (or RT) range
     */
    std::vector<SplinePackage> packages_;

    /**
     * @brief section common for all constructors
     */
    void init_(const std::vector<double>& pos, const std::vector<double>& intensity);


};

}

