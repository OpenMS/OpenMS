// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_HASHEDSPECTRUM_H
#define OPENMS_FILTERING_DATAREDUCTION_HASHEDSPECTRUM_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/MATH/MISC/Spline2d.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
/**
 * @brief Data structure for spline interpolation of MS1 spectra
 *
 * The data structure consists of a set of splines, each interpolating the MS1 spectrum in a certain m/z range.
 * Between these splines no raw data points exist and the spectrum intensity is identical to zero.
 *
 * @see SplinePackage
 * @see MSSpectrum<Peak1D>
 */
class OPENMS_DLLAPI HashedSpectrum
{
  public:
  
    /**
     * @brief m/z interval containing pointers to the first and last data points
     * If there are no data points in this interval, both pointers are null.
     * Boundaries of the interval implicitly given by mz_min, mz_bin and hash table index.
     */
    struct MzInterval
    {
		MSSpectrum<Peak1D>::pointer first;
		MSSpectrum<Peak1D>::pointer last;
	};
  
    /**
     * @brief constructor taking an MSSpectrum
     */
    HashedSpectrum(MSSpectrum<Peak1D>& raw_spectrum);

    /**
     * @brief destructor
     */
    ~HashedSpectrum();

    /**
     * @brief returns the minimum m/z of the spectrum
     */
    double getMzMin() const;

    /**
     * @brief returns the maximum m/z of the spectrum
     */
    double getMzMax() const;

    /**
     * @brief returns the m/z bin size
     */
    double getMzBin() const;
    
  private:

    /// hide default C'tor
    HashedSpectrum();
    
    /**
     * @brief m/z limits of the spectrum
     */
    double mz_min_;
    double mz_max_;

    /**
     * @brief m/z bin size [Th]
     */
    double mz_bin_;
        
    /**
     * @brief vector of all m/z intervals (hash table)
     */
    std::vector<MzInterval> intervals_;

};

}

#endif /* HASHEDSPECTRUM_H_ */
