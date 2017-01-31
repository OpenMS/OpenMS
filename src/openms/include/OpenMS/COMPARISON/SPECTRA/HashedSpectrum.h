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
     * 
     * @param raw_spectrum    
     * @param mz_bin
     * @param mz_tolerance
     */
    HashedSpectrum(MSSpectrum<Peak1D>& raw_spectrum, double mz_bin, double mz_tolerance, bool mz_unit_ppm);

    /**
     * @brief destructor
     */
    ~HashedSpectrum();

    /**
     * @brief returns the m/z bin size
     */
    double getMzBin() const;
    
    /**
     * @brief returns the m/z tolerance
     */
    double getMzTolerance() const;
    
    /**
     * @brief returns whether m/z bin unit ppm or Th
     */
    bool getMzUnitPpm() const;
    
    /**
     * @brief returns pointer to the peak closest to m/z
     * 
     * @param mz    m/z position
     * @return pointer to the peak closest to mz or NULL if difference between mz and closest peak m/z exceeds mz_tolerance_
     */
    MSSpectrum<Peak1D>::pointer getPeak(double mz);
    
  private:

    /// hide default C'tor
    HashedSpectrum();
    
    /**
     * @brief m/z bin size
     */
    double mz_bin_;
    
    /**
     * @brief m/z tolerance
     */
    double mz_tolerance_;
    
    /**
     * @brief whether unit for mz_bin_ and mz_tolerance_ is ppm or Th
     */
    bool mz_unit_ppm_;
        
    /**
     * @brief minimal m/z in the spectrum
     */
    double mz_min_;
        
    /**
     * @brief vector of all m/z intervals (hash table)
     */
    std::vector<MzInterval> intervals_;

    /**
     * @brief return index of the interval containing m/z
     */
    int getIndex_(double mz);

    /**
     * @brief return distance between two m/z in either Th or ppm (depending on mz_unit_ppm_)
     */
    double getDistance_(double mz1, double mz2);

};

}

#endif /* HASHEDSPECTRUM_H_ */
