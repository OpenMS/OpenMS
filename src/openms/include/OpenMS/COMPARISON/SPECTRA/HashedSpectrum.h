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

#include <vector>

namespace OpenMS
{
/**
 * @brief Auxiliary data structure for the quick access of peaks in MSSpectrum<Peak1D>
 *
 * Finding a peak with a particular m/z in MSSpectrum<Peak1D> can be slow. This data structure
 * provides a hash index (bins_) for quick access. Given an m/z of interest, one can
 * quickly jump in the relevant neighbouthood i.e. bin(s) where such a peak might be located.
 * The getPeak(double mz) method returns a pointer to the peak in MSSpectrum<Peak1D>.
 *
 * @see MSSpectrum<Peak1D>
 */
class OPENMS_DLLAPI HashedSpectrum
{
  public:
  
    /**
     * @brief m/z interval containing pointers to the first and last data points
     * If there are no data points in this interval, both pointers are null.
     * Boundaries of the interval implicitly given by mz_min, mz_bin, mz_unit_ppm_
     * and the hash table index.
     */
    struct MzInterval
    {
		MSSpectrum<Peak1D>::Iterator begin;
		MSSpectrum<Peak1D>::Iterator end;
	};
    
    /**
     * @brief constructor taking an MSSpectrum
     * 
     * @param spectrum    MS spectrum to be indexed
     * @param mz_bin    size of the m/z bins aka intervals
     * @param mz_unit_ppm    unit for mz_bin 
     */
    HashedSpectrum(MSSpectrum<Peak1D>& spectrum, double mz_bin, bool mz_unit_ppm);

    /**
     * @brief destructor
     */
    ~HashedSpectrum();

    /**
     * @brief returns the m/z bin size
     */
    double getMzBin() const;
    
    /**
     * @brief returns whether m/z bin unit ppm or Th
     */
    bool getMzUnitPpm() const;
    
    /**
     * @brief find the peak closest to m/z
     * 
     * @param mz    m/z position
     * @param mz_tolerance    m/z tolerance between mz and peak
     * @param mz_unit_ppm    unit for mz_tolerance (ppm or Th)
     * @return pointer to the peak closest to mz or NULL if difference between mz and closest peak exceeds the tolerance
     */
    MSSpectrum<Peak1D>::Iterator findNearest(const double mz, const double mz_tolerance, const bool mz_unit_ppm) const;
    
  private:

    /// hide default C'tor
    HashedSpectrum();
    
    /**
     * @brief m/z bin size
     */
    double mz_bin_;
    
    /**
     * @brief unit for mz_bin_ (ppm or Th)
     */
    bool mz_unit_ppm_;
        
    /**
     * @brief minimal m/z in the spectrum
     */
    double mz_min_;
        
    /**
     * @brief MS spectrum to be indexed
     */
    MSSpectrum<Peak1D>& spectrum_;

    /**
     * @brief vector of all m/z intervals (hash table)
     */
    std::vector<MzInterval> bins_;

    /**
     * @brief return index of the interval containing m/z
     */
    int getIndex_(const double mz) const;

    /**
     * @brief return distance between two m/z in either Th or ppm (depending on mz_unit_ppm_)
     */
    double getDistance_(const double mz1, const double mz2) const;

};

}

#endif /* HASHEDSPECTRUM_H_ */
