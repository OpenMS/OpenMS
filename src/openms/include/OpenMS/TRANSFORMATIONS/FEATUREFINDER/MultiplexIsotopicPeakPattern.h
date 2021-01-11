// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief data structure for pattern of isotopic peaks
   * 
   * Groups of peptides appear as characteristic patterns of isotopic peaks
   * in MS1 spectra. For example, for an Arg6 labeled SILAC peptide pair
   * of charge 2+ with three isotopic peaks we expect peaks
   * at relative m/z shifts of 0, 0.5, 1, 3, 3.5 and 4 Th.
   */
  class OPENMS_DLLAPI MultiplexIsotopicPeakPattern
  {
    public:

    /**
     * @brief constructor
     */
    MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi);
    
    /**
     * @brief returns charge
     */
    int getCharge() const;
     
    /**
     * @brief returns peaks per peptide
     */
    int getPeaksPerPeptide() const;
    
    /**
     * @brief returns mass shifts
     */
    MultiplexDeltaMasses getMassShifts() const;
    
    /**
     * @brief returns mass shift index
     */
    int getMassShiftIndex() const;
    
    /**
     * @brief returns number of mass shifts i.e. the number of peptides in the multiplet
     */
    unsigned getMassShiftCount() const;
   
    /**
     * @brief returns mass shift at position i
     */
    double getMassShiftAt(size_t i) const;
    
    /**
     * @brief returns m/z shift at position i
     */
    double getMZShiftAt(size_t i) const;
    
    /**
     * @brief returns number of m/z shifts
     */
    unsigned getMZShiftCount() const;
    
    private:

    /**
     * @brief m/z shifts between isotopic peaks
     * (number of mz_shifts_ = peaks_per_peptide_ * number of mass_shifts_)
     */
    std::vector<double> mz_shifts_;

    /**
     * @brief charge
     */
    int charge_;

    /**
     * @brief number of isotopic peaks in each peptide
     */
    int peaks_per_peptide_;
   
    /**
     * @brief mass shifts between peptides
     * (including zero mass shift for first peptide)
     */
    MultiplexDeltaMasses mass_shifts_;

    /**
     * @brief index in mass shift list
     */
    int mass_shift_index_;
      
 };
  
}

