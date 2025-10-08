// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>

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

