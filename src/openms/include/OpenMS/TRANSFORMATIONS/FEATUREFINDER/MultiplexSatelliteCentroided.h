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

#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief data structure storing a single satellite peak
   *
   * The satellite peak is part of a centroided MSExperiment.
   * Hence indices rt_idx_ and mz_idx_ are sufficient to specify RT, m/z and intensity.
   * 
   * @see MultiplexFilteredPeak, MultiplexSatelliteProfile
   */
  class OPENMS_DLLAPI MultiplexSatelliteCentroided
  {
    public:

    /**
     * @brief constructor
     */
    MultiplexSatelliteCentroided(size_t rt_idx, size_t mz_idx);
    
    /**
     * @brief returns the m/z index of the satellite peak
     */
    size_t getMZidx() const;
     
    /**
     * @brief returns the RT index of the satellite peak
     */
    size_t getRTidx() const;
    
    private:
     
    /**
     * @brief indices of the satellite peak position in the centroided experiment
     * 
     * Spectral index and peak index within the spectrum of the satellite peak.
     */
    size_t rt_idx_;
    size_t mz_idx_;
    
  };
  
}

