// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once


#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <OpenMS/OpenMSConfig.h>

#include <map>
#include <vector>

namespace OpenMS
{
  class String;
}

namespace OpenSwath
{

  /**
    @brief Quality Control function for OpenSwath

    All functions are static, so there is no need to instanciate the class.

  */
  class OPENMS_DLLAPI SwathQC
  {
  public:
    typedef std::map<int,int> ChargeDistribution;
    SwathQC() = delete;

    /**
      @brief Sample the spectra in all Swath maps of a particular ms-level (usually MS1) and determine all charge states.

      From all swath maps in @p swath_maps which match the desired @p ms_level we extract @p nr_samples spectra (subsampling),
      determine the charge states of all isotopic envelopes and return their total counts.

      PeakPicking is performed internally if the data is estimated to be profile data.
      Uses Deisotoper for charge determination.

      @param [IN] swath_maps Swath maps of mixed ms-level
      @param [IN] ms_level   Filter swath maps by this ms level (allowed is 1, 2; Exception otherwise)
      @param [IN] nr_samples Number of spectra to sample per Swath map. To sample all spectra, set to -1
      @param [IN] mz_tol     Error tolerance in Th in which an isotopic peak is expected (assuming C12-C13 distance)
      @return Distribution of charge (key = charge, value = counts)

      @throw Exception::InvalidParameter if ms_level <1 or >2
      @throw Exception::Postcondition if Deisotoper did not return charge data

    */
    static ChargeDistribution getChargeDistribution(const std::vector<SwathMap>& swath_maps, const int level, const size_t nr_samples, const double mz_tol);


    static void storeJSON(const OpenMS::String& filename, const ChargeDistribution& cd);

  };
    
}

