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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>

namespace OpenMS
{

/**
 *  @brief An implementation of the Morpheus PSM scoring function
 *  Inspired by a C# implementation by C. Wenger released under MIT license
 */               
struct OPENMS_DLLAPI MorpheusScore
{

  /// score and subscores
  struct OPENMS_DLLAPI Result
  {
    Size matches = 0; //< matched theoretical peaks
    Size n_peaks = 0; //< number of theoretical peaks
    float score = 0; //< Morpheus score (matched peaks + matched ion current / TIC)
    float MIC = 0; //< ion current of matches (experimental peaks)
    float TIC = 0; //< total ion current (experimental peak) 
    float err = 0; //< average absolute mass error of matched fragments (in Da)
  };

  /// returns Morpheus Score, #matched ions, #total ions, #matched intensities, #total fragment intensities (TIC)
  static Result compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const PeakSpectrum& theo_spectrum);
};

}


