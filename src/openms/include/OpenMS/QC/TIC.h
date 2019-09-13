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
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>
#include <OpenMS/KERNEL/MSExperiment.h>

/**
 * @brief Total Ion Count (TIC) as a QC metric
 *
 * Simple class to calculate the TIC of an MSExperiment.
 * Allows for multiple usage, because each calculated TIC is
 * stored internally. Those results can then be returned using
 * getResults().
 *
 */

namespace OpenMS
{
  class OPENMS_DLLAPI TIC : public QCBase
  {
  public:
    /// Constructor
    TIC() = default;

    /// Destructor
    virtual ~TIC() = default;
    void clear();

    /**
    @brief Compute Total Ion Count and applies the resampling algorithm, if a bin size in RT seconds greater than 0 is given.

    All MS1 TICs within a bin are summed up.

    @param exp Peak map to compute the MS1 tick from
    @param bin_size RT bin size in seconds
    @return TIC Chromatogram
    **/
    void compute(const MSExperiment &exp, float bin_size=0);

    const String& getName() const override;

    const std::vector<MSChromatogram>& getResults() const ;

    QCBase::Status requires() const override;

  private:
    const String name_ = "TIC";
    std::vector<MSChromatogram> results_;
  };
}
