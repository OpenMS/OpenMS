// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  /**
    @brief SpectrumCheapDPCorr calculates an optimal alignment on stick spectra

    to save computing time only Peak Pairs that could get a score > 0 are considered
    which Peak Pairs could get scores > 0 ? <br>
    Peaks get a score depending on the difference in position and the heights of the peaks <br>
    pairs with positions that differ more than some limit get score 0

    @htmlinclude OpenMS_SpectrumCheapDPCorr.parameters

    @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI SpectrumCheapDPCorr :
    public PeakSpectrumCompareFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    SpectrumCheapDPCorr();

    /// copy constructor
    SpectrumCheapDPCorr(const SpectrumCheapDPCorr & source);

    /// destructor
    ~SpectrumCheapDPCorr() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    SpectrumCheapDPCorr & operator=(const SpectrumCheapDPCorr & source);

    double operator()(const PeakSpectrum & a, const PeakSpectrum & b) const override;

    double operator()(const PeakSpectrum & a) const override;
    // @}

    // @name Accessors
    // @{
    ///
    static PeakSpectrumCompareFunctor * create() { return new SpectrumCheapDPCorr(); }

    ///
    static const String getProductName()
    {
      return "SpectrumCheapDPCorr";
    }

    /// return consensus spectrum from last function call operator
    const PeakSpectrum & lastconsensus() const;

    ///
    Map<UInt, UInt> getPeakMap() const;

    /// set weighting of the second spectrum for consensus from next function call operator
    void setFactor(double f);
    // @}

private:

    /// O(n^2) dynamical programming
    double dynprog_(const PeakSpectrum &, const PeakSpectrum &, int, int, int, int) const;

    /// similarity of two peaks
    double comparepeaks_(double posa, double posb, double inta, double intb) const;

    static const String info_;

    /// consensus spectrum of the last comparison
    mutable PeakSpectrum lastconsensus_;

    /// should peaks with no alignment partner be kept in the consensus?
    bool keeppeaks_;

    /// weighting factor for the next consensus spectrum
    mutable double factor_;

    /// last peak map
    mutable Map<UInt, UInt> peak_map_;
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMCHEAPDPCORR_H
