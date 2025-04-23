// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>

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

    /// return consensus spectrum from last function call operator
    const PeakSpectrum & lastconsensus() const;

    ///
    std::map<UInt, UInt> getPeakMap() const;

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
    mutable std::map<UInt, UInt> peak_map_;
  };

}
