// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <fstream>

#include <boost/numeric/conversion/cast.hpp>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

namespace OpenMS
{
  /**
    @brief Several helpers to convert OpenMS datastructures to structures that
           implement the OpenSWATH interfaces.
  */
  class OPENMS_DLLAPI OpenSwathDataAccessHelper
  {
public:
    /// Convert a SpectrumPtr to an OpenMS Spectrum
    static void convertToOpenMSSpectrum(const OpenSwath::SpectrumPtr& sptr, OpenMS::MSSpectrum & spectrum);

    /// Convert an OpenMS Spectrum to an SpectrumPtr
    static OpenSwath::SpectrumPtr convertToSpectrumPtr(const OpenMS::MSSpectrum & spectrum);

    /// Convert a ChromatogramPtr to an OpenMS Chromatogram
    static void convertToOpenMSChromatogram(const OpenSwath::ChromatogramPtr& cptr, OpenMS::MSChromatogram & chromatogram);

    /// Convert an OpenMS Chromatogram to an ChromatogramPtr
    static OpenSwath::ChromatogramPtr convertToChromatogramPtr(const OpenMS::MSChromatogram & chromatogram);

    static void convertToOpenMSChromatogramFilter(OpenMS::MSChromatogram & chromatogram,
                                                  const OpenSwath::ChromatogramPtr& cptr,
                                                  double rt_min,
                                                  double rt_max);

    /// convert from the OpenMS TargetedExperiment to the LightTargetedExperiment
    static void convertTargetedExp(const OpenMS::TargetedExperiment & transition_exp_, OpenSwath::LightTargetedExperiment & transition_exp);

    /// convert from the OpenMS TargetedExperiment Peptide to the LightTargetedExperiment Peptide
    static void convertTargetedCompound(const TargetedExperiment::Peptide& pep, OpenSwath::LightCompound& comp);

    /// convert from the OpenMS TargetedExperiment Compound to the LightTargetedExperiment Compound
    static void convertTargetedCompound(const TargetedExperiment::Compound& compound, OpenSwath::LightCompound& comp);

    /// convert from the LightCompound to an OpenMS AASequence (with correct modifications)
    static void convertPeptideToAASequence(const OpenSwath::LightCompound & peptide, AASequence & aa_sequence);

  };

} //end namespace OpenMS


