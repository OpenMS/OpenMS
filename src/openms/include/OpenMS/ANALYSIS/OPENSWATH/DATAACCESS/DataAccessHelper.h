// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_DATAACCESSHELPER_H
#define OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_DATAACCESSHELPER_H

#include <fstream>

#include <boost/numeric/conversion/cast.hpp>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

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
    static void convertToOpenMSSpectrum(const OpenSwath::SpectrumPtr sptr, OpenMS::MSSpectrum<> & spectrum);

    /// Convert an OpenMS Spectrum to an SpectrumPtr
    static OpenSwath::SpectrumPtr convertToSpectrumPtr(const OpenMS::MSSpectrum<> & spectrum);

    /// Convert a ChromatogramPtr to an OpenMS Chromatogram
    static void convertToOpenMSChromatogram(OpenMS::MSChromatogram<> & chromatogram, const OpenSwath::ChromatogramPtr cptr);

    /// convert from the OpenMS TargetedExperiment to the LightTargetedExperiment
    static void convertTargetedExp(const OpenMS::TargetedExperiment & transition_exp_, OpenSwath::LightTargetedExperiment & transition_exp);

    /// convert from the OpenMS TargetedExperiment Peptide to the LightTargetedExperiment Peptide
    static void convertTargetedPeptide(const TargetedExperiment::Peptide& pep, OpenSwath::LightPeptide & p);

    /// convert from the LightPeptide to an OpenMS AASequence (with correct modifications)
    static void convertPeptideToAASequence(const OpenSwath::LightPeptide & peptide, AASequence & aa_sequence);

  };

} //end namespace OpenMS

#endif
