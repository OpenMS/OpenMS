// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSOPENMS_H
#define OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSOPENMS_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /**
    @brief An implementation of the OpenSWATH Spectrum Access interface using OpenMS

  */
  class OPENMS_DLLAPI SpectrumAccessOpenMS :
    public OpenSwath::ISpectrumAccess
  {
public:
    typedef OpenMS::MSExperiment<Peak1D> MSExperimentType;
    typedef OpenMS::MSSpectrum<Peak1D> MSSpectrumType;
    typedef OpenMS::MSChromatogram<ChromatogramPeak> MSChromatogramType;

    /// Constructor
    explicit SpectrumAccessOpenMS(boost::shared_ptr<MSExperimentType> ms_experiment);

    /// Destructor
    ~SpectrumAccessOpenMS();

    /// Copy constructor
    SpectrumAccessOpenMS(const SpectrumAccessOpenMS & rhs);

    /// Light clone operator (actual data will not get copied)
    boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const;

    OpenSwath::SpectrumPtr getSpectrumById(int id);

    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const;

    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const;

    size_t getNrSpectra() const;

    SpectrumSettings getSpectraMetaInfo(int id) const;

    OpenSwath::ChromatogramPtr getChromatogramById(int id);

    // FEATURE ?
    // ChromatogramPtr getChromatogramByPrecursorMZ(double mz, double deltaMZ);

    size_t getNrChromatograms() const;

    ChromatogramSettings getChromatogramMetaInfo(int id) const;

    std::string getChromatogramNativeID(int id) const;

private:
    boost::shared_ptr<MSExperimentType> ms_experiment_;

  };
} //end namespace OpenMS

#endif
