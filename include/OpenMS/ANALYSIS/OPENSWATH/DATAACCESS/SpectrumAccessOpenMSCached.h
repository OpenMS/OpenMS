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

#ifndef OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSOPENMSCACHED_H
#define OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSOPENMSCACHED_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h>

namespace OpenMS
{

  /**
    @brief An implementation of the Spectrum Access interface using on-disk caching

    This class implements the OpenSWATH Spectrum Access interface
    (ISpectrumAccess) using the CachedmzML class which is able to read and
    write a cached mzML file.

    @note This implementation is @a not thread-safe since it keeps internally a
    single file access pointer which it moves when accessing a specific
    data item. The caller is responsible to ensure that access is performed
    atomically.

  */
  class OPENMS_DLLAPI SpectrumAccessOpenMSCached :
    public OpenSwath::ISpectrumAccess
  {

public:
    typedef OpenMS::MSExperiment<Peak1D> MSExperimentType;
    typedef OpenMS::MSSpectrum<Peak1D> MSSpectrumType;

    /**
      @brief Constructor, opens the file stream

      @param filename The filename of the .mzML file (it is assumed a second
      file .mzML.cached exists).

      @throws Exception::FileNotFound is thrown if the file is not found
      @throws Exception::ParseError is thrown if the file cannot be parsed
    */
    explicit SpectrumAccessOpenMSCached(String filename);

    /**
      @brief Destructor
    */
    ~SpectrumAccessOpenMSCached();

    OpenSwath::SpectrumPtr getSpectrumById(int id);

    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const;

    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const;

    size_t getNrSpectra() const;

    SpectrumSettings getSpectraMetaInfo(int id) const;

    OpenSwath::ChromatogramPtr getChromatogramById(int id);

    size_t getNrChromatograms() const;

    ChromatogramSettings getChromatogramMetaInfo(int id) const;

    std::string getChromatogramNativeID(int id) const;

private:

    /// Meta data
    MSExperimentType meta_ms_experiment_;

    /// Internal filestream 
    std::ifstream ifs_;

    /// Name of the mzML file
    String filename_;

    /// Name of the cached mzML file
    String filename_cached_;

    /// Indices
    std::vector<Size> spectra_index_;
    std::vector<Size> chrom_index_;
  };

} //end namespace

#endif
