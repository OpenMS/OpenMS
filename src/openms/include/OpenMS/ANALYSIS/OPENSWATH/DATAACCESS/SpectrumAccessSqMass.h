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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSSQMASS_H
#define OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSSQMASS_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <boost/shared_ptr.hpp>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort

namespace OpenMS
{
  /**
   * @brief An implementation of the Spectrum Access interface using SQL files
   *
   * The interface takes an MzMLSqliteHandler object to access spectra and
   * chromatograms from a sqlite file (sqMass). Currently access to individual
   * spectra and chromatograms are not supported due to large overhead of
   * opening a DB connection and performing a single query.
   *
   * Instead, the users should use getAllSpectra which returns all available
   * spectra together.
   *
   * The interface allows to be constructed in a way as to only provide access
   * to a subset of spectra / chromatograms by supplying a set of indices which
   * are then used to provide a transparent interface to any consumer who will
   * get access to the described subset of spectra / chromatograms. This can be
   * useful to provide a specific interface to MS1 or MS2  spectra only or to
   * different DIA / SWATH-MS windows.
   *
  */
  class OPENMS_DLLAPI SpectrumAccessSqMass :
    public OpenSwath::ISpectrumAccess
  {

public:
    typedef OpenMS::MSSpectrum<Peak1D> MSSpectrumType;
    typedef OpenMS::MSChromatogram<ChromatogramPeak> MSChromatogramType;

    /// Constructor
    SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler) :
      handler_(handler)
    {}

    SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler, std::vector<int> indices) :
      handler_(handler),
      sidx_(indices)
    {}

    /// Destructor
    virtual ~SpectrumAccessSqMass() {}

    /// Copy constructor
    SpectrumAccessSqMass(const SpectrumAccessSqMass & rhs) :
      handler_(rhs.handler_),
      sidx_(rhs.sidx_)
    {
    }

    /// Light clone operator (actual data will not get copied)
    boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const
    {
      return boost::shared_ptr<SpectrumAccessSqMass>(new SpectrumAccessSqMass(*this));
    }

    OpenSwath::SpectrumPtr getSpectrumById(int /* id */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    OpenSwath::SpectrumMeta getSpectrumMetaById(int /* id */) const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    void getAllSpectra(std::vector< OpenSwath::SpectrumPtr > & spectra, std::vector< OpenSwath::SpectrumMeta > & spectra_meta) 
    {
      // read MSSpectra and prepare for conversion
      std::vector<MSSpectrum<> > tmp_spectra;
      if (sidx_.empty())
      {
        MSExperiment exp;
        handler_.readExperiment(exp, false);
        tmp_spectra = exp.getSpectra();
      }
      else
      {
        handler_.readSpectra(tmp_spectra, sidx_, false);
      }
      spectra.reserve(tmp_spectra.size());
      spectra_meta.reserve(tmp_spectra.size());
    
      for (Size k = 0; k < tmp_spectra.size(); k++)
      {
        const MSSpectrumType& spectrum = tmp_spectra[k];
        OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
        OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
        for (MSSpectrumType::const_iterator it = spectrum.begin(); it != spectrum.end(); ++it)
        {
          mz_array->data.push_back(it->getMZ());
          intensity_array->data.push_back(it->getIntensity());
        }

        OpenSwath::SpectrumMeta m;
        m.id = spectrum.getNativeID();
        m.RT = spectrum.getRT();
        m.ms_level = spectrum.getMSLevel();
        spectra_meta.push_back(m);

        OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
        sptr->setMZArray(mz_array);
        sptr->setIntensityArray(intensity_array);
        spectra.push_back(sptr);
      }
    }

    std::vector<std::size_t> getSpectraByRT(double /* RT */, double /* deltaRT */) const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    size_t getNrSpectra() const
    {
      if (sidx_.empty())
      {
        return handler_.getNrSpectra();
      }
      else
      {
        return sidx_.size();
      }
    }

    OpenSwath::ChromatogramPtr getChromatogramById(int /* id */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    size_t getNrChromatograms() const
    {
      // TODO: currently chrom indices are not supported
      return handler_.getNrChromatograms();
    }

    std::string getChromatogramNativeID(int /* id */) const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

private:

    /// Access to underlying sqMass file
    OpenMS::Internal::MzMLSqliteHandler handler_;
    /// Optional subset of spectral indices
    std::vector<int> sidx_;
  };
} //end namespace OpenMS

#endif // OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_SPECTRUMACCESSSQMASS_H


