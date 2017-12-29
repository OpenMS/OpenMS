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

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h>

namespace OpenMS
{

    /// Constructor
  SpectrumAccessSqMass::SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler) :
      handler_(handler)
    {}

    SpectrumAccessSqMass::SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler, std::vector<int> indices) :
      handler_(handler),
      sidx_(indices)
    {}


    SpectrumAccessSqMass::SpectrumAccessSqMass(SpectrumAccessSqMass sp, std::vector<int> indices) :
      handler_(sp.handler_)
    {
      if (indices.empty())
      {
        sidx_ = sp.sidx_;
      }
      else if (sp.sidx_.empty())
      {
        sidx_ = indices;
      }
      else
      {
        // we only want to select a subset of the currently selected indices
        for (Size k = 0; k < indices.size(); k++)
        {
          sidx_.push_back( sp.sidx_[ indices[k] ] );
        }
      }
    }

    /// Destructor
    SpectrumAccessSqMass::~SpectrumAccessSqMass() {}

    /// Copy constructor
    SpectrumAccessSqMass::SpectrumAccessSqMass(const SpectrumAccessSqMass & rhs) :
      handler_(rhs.handler_),
      sidx_(rhs.sidx_)
    {
    }

    /// Light clone operator (actual data will not get copied)
    boost::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessSqMass::lightClone() const
    {
      return boost::shared_ptr<SpectrumAccessSqMass>(new SpectrumAccessSqMass(*this));
    }

    OpenSwath::SpectrumPtr SpectrumAccessSqMass::getSpectrumById(int id)
    {
      std::vector<int> indices;
      if (sidx_.empty())
      {
        indices.push_back(id);
      }
      else
      {
        indices.push_back(sidx_[id]);
      }

      // read MSSpectra and prepare for conversion
      std::vector<MSSpectrum> tmp_spectra;
      handler_.readSpectra(tmp_spectra, indices, false);

      const MSSpectrumType& spectrum = tmp_spectra[0];
      OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
      OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
      for (MSSpectrumType::const_iterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        mz_array->data.push_back(it->getMZ());
        intensity_array->data.push_back(it->getIntensity());
      }

      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      sptr->setMZArray(mz_array);
      sptr->setIntensityArray(intensity_array);
      return sptr;
    }

    OpenSwath::SpectrumMeta SpectrumAccessSqMass::getSpectrumMetaById(int id) const
    {
      std::vector<int> indices;
      if (sidx_.empty())
      {
        indices.push_back(id);
      }
      else
      {
        indices.push_back(sidx_[id]);
      }

      // read MSSpectra and prepare for conversion
      std::vector<MSSpectrum> tmp_spectra;
      handler_.readSpectra(tmp_spectra, indices, false);

      const MSSpectrumType& spectrum = tmp_spectra[0];
      OpenSwath::SpectrumMeta m;
      m.id = spectrum.getNativeID();
      m.RT = spectrum.getRT();
      m.ms_level = spectrum.getMSLevel();
      return m;
    }

    void SpectrumAccessSqMass::getAllSpectra(std::vector< OpenSwath::SpectrumPtr > & spectra, std::vector< OpenSwath::SpectrumMeta > & spectra_meta) const
    {
      // read MSSpectra and prepare for conversion
      std::vector<MSSpectrum> tmp_spectra;

      if (sidx_.empty())
      {
        MSExperiment exp;
        {
          handler_.readExperiment(exp, false);
        }

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

    std::vector<std::size_t> SpectrumAccessSqMass::getSpectraByRT(double RT, double deltaRT) const
    {
      OPENMS_PRECONDITION(deltaRT >= 0, "Delta RT needs to be a positive number");
      std::cout << "std::vector<std::size_t> SpectrumAccessSqMass::getSpectraByRT(double RT, double deltaRT) const " << std::endl;
      std::vector<std::size_t> res = handler_.getSpectraIndicesbyRT(RT, deltaRT, sidx_);

      if (sidx_.empty())
      {
        return res;
      }
      else
      {
        // we need to map the resulting indices back to the external indices
        std::vector<std::size_t> res_mapped;
        for (Size k = 0; k > res.size(); k++)
        {
          for (Size s_it = 0; s_it > sidx_.size(); s_it++)
          {
            if (res[k] == (size_t)sidx_[s_it]) {res_mapped.push_back(s_it);}
          }
        }
        return res_mapped;
      }
    }

    size_t SpectrumAccessSqMass::getNrSpectra() const
    {
      size_t res;
      if (sidx_.empty())
      {
        res = handler_.getNrSpectra();
      }
      else
      {
        res = sidx_.size();
      }
      return res;
    }

    OpenSwath::ChromatogramPtr SpectrumAccessSqMass::getChromatogramById(int /* id */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

    size_t SpectrumAccessSqMass::getNrChromatograms() const
    {
      size_t res;
      // TODO: currently chrom indices are not supported
      res = handler_.getNrChromatograms();
      return res;
    }

    std::string SpectrumAccessSqMass::getChromatogramNativeID(int /* id */) const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }

} //end namespace OpenMS
