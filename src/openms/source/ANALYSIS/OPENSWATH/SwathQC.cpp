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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/SwathQC.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <nlohmann/json.hpp>

#include <cmath>       // round
#include <iomanip>     // setw

namespace OpenSwath
{
  using namespace OpenMS;

  SwathQC::SwathQC(const size_t cd_spectra, const double decon_ms1_mz_tol)
    : cd_(),
      nr_ms1_spectra_(0),
      cd_spectra_(cd_spectra),
      decon_ms1_mz_tol_(decon_ms1_mz_tol),
      ms1_spectra_seen_(0)
  {
  }

  std::function<void(const ExperimentalSettings&)> SwathQC::getExpSettingsFunc()
  {
    auto f = [this](const ExperimentalSettings& es)
    {
      // if member is set, we already have what we want. Besides, some parsers might call this function
      // during parse, where the information is probably missing
      if (nr_ms1_spectra_ > 0) return;

      if (!es.metaValueExists("nr_ms1_spectra"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected meta value 'nr_ms1_spectra'");
      }
      // change the member, when the lambda gets called (acts like a callback function)
      nr_ms1_spectra_ = es.getMetaValue("nr_ms1_spectra");
    };
    return f;
  }

  std::function<void(const MSSpectrum&)> SwathQC::getSpectraProcessingFunc()
  {
    auto f = [this](const MSSpectrum& spec)
    {
      // only look at MS1 spectra (for now)
      if (spec.getMSLevel() != 1) return;

      if (!isSubsampledSpectrum_(nr_ms1_spectra_, cd_spectra_, ms1_spectra_seen_))
      { 
        return;
      }

      ++ms1_spectra_seen_;

      PeakPickerHiRes pp;
      auto t = spec.getType(true);
      MSSpectrum tmp;
      if (t == MSSpectrum::SpectrumSettings::PROFILE) 
      {
        pp.pick(spec, tmp);
      }
      else if (t == MSSpectrum::SpectrumSettings::CENTROID)
      {
        tmp = spec; // make a copy, since deisotopeAndSingleCharge() will modify 
      }
      else
      {
        return; // unknown: too dangerous to analyse
      }

      if (tmp.empty())
      {
        return; // something went wrong with the spectrum after peak picking (e.g. returned empty spectrum)
      }

      // Note: this will pick up also non-peptide signals; filtering by averagine might yield better results
      Deisotoper::deisotopeAndSingleCharge(tmp, this->decon_ms1_mz_tol_, false, 1, 10, true, 3, 10, false, true);
      if (tmp.getIntegerDataArrays().empty())
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IntegerDataArray must not be empty!");
      }
      const auto& ida = tmp.getIntegerDataArrays().back();
      if (ida.getName() != "charge")
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IntegerDataArray.back().getName() != \"charge\"");
      }
      // add charges to output map
      for (const auto& q : ida)
      {
        ++cd_[q];
      }
    };
    return f;
  }

  SwathQC::ChargeDistribution SwathQC::getChargeDistribution(const std::vector<SwathMap>& swath_maps, const size_t nr_samples, const double mz_tol)
  {
    ChargeDistribution cd;
    
    SwathQC sq(nr_samples, mz_tol);
    sq.setNrMS1Spectra(0); // leave at 0, such that all incoming spectra are sampled
    auto f_spec = sq.getSpectraProcessingFunc();

    for (const SwathMap& m : swath_maps)
    {
      // only look at MS1 swath maps
      if (!m.ms1) continue;
      
      MSSpectrum s;
      size_t nr_spec = m.sptr->getNrSpectra();
      for (size_t i = 0; i < nr_spec; ++i)
      {
        // we do not convert all spectra from SWATHMap (hence not using the sampling build into getSpectraProcessingFunc())
        // , since this is potentially expensive, but rather only take the ones we need
        if (!isSubsampledSpectrum_(nr_spec, nr_samples, i)) continue;

        OpenMS::OpenSwathDataAccessHelper::convertToOpenMSSpectrum(m.sptr->getSpectrumById(int(i)), s);
        f_spec(s);
      }
    }

    return sq.getChargeDistribution();
  }


  void SwathQC::storeJSON(const OpenMS::String& filename)
  {
    using json = nlohmann::json;

    json out;
    out["ChargeDistributionMS1"] = cd_;
    
    std::ofstream o(filename);
    o << std::setw(2) << out << std::endl;
    // check after writing, to include check for full disk
    if (!o) // fail || bad 
    { 
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    o.close();
  }

  const SwathQC::ChargeDistribution& SwathQC::getChargeDistribution() const
  {
    return cd_;
  }

  void SwathQC::setNrMS1Spectra(size_t nr)
  {
    nr_ms1_spectra_ = nr;
  }

  bool SwathQC::isSubsampledSpectrum_(const size_t total_spec_count, const size_t subsample_count, const size_t idx)
  {
    // if number of MS1 spectra is unknown, we sample everything
    if (total_spec_count == 0) return true;

    if (idx >= total_spec_count) return false;
    if (subsample_count == 0) return false;
    
    // use floating points step-size to ensure uniform sampling from spectra range
    double spec_count = (double)total_spec_count;
    double step_size = spec_count / std::min(spec_count, (double)subsample_count); // guaranteed >= 1

    // estimate the number of steps we need to get to 'idx'
    double steps = idx / step_size;
    // but the number of steps can only be integral ... try both possibilities
    double steps_low = std::floor(steps) * step_size;
    double steps_high = std::ceil(steps) * step_size;

    return (std::lround(steps_low) == (long)idx || std::lround(steps_high) == (long)idx);
  }

}
