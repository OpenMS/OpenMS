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
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <nlohmann/json.hpp>

#include <cmath>       // round

namespace OpenSwath
{
  using namespace OpenMS;

  SwathQC::ChargeDistribution SwathQC::getChargeDistribution(const std::vector<SwathMap>& swath_maps, const int level, const size_t nr_samples, const double mz_tol)
  {
    bool ms1;
    switch (level)
    {
      case 1:
        ms1 = true;
        break;
      case 2:
        ms1 = false;
        break;
      default:
        throw OpenMS::Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Level must be 1 or 2");
    }

    ChargeDistribution cd;
    PeakPickerHiRes pp;

    for (const SwathMap& m : swath_maps)
    {
      // only look at swath maps with the desired level
      if (m.ms1 != ms1) continue;
        
      // use float to ensure uniform sampling from spectra range
      double spec_count = (double)m.sptr->getNrSpectra();
      double step_size = spec_count / std::min(spec_count, (double)nr_samples); // guaranteed >= 1
      MSSpectrum tmp;
      for (double steps = 0; steps < spec_count - 0.5; steps+=step_size) // use -0.5 to avoid summation errors
      {
        MSSpectrum s;
        int id = (int)round(steps);
        OpenMS::OpenSwathDataAccessHelper::convertToOpenMSSpectrum(m.sptr->getSpectrumById(id), s);
        auto t = s.getType(true);
        MSSpectrum* spec;
        if (t == MSSpectrum::SpectrumSettings::PROFILE) 
        {
          pp.pick(s, tmp);
          spec = &tmp;
        }
        else if (t == MSSpectrum::SpectrumSettings::CENTROID)
        {
          spec = &s;
        }
        else
        {
          continue; // unknown: too dangerous to analyse
        }

        // Note: this will pick up also non-peptide signals; filtering by averagine might yield better results
        Deisotoper::deisotopeAndSingleCharge(*spec, mz_tol, false, 1, 10, true, 3, 10, false, true);
        if (spec->getIntegerDataArrays().empty())
        {
          throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IntegerDataArray must not be empty!");
        }
        const auto& ida = spec->getIntegerDataArrays().back();
        if (ida.getName() != "charge")
        {
          throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IntegerDataArray.back().getName() != \"charge\"");
        }
        // add charges to output map
        for (const auto& q : ida)
        {
          ++cd[q];
        }
      }
    }

    return cd;
  }


  void SwathQC::storeJSON(const OpenMS::String& filename, const ChargeDistribution& cd)
  {
    using json = nlohmann::json;

    json out;
    out["ChargeDistributionMS1"] = cd;
    
    std::ofstream o(filename.c_str());
    o << std::setw(2) << out << std::endl;
    o.close();
  }
}