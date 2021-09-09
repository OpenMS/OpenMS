// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/IONMOBILITY/IMTypes.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  const std::string NamesOfDriftTimeUnit[] = {"<NONE>", "ms", "1/K0", "FAIMS_CV"};
  const std::string NamesOfIMFormat[] = {"none", "concatenated", "multiple_spectra", "mixed"};


 DriftTimeUnit toDriftTimeUnit(const std::string& dtu_string)
  {
    auto first = &NamesOfDriftTimeUnit[0];
    auto last = &NamesOfDriftTimeUnit[(size_t) DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT];
    const auto it = std::find(first, last, dtu_string);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", dtu_string);
    }
    return DriftTimeUnit(it - first);
  }

  const std::string& toString(const DriftTimeUnit value)
  {
    if (value == DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_DRIFTTIMEUNIT");
    }
    return NamesOfDriftTimeUnit[(size_t) value];
  }

  IMFormat toIMFormat(const std::string& IM_format)
  {
    auto first = &NamesOfIMFormat[0];
    auto last = &NamesOfIMFormat[(size_t) IMFormat::SIZE_OF_IMFORMAT];
    const auto it = std::find(first, last, IM_format);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", IM_format);
    }
    return IMFormat(it - first);
  }

  const std::string& toString(const IMFormat value)
  {
    if (value == IMFormat::SIZE_OF_IMFORMAT)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_IMFORMAT");
    }
    return NamesOfIMFormat[(size_t)value];
  }

  IMFormat IMTypes::determineIMFormat(const MSExperiment& exp)
  {
    std::set<IMFormat> occs;
    for (const auto& spec : exp.getSpectra())
    {
      occs.insert(determineIMFormat(spec));
    }
    occs.erase(IMFormat::NONE); // ignore NONE (i.e. normal spectra)

    if (occs.empty())
    {
      return IMFormat::NONE;
    }
    if (occs.size() == 1 && (occs.find(IMFormat::CONCATENATED) != occs.end() || occs.find(IMFormat::MULTIPLE_SPECTRA) != occs.end()))
    {
      return *occs.begin();
    }
    if (occs.size() == 2 && occs.find(IMFormat::CONCATENATED) != occs.end() && occs.find(IMFormat::MULTIPLE_SPECTRA) != occs.end())
    {
      return IMFormat::MIXED;
    }

    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "subfunction returned invalid value(s)", "Number of different values: " + String(occs.size()));
  }

  IMFormat IMTypes::determineIMFormat(const MSSpectrum& spec)
  {
    bool has_float_data = spec.containsIMData(); // cache value; query is 'expensive'
    bool has_drift_time = spec.getDriftTime() != DRIFTTIME_NOT_SET;
    if (has_float_data && has_drift_time)
    {
      const auto& fda = spec.getFloatDataArrays()[spec.getIMData().first];
      String array_val = fda.empty() ? "[empty]" : String(fda[0]);
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSSpectrum contains both an float-data-array and a single drift time. At most one is allowed per spectrum!", String("Array: ") + array_val + ", ... <> Spec: " + spec.getDriftTime());
    }

    if (has_float_data)
    {
      return IMFormat::CONCATENATED;
    }
    else if (has_drift_time)
    {
      if (spec.getDriftTimeUnit() == DriftTimeUnit::NONE)
      {
        OPENMS_LOG_WARN << "Warning: no drift time unit set for spectrum " << spec.getNativeID() << "\n";
      }
      return IMFormat::MULTIPLE_SPECTRA;
    }
    return IMFormat::NONE;
  }



}// namespace OpenMS
