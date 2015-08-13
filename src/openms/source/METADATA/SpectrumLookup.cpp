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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumLookup.h>

using namespace std;

namespace OpenMS
{
  SpectrumLookup::SpectrumLookup(): spectra_(0), n_spectra_(0)
  {
  }


  SpectrumLookup::~SpectrumLookup()
  {}
  

  void SpectrumLookup::setSpectra(vector<MSSpectrum<> >& spectra,
                                  const String& id_regexp_match,
                                  const String& id_regexp_replace)
  {
    spectra_ = &spectra;
    n_spectra_ = spectra.size();
    if (id_regexp_match.empty())
    {
      for (Size i = 0; i < n_spectra_; ++i)
      {
        MSSpectrum<>& spec = spectra[i];
        rts_[spec.getRT()] = i;
        ids_[spec.getNativeID()] = i;
      }
    }
    else
    {
      boost::regex re(id_regexp_match);
      for (Size i = 0; i < n_spectra_; ++i)
      {
        MSSpectrum<>& spec = spectra[i];
        rts_[spec.getRT()] = i;
        String derived_id = boost::regex_replace(spec.getNativeID(), re,
                                                 id_regexp_replace,
                                                 boost::format_no_copy);
        ids_[derived_id] = i;
      }
    }
  }

  
  MSSpectrum<>& SpectrumLookup::findByRT(double rt, double tolerance) const
  {
    double upper_diff = numeric_limits<double>::infinity();
    map<double, Size>::const_iterator upper = rts_.upper_bound(rt);
    if (upper != rts_.end())
    {
      upper_diff = upper->first - rt;
    }
    double lower_diff = numeric_limits<double>::infinity();
    map<double, Size>::const_iterator lower = upper;
    if (lower != rts_.begin())
    {
      --lower;
      lower_diff = rt - lower->first;
    }
    if ((lower_diff < upper_diff) && (lower_diff <= tolerance))
    {
      return (*spectra_)[lower->second];
    }
    if (upper_diff <= tolerance) return (*spectra_)[upper->second];

    String element = "spectrum with RT " + String(rt);
    throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                     element);
  }


  MSSpectrum<>& SpectrumLookup::findByNativeID(const String& native_id) const
  {
    map<String, Size>::const_iterator pos = ids_.find(native_id);
    if (pos == ids_.end())
    {
      String element = "spectrum with native ID '" + native_id + "'";
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       element);
    }
    return (*spectra_)[pos->second];
  }


  MSSpectrum<>& SpectrumLookup::findByIndex(Size index, bool count_from_one)
    const
  {
    Size adjusted_index = index;
    if (count_from_one) --adjusted_index; // overflow (index = 0) handled below
    if (adjusted_index >= n_spectra_)
    {
      String element = "spectrum with index " + String(index);
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       element);
    }
    return (*spectra_)[adjusted_index];
  }

  
  void SpectrumLookup::getSpectrumMetaData(const MSSpectrum<>& spectrum,
                                           SpectrumMetaData& metadata,
                                           MetaDataFlags flags)
  {
    if ((flags & METADATA_RT) == METADATA_RT)
    {
      metadata.rt = spectrum.getRT();
    }
    if ((flags & METADATA_MZ) == METADATA_MZ)
    {
      if (spectrum.getPrecursors().empty()) metadata.mz = 0.0;
      else metadata.mz = spectrum.getPrecursors()[0].getMZ();
    }
    if ((flags & METADATA_CHARGE) == METADATA_CHARGE)
    {
      if (spectrum.getPrecursors().empty()) metadata.charge = 0;
      else metadata.charge = spectrum.getPrecursors()[0].getCharge();
    }
    if ((flags & METADATA_NATIVEID) == METADATA_NATIVEID)
    {
      metadata.native_ID = spectrum.getNativeID();
    }
  }


  void SpectrumLookup::addReferenceFormat(const String& regexp,
                                          bool count_from_one,
                                          double rt_tolerance)
  {
    ReferenceFormat format;
    format.re.assign(regexp);
    format.count_from_one = count_from_one;
    format.rt_tolerance = rt_tolerance;
    reference_formats_.push_back(format);
  }


  MSSpectrum<>& SpectrumLookup::findByRegExpMatch_(const String& spectrum_ref,
                                                   const String& regexp,
                                                   const boost::smatch& match,
                                                   bool count_from_one,
                                                   double rt_tolerance) const
  {
    if (match["SCAN"].matched)
    {
      String value = match["SCAN"].str();
      if (!value.empty()) 
      {
        Size index = value.toInt();
        return findByIndex(index, count_from_one);
      }
    }
    if (match["ID"].matched)
    {
      String value = match["ID"].str();
      if (!value.empty()) 
      {
        return findByNativeID(value);
      }
    }
    if (match["RT"].matched)
    {
      String value = match["RT"].str();
      if (!value.empty())
      {
        double rt = value.toDouble();
        return findByRT(rt, rt_tolerance);
      }
    }
    String msg = "Unexpected format of spectrum reference '" + spectrum_ref +
      "'. The regular expression '" + regexp + "' matched, but no usable "
      "information could be extracted.";
    throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                        msg);
  }


  MSSpectrum<>& SpectrumLookup::findByReference(const String& spectrum_ref)
    const
  {
    for (vector<ReferenceFormat>::const_iterator it =
           reference_formats_.begin(); it != reference_formats_.end(); ++it)
    {
      boost::smatch match;
      bool found = boost::regex_search(spectrum_ref, match, it->re);
      if (found)
      {
        return findByRegExpMatch_(spectrum_ref, it->re.str(), match,
                                  it->count_from_one, it->rt_tolerance);
      }
    }
    String msg = "Spectrum reference doesn't match any known format";
    throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                spectrum_ref, msg);
  }


  void SpectrumLookup::getSpectrumMetaDataByReference(
    const String& spectrum_ref, SpectrumMetaData& metadata, MetaDataFlags flags)
    const
  {
    for (vector<ReferenceFormat>::const_iterator it = 
           reference_formats_.begin(); it != reference_formats_.end(); ++it)
    {
      boost::smatch match;
      bool found = boost::regex_search(spectrum_ref, match, it->re);
      if (found)
      {
        // first try to extract the requested meta data from the reference:
        if (((flags & METADATA_RT) == METADATA_RT) && match["RT"].matched)
        {
          String value = match["RT"].str();
          if (!value.empty())
          {
            metadata.rt = value.toDouble();
            flags &= ~METADATA_RT; // unset flag
          }
        }
        if (((flags & METADATA_MZ) == METADATA_MZ) && match["MZ"].matched)
        {
          String value = match["MZ"].str();
          if (!value.empty())
          {
            metadata.mz = value.toDouble();
            flags &= ~METADATA_MZ; // unset flag
          }
        }
        if (((flags & METADATA_CHARGE) == METADATA_CHARGE) &&
            match["CHARGE"].matched)
        {
          String value = match["CHARGE"].str();
          if (!value.empty())
          {
            metadata.charge = value.toDouble();
            flags &= ~METADATA_CHARGE; // unset flag
          }
        }
        if (((flags & METADATA_NATIVEID) == METADATA_NATIVEID) &&
            match["ID"].matched)
        {
          metadata.native_ID = match["ID"].str();
          if (!metadata.native_ID.empty())
          {
            flags &= ~METADATA_NATIVEID; // unset flag
          }
        }
        if (flags) // not all requested values have been found -> look them up
        {
          MSSpectrum<>& spec = findByRegExpMatch_(spectrum_ref, it->re.str(),
                                                  match, it->count_from_one,
                                                  it->rt_tolerance);
          getSpectrumMetaData(spec, metadata, flags);
        }
      }
    }
  }

} // namespace OpenMS
