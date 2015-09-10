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

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace std;

namespace OpenMS
{
  const String& SpectrumLookup::regexp_names_ = "INDEX0 INDEX1 SCAN ID RT";

  SpectrumLookup::SpectrumLookup(): 
    rt_tolerance(0.01), n_spectra_(0),
    regexp_name_list_(ListUtils::create<String>(regexp_names_, ' '))
  {}


  SpectrumLookup::~SpectrumLookup()
  {}


  bool SpectrumLookup::empty() const
  {
    return n_spectra_ == 0;
  }
  

  template <typename SpectrumContainer>
  void SpectrumLookup::readSpectra(const SpectrumContainer& spectra, 
                                   const String& scan_regexp)
  {
    n_spectra_ = spectra.size();
    rts_.clear();
    ids_.clear();
    scans_.clear();
    boost::regex re;
    if (!scan_regexp.empty())
    {
      if (!scan_regexp.hasSubstring("?<SCAN>"))
      {
        String msg = "The regular expression for extracting scan numbers from native IDs must contain a named group '?<SCAN>'.";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         __PRETTY_FUNCTION__, msg);
      }
      re.assign(scan_regexp);
    }
    for (Size i = 0; i < n_spectra_; ++i)
    {
      const MSSpectrum<>& spec = spectra[i];
      String native_id = spec.getNativeID();
      rts_[spec.getRT()] = i;
      ids_[native_id] = i;
      if (!scan_regexp.empty())
      {
        boost::smatch match;
        bool found = boost::regex_search(native_id, match, re);
        if (found && match["SCAN"].matched)
        {
          String value = match["SCAN"].str();
          try
          {
            Size scan_no = value.toInt();
            scans_[scan_no] = i;
            continue;
          }
          catch(Exception::ConversionError& e)
          {
          }
        }
        LOG_WARN << "Warning: Could not extract scan number from spectrum native ID '" + native_id + "' using regular expression '" + scan_regexp + "'. Look-up by scan number may not work properly." << endl;
      }
    }
  }

  
  Size SpectrumLookup::findByRT(double rt) const
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
    if ((lower_diff < upper_diff) && (lower_diff <= rt_tolerance))
    {
      return lower->second;
    }
    if (upper_diff <= rt_tolerance) return upper->second;

    String element = "spectrum with RT " + String(rt);
    throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                     element);
  }


  Size SpectrumLookup::findByNativeID(const String& native_id) const
  {
    map<String, Size>::const_iterator pos = ids_.find(native_id);
    if (pos == ids_.end())
    {
      String element = "spectrum with native ID '" + native_id + "'";
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       element);
    }
    return pos->second;
  }


  Size SpectrumLookup::findByIndex(Size index, bool count_from_one) const
  {
    Size adjusted_index = index;
    if (count_from_one) --adjusted_index; // overflow (index = 0) handled below
    if (adjusted_index >= n_spectra_)
    {
      String element = "spectrum with index " + String(index);
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       element);
    }
    return adjusted_index;
  }


  Size SpectrumLookup::findByScanNumber(Size scan_number) const
  {
    map<Size, Size>::const_iterator pos = scans_.find(scan_number);
    if (pos == scans_.end())
    {
      String element = "spectrum with scan number " + String(scan_number);
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       element);
    }
    return pos->second;
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


  void SpectrumLookup::addReferenceFormat(const String& regexp)
  {
    // does the reg. exp. contain any of the recognized group names?
    bool found = false;
    for (vector<String>::iterator it = regexp_name_list_.begin();
         it != regexp_name_list_.end(); ++it)
    {
      if (regexp.hasSubstring("?<" + (*it) + ">"))
      {
        found = true;
        break;
      }
    }
    if (!found)
    {
      String msg = "The regular expression describing the reference format must contain at least one of the following named groups (in the format '?<GROUP>'): " + regexp_names_;
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       msg);
    }

    boost::regex re(regexp);
    reference_formats.push_back(re);
  }


  Size SpectrumLookup::findByRegExpMatch_(const String& spectrum_ref,
                                          const String& regexp, 
                                          const boost::smatch& match) const
  {
    if (match["INDEX0"].matched)
    {
      String value = match["INDEX0"].str();
      if (!value.empty()) 
      {
        Size index = value.toInt();
        return findByIndex(index, false);
      }
    }
    if (match["INDEX1"].matched)
    {
      String value = match["INDEX1"].str();
      if (!value.empty()) 
      {
        Size index = value.toInt();
        return findByIndex(index, true);
      }
    }
    if (match["SCAN"].matched)
    {
      String value = match["SCAN"].str();
      if (!value.empty()) 
      {
        Size scan_number = value.toInt();
        return findByScanNumber(scan_number);
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
        return findByRT(rt);
      }
    }
    String msg = "Unexpected format of spectrum reference '" + spectrum_ref +
      "'. The regular expression '" + regexp + "' matched, but no usable "
      "information could be extracted.";
    throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                        msg);
  }


  Size SpectrumLookup::findByReference(const String& spectrum_ref) const
  {
    for (vector<boost::regex>::const_iterator it = reference_formats.begin();
         it != reference_formats.end(); ++it)
    {
      boost::smatch match;
      bool found = boost::regex_search(spectrum_ref, match, *it);
      if (found)
      {
        return findByRegExpMatch_(spectrum_ref, it->str(), match);
      }
    }
    String msg = "Spectrum reference doesn't match any known format";
    throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                spectrum_ref, msg);
  }


  template <typename SpectrumContainer>
  void SpectrumLookup::getSpectrumMetaDataByReference(
    const SpectrumContainer& spectra, const String& spectrum_ref,
    SpectrumMetaData& metadata, MetaDataFlags flags) const
  {
    for (vector<boost::regex>::const_iterator it = reference_formats.begin();
         it != reference_formats.end(); ++it)
    {
      boost::smatch match;
      bool found = boost::regex_search(spectrum_ref, match, *it);
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
          Size index = findByRegExpMatch_(spectrum_ref, it->str(), match);
          const MSSpectrum<>& spectrum = spectra[index];
          getSpectrumMetaData(spectrum, metadata, flags);
        }
      }
    }
  }


  bool SpectrumLookup::addMissingRTsToPeptideIDs(
    vector<PeptideIdentification>& peptides, const String& filename,
    bool stop_on_error)
  {
    SpectrumLookup lookup;
    MSExperiment<> exp;
    bool success = true;
    for (vector<PeptideIdentification>::iterator it = peptides.begin();
         it != peptides.end(); ++it)
    {
      if (boost::math::isnan(it->getRT()))
      {
        if (lookup.empty()) // load raw data only if we have to
        {
          FileHandler().loadExperiment(filename, exp);
          lookup.readSpectra(exp);
        }
        String spectrum_id = it->getMetaValue("spectrum_reference");
        try
        {
          Size index = lookup.findByNativeID(spectrum_id);
          it->setRT(exp[index].getRT());
        }
        catch(Exception::ElementNotFound& e)
        {
          LOG_ERROR << "Error: Failed to look up retention time for peptide ID with spectrum reference '" + spectrum_id + "' - no spectrum with corresponding native ID found." << endl;
          success = false;
          if (stop_on_error) break;
        }
      }
    }
    return success;
  }

} // namespace OpenMS
