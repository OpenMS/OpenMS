// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <boost/regex.hpp>

using namespace std;

namespace OpenMS
{
  const String& SpectrumLookup::default_scan_regexp = R"(=(?<SCAN>\d+)$)";

  const String& SpectrumLookup::regexp_names_ = "INDEX0 INDEX1 SCAN ID RT";

  SpectrumLookup::SpectrumLookup(): 
    rt_tolerance(0.01), n_spectra_(0),
    regexp_name_list_(ListUtils::create<String>(regexp_names_, ' '))
  {}


  SpectrumLookup::~SpectrumLookup() = default;


  bool SpectrumLookup::empty() const
  {
    return n_spectra_ == 0;
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
    if (upper_diff <= rt_tolerance)
    {
      return upper->second;
    }
    String element = "spectrum with RT " + String(rt);
    throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     element);
  }


  Size SpectrumLookup::findByNativeID(const String& native_id) const
  {
    map<String, Size>::const_iterator pos = ids_.find(native_id);
    if (pos == ids_.end())
    {
      String element = "spectrum with native ID '" + native_id + "'";
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       element);
    }
    return pos->second;
  }


  Size SpectrumLookup::findByIndex(Size index, bool count_from_one) const
  {
    Size adjusted_index = index;
    if (count_from_one)
    {
      --adjusted_index; // overflow (index = 0) handled below
    }
    if (adjusted_index >= n_spectra_)
    {
      String element = "spectrum with index " + String(index);
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
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
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       element);
    }
    return pos->second;
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
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
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
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        msg);
  }

  Size SpectrumLookup::findByReference(const String& spectrum_ref) const
  {
    for (const boost::regex& reg : reference_formats)
    {
      boost::smatch match;
      bool found = boost::regex_search(spectrum_ref, match, reg);
      if (found)
      {
        return findByRegExpMatch_(spectrum_ref, reg.str(), match);
      }
    }
    String msg = "Spectrum reference doesn't match any known format";
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                spectrum_ref, msg);
  }


  bool SpectrumLookup::isNativeID(const String& id)
  {
    return SpectrumNativeIDParser::isNativeID(id);
  }

  std::string SpectrumLookup::getRegExFromNativeID(const String& id)
  {
    return SpectrumNativeIDParser::getRegExFromNativeID(id);
  }

  Int SpectrumLookup::extractScanNumber(const String& native_id,
                                        const boost::regex& scan_regexp,
                                        bool no_error)
  {
    return SpectrumNativeIDParser::extractScanNumber(native_id, scan_regexp, no_error);
  }

  Int SpectrumLookup::extractScanNumber(const String& native_id,
                                        const String& native_id_type_accession)
  {
    return SpectrumNativeIDParser::extractScanNumber(native_id, native_id_type_accession);
  } 

  void SpectrumLookup::addEntry_(Size index, double rt, Int scan_number,
                                 const String& native_id)
  {
    rts_[rt] = index;
    ids_[native_id] = index;
    if (scan_number != -1)
    {
      scans_[scan_number] = index;
    }
  }

  void SpectrumLookup::setScanRegExp_(const String& scan_regexp)
  {
    if (!scan_regexp.empty())
    {
      if (!scan_regexp.hasSubstring("?<SCAN>"))
      {
        String msg = "The regular expression for extracting scan numbers from native IDs must contain a named group '?<SCAN>'.";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      scan_regexp_.assign(scan_regexp);
    }
  }

} // namespace OpenMS

