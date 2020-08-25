// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
  const String& SpectrumLookup::default_scan_regexp = "=(?<SCAN>\\d+)$";

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
    if (count_from_one) --adjusted_index; // overflow (index = 0) handled below
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
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                spectrum_ref, msg);
  }


  Int SpectrumLookup::extractScanNumber(const String& native_id,
                                        const boost::regex& scan_regexp, 
                                        bool no_error)
  {
    vector<string> matches;
    boost::sregex_token_iterator current_begin(native_id.begin(), native_id.end(), scan_regexp, 1);
    boost::sregex_token_iterator current_end(native_id.end(), native_id.end(), scan_regexp, 1);
    matches.insert(matches.end(), current_begin, current_end);
    if (!matches.empty())
    {
      // always use the last possible matching subgroup 
      String last_value = String(matches.back());
      try
      {
        return last_value.toInt();
      }
      catch (Exception::ConversionError&)
      {
      }
    }
    if (!no_error)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  native_id, "Could not extract scan number");
    }
    return -1;
  }

  Int SpectrumLookup::extractScanNumber(const String& native_id,
                                        const String& native_id_type_accession)
  {
    // check accession for data type to extract (e.g. MS:1000768 - Thermo nativeID format - scan=xsd:positiveInteger)
    boost::regex regexp;
    // list of CV accessions with native id format "scan=NUMBER"
    std::vector<String> scan = {"MS:1000768","MS:1000769","MS:1000771","MS:1000772","MS:1000776"};
    // list of CV accession with native id format "file=NUMBER"
    std::vector<String> file = {"MS:1000773","MS:1000775"};
    // expected number of subgroups
    vector<int> subgroups = {1};
    
    // "scan=NUMBER" 
    if (std::find(scan.begin(), scan.end(), native_id_type_accession) != scan.end())
    {
      regexp = std::string("scan=(?<GROUP>\\d+)");
    }
    // id="sample=1 period=1 cycle=96 experiment=1" - this will be described by a combination of (cycle * 1000 + experiment)
    else if (native_id_type_accession == "MS:1000770") // WIFF nativeID format
    {
      regexp = std::string("cycle=(?<GROUP>\\d+).experiment=(?<GROUP>\\d+)");
      subgroups = {1, 2};
    }
    // "file=NUMBER"
    else if (std::find(file.begin(), file.end(), native_id_type_accession) != file.end())
    {
      regexp = std::string("file=(?<GROUP>\\d+)");
    }
    // "index=NUMBER"
    else if (native_id_type_accession == "MS:1000774")
    {
      regexp = std::string("index=(?<GROUP>\\d+)");
    }
    // "spectrum=NUMBER"
    else if (native_id_type_accession == "MS:1000777")
    {
      regexp = std::string("spectrum=(?<GROUP>\\d+)");
    }
    // NUMBER 
    else if (native_id_type_accession == "MS:1001530")  
    {
      regexp = std::string("(?<GROUP>\\d+)");
    }
    else
    {
      OPENMS_LOG_WARN << "native_id: " << native_id << " accession: " << native_id_type_accession << " Could not extract scan number - no valid native_id_type_accession was provided" << std::endl;
    }

    if (!regexp.empty()) 
    {
      vector<string> matches;
      boost::sregex_token_iterator current_begin(native_id.begin(), native_id.end(), regexp, subgroups);
      boost::sregex_token_iterator current_end(native_id.end(), native_id.end(), regexp, subgroups);
      matches.insert(matches.end(), current_begin, current_end);
      if (matches.size() == 1) // default case: one native identifier
      {
        try
        {
          String value = String(matches[0]);
          return value.toInt();
        }
        catch (Exception::ConversionError&)
        {
          OPENMS_LOG_WARN << "Value: '" << String(matches[0]) << "' could not be converted to int in string. Native ID='" << native_id << "'" << std::endl;
          return -1;
        }
      }
      else if (matches.size() == 2) // special case: wiff file with two native identifiers
      {
        try
        {
          if (String(matches[1]).toInt() < 1000) // checks if value of experiment is smaller than 1000 (cycle * 1000 + experiment)
          {
            int value = String(matches[0]).toInt() * 1000 + String(matches[1]).toInt();
            return value; 
          }
          else
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The value of experiment is too large and can not be handled properly.", String(matches[1]));
          }
        }
        catch (Exception::ConversionError&)
        {
          OPENMS_LOG_WARN << "Value: '" << String(matches[0]) << "' could not be converted to int in string. Native ID='" 
            << native_id << "' accession='" << native_id_type_accession << "'" << std::endl;
          return -1;
        }
      }
      else
      {
        return -1;
      }
    }
    return -1;
  } 

  void SpectrumLookup::addEntry_(Size index, double rt, Int scan_number,
                                 const String& native_id)
  {
    rts_[rt] = index;
    ids_[native_id] = index;
    if (scan_number != -1) scans_[scan_number] = index;
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

