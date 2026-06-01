// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <vector>

using namespace std;

namespace OpenMS
{
  bool SpectrumNativeIDParser::isNativeID(const String& id)
  {
    return id.hasPrefix("scan=") || id.hasPrefix("scanId=") || id.hasPrefix("scanID=")
        || id.hasPrefix("controllerType=") || id.hasPrefix("function=") || id.hasPrefix("sample=")
        || id.hasPrefix("index=") || id.hasPrefix("spectrum=") || id.hasPrefix("file=")
        || id.hasPrefix("frame=");
  }

  std::string SpectrumNativeIDParser::getRegExFromNativeID(const String& id)
  {
    // "scan=NUMBER" e.g. Bruker/Agilent
    // "controllerType=0 controllerNumber=1 scan=NUMBER" for Thermo
    // "function= process= scan=NUMBER" for Waters
    // "frame=FRAME_ID scan=SCAN_ID [precursor=PREC_ID]" for Bruker TDF
    //   (MS:1002818). Extra trailing tokens (precursor=, windowGroup=,
    //   scanStart=, scanEnd=, merged=) are used by OpenMS DDA/DIA PASEF
    //   and pwiz combined mode; the regex targets the "scan=<int>" token
    //   which remains the most meaningful scan-number proxy.
    if (id.hasPrefix("scan=")
     || id.hasPrefix("controllerType=")
     || id.hasPrefix("function=")
     || id.hasPrefix("frame=")) return std::string(R"(scan=(?<GROUP>\d+))");

    // "index=NUMBER"
    if (id.hasPrefix("index=")) return std::string(R"(index=(?<GROUP>\d+))");

    // "scanId=NUMBER" or "scanID=NUMBER" - MS_Agilent_MassHunter_nativeID_format
    if (id.hasPrefix("scanId=")) return std::string(R"(scanId=(?<GROUP>\d+))");
    if (id.hasPrefix("scanID=")) return std::string(R"(scanID=(?<GROUP>\d+))");

    // "spectrum=NUMBER"
    if (id.hasPrefix("spectrum=")) return std::string(R"(spectrum=(?<GROUP>\d+))");

    // "file=NUMBER" Bruker FID or single peak list
    if (id.hasPrefix("file=")) return std::string(R"(file=(?<GROUP>\d+))");

    // NUMBER
    return std::string(R"((?<GROUP>\d+))");
  }

  Int SpectrumNativeIDParser::extractScanNumber(const String& native_id,
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

  Int SpectrumNativeIDParser::extractScanNumber(const String& native_id,
                                        const String& native_id_type_accession)
  {
    // check accession for data type to extract (e.g. MS:1000768 - Thermo nativeID format - scan=xsd:positiveInteger)
    boost::regex regexp;
    // list of CV accessions with native id format "scan=NUMBER".
    // MS:1002818 = "Bruker TDF nativeID format" (pattern "frame=<int> scan=<int>",
    // we extend it with a trailing "precursor=<int>" for aggregated DDA
    // output — see BrukerTimsFile::loadDDA_ for the rationale). The regex
    // targets the "scan=<int>" token in all cases; trailing tokens are
    // ignored.
    std::vector<String> scan = {"MS:1000768","MS:1000769","MS:1000771","MS:1000772","MS:1000776","MS:1002818"};
    // list of CV accession with native id format "file=NUMBER"
    std::vector<String> file = {"MS:1000773","MS:1000775"};
    // expected number of subgroups
    vector<int> subgroups = {1};

    // "scan=NUMBER"
    if (std::find(scan.begin(), scan.end(), native_id_type_accession) != scan.end())
    {
      regexp = std::string(R"(scan=(?<GROUP>\d+))");
    }
    // id="sample=1 period=1 cycle=96 experiment=1" - this will be described by a combination of (cycle * 1000 + experiment)
    else if (native_id_type_accession == "MS:1000770") // WIFF nativeID format
    {
      regexp = std::string(R"(cycle=(?<GROUP>\d+)\s+experiment=(?<GROUP>\d+))");
      subgroups = {1, 2};
    }
    // "file=NUMBER"
    else if (std::find(file.begin(), file.end(), native_id_type_accession) != file.end())
    {
      regexp = std::string(R"(file=(?<GROUP>\d+))");
    }
    // "index=NUMBER"
    else if (native_id_type_accession == "MS:1000774")
    {
      regexp = std::string(R"(index=(?<GROUP>\d+))");
    }
    // "scanId=NUMBER" - MS_Agilent_MassHunter_nativeID_format
    else if (native_id_type_accession == "MS:1001508")
    {
      regexp = std::string(R"(scanId=(?<GROUP>\d+))");
    }
    // "spectrum=NUMBER"
    else if (native_id_type_accession == "MS:1000777")
    {
      regexp = std::string(R"(spectrum=(?<GROUP>\d+))");
    }
    // NUMBER
    else if (native_id_type_accession == "MS:1001530")
    {
      regexp = std::string(R"((?<GROUP>\d+))");
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

      if (matches.size() < subgroups.size()) {
          OPENMS_LOG_WARN << "native_id '" << native_id <<"' is invalid. Could not extract scan number." << std::endl;
          return -1;
      }

      if (subgroups.size() == 1) // default case: one native identifier
      {
        try
        {
          // In case of merged spectra the last native id matches the scan number of the merged scan.
          String value = String(matches[matches.size() - 1]);
          if (native_id_type_accession == "MS:1000774")
          {
            return value.toInt() + 1; // if the native ID is index=.., the scan number is usually considered index+1 (especially for pepXML)
          }
          else
          {
            return value.toInt();
          }
        }
        catch (Exception::ConversionError&)
        {
          OPENMS_LOG_WARN << "Value: '" << String(matches[matches.size() - 1]) << "' could not be converted to int in string. Native ID='" << native_id << "'" << std::endl;
          return -1;
        }
      }
      else if (subgroups.size() == 2) // special case: wiff file with two native identifiers
      {
        try
        {
          // In case of merged spectra the last native id matches the scan number of the merged scan.
          String cycle_str = matches[matches.size() - 2];
          String experiment_str = matches[matches.size() - 1];

          if (experiment_str.toInt() < 1000) // checks if value of experiment is smaller than 1000 (cycle * 1000 + experiment)
          {
            int value = cycle_str.toInt() * 1000 + experiment_str.toInt();
            return value;
          }
          else
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The value of experiment is too large and can not be handled properly.", experiment_str);
          }
        }
        catch (Exception::ConversionError&)
        {
          OPENMS_LOG_WARN << "Values: '" << matches[matches.size() - 2] << "', '" << matches[matches.size() - 1] << "' could not be converted to int in string. Native ID='"
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

} // namespace OpenMS
