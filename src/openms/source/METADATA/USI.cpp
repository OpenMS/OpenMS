// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/USI.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <boost/regex.hpp>
#include <sstream>

namespace OpenMS
{
  // Static member definitions
  const std::string USI::USI_PREFIX = "mzspec:";
  const std::string USI::CV_ACCESSION = "MS:1003063";
  const std::string USI::CV_NAME = "universal spectrum identifier";

  USI::USI() :
    collection_(),
    ms_run_(),
    index_type_(IndexType::SCAN),
    index_(),
    interpretation_()
  {
  }

  USI::USI(const std::string& collection,
           const std::string& ms_run,
           IndexType index_type,
           const std::string& index,
           const std::string& interpretation) :
    collection_(collection),
    ms_run_(ms_run),
    index_type_(index_type),
    index_(index),
    interpretation_(interpretation)
  {
  }

  USI::USI(const std::string& usi_string) :
    collection_(),
    ms_run_(),
    index_type_(IndexType::SCAN),
    index_(),
    interpretation_()
  {
    if (!fromString(usi_string))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  usi_string, "Invalid USI format");
    }
  }

  bool USI::operator==(const USI& rhs) const
  {
    return collection_ == rhs.collection_ &&
           ms_run_ == rhs.ms_run_ &&
           index_type_ == rhs.index_type_ &&
           index_ == rhs.index_ &&
           interpretation_ == rhs.interpretation_;
  }

  bool USI::operator!=(const USI& rhs) const
  {
    return !(*this == rhs);
  }

  bool USI::operator<(const USI& rhs) const
  {
    if (collection_ != rhs.collection_)
      return collection_ < rhs.collection_;
    if (ms_run_ != rhs.ms_run_)
      return ms_run_ < rhs.ms_run_;
    if (index_type_ != rhs.index_type_)
      return index_type_ < rhs.index_type_;
    if (index_ != rhs.index_)
      return index_ < rhs.index_;
    return interpretation_ < rhs.interpretation_;
  }

  bool USI::isValid() const
  {
    return !collection_.empty() && !ms_run_.empty() && !index_.empty();
  }

  bool USI::isValidUSI(const std::string& usi_string)
  {
    return tryParse(usi_string).has_value();
  }

  std::optional<USI> USI::tryParse(const std::string& usi_string)
  {
    // Quick prefix check
    if (!StringUtils::hasPrefix(usi_string, "mzspec:"))
    {
      return std::nullopt;
    }

    USI usi;
    if (usi.fromString(usi_string))
    {
      return usi;
    }
    return std::nullopt;
  }

  const std::string& USI::getCollection() const
  {
    return collection_;
  }

  void USI::setCollection(const std::string& collection)
  {
    collection_ = collection;
  }

  const std::string& USI::getMSRun() const
  {
    return ms_run_;
  }

  void USI::setMSRun(const std::string& ms_run)
  {
    ms_run_ = ms_run;
  }

  USI::IndexType USI::getIndexType() const
  {
    return index_type_;
  }

  void USI::setIndexType(IndexType index_type)
  {
    index_type_ = index_type;
  }

  const std::string& USI::getIndex() const
  {
    return index_;
  }

  void USI::setIndex(const std::string& index)
  {
    index_ = index;
  }

  const std::string& USI::getInterpretation() const
  {
    return interpretation_;
  }

  void USI::setInterpretation(const std::string& interpretation)
  {
    interpretation_ = interpretation;
  }

  bool USI::hasInterpretation() const
  {
    return !interpretation_.empty();
  }

  std::string USI::toString() const
  {
    if (!isValid())
    {
      return "";
    }

    std::ostringstream oss;
    oss << USI_PREFIX << collection_ << ":" << ms_run_ << ":"
        << indexTypeToString(index_type_) << ":" << index_;

    if (!interpretation_.empty())
    {
      oss << ":" << interpretation_;
    }

    return oss.str();
  }

  bool USI::fromString(const std::string& usi_string)
  {
    // Reset current values
    collection_.clear();
    ms_run_.clear();
    index_type_ = IndexType::SCAN;
    index_.clear();
    interpretation_.clear();

    // Check prefix
    if (!StringUtils::hasPrefix(usi_string, USI_PREFIX))
    {
      return false;
    }

    // Remove prefix and split by colon
    std::string remainder = StringUtils::substr(usi_string, USI_PREFIX.size());
    std::vector<std::string> parts;
    StringUtils::split(remainder, ':', parts);

    // Minimum required: collection, ms_run, index_type, index (4 parts)
    if (parts.size() < 4)
    {
      return false;
    }

    collection_ = parts[0];
    ms_run_ = parts[1];

    // Parse index type
    try
    {
      index_type_ = indexTypeFromString(parts[2]);
    }
    catch (const Exception::InvalidValue&)
    {
      return false;
    }

    index_ = parts[3];

    // Optional interpretation (remaining parts joined by colon)
    if (parts.size() > 4)
    {
      std::ostringstream oss;
      for (size_t i = 4; i < parts.size(); ++i)
      {
        if (i > 4)
        {
          oss << ":";
        }
        oss << parts[i];
      }
      interpretation_ = oss.str();
    }

    return isValid();
  }

  USI USI::createFromScanNumber(const std::string& dataset_id,
                                const std::string& filename,
                                int scan_number,
                                const std::string& interpretation)
  {
    return USI(dataset_id, filename, IndexType::SCAN,StringUtils::toStr(scan_number), interpretation);
  }

  USI USI::createFromNativeID(const std::string& dataset_id,
                              const std::string& filename,
                              const std::string& native_id)
  {
    // Try to extract scan number from native ID
    auto scan_num = extractScanNumberFromNativeID(native_id);
    
    if (scan_num.has_value())
    {
      return USI(dataset_id, filename, IndexType::SCAN,StringUtils::toStr(scan_num.value()));
    }
    else
    {
      // Fall back to using nativeId index type
      return USI(dataset_id, filename, IndexType::NATIVEID, native_id);
    }
  }

  std::optional<int> USI::extractScanNumberFromNativeID(const std::string& native_id)
  {
    // Use SpectrumNativeIDParser for comprehensive native ID parsing
    // It handles various vendor formats: Thermo, Waters, Bruker, Agilent, etc.
    if (SpectrumNativeIDParser::isNativeID(native_id))
    {
      std::string regex_str = SpectrumNativeIDParser::getRegExFromNativeID(native_id);
      if (!regex_str.empty())
      {
        boost::regex scan_regexp(regex_str);
        Int scan = SpectrumNativeIDParser::extractScanNumber(native_id, scan_regexp, true);
        if (scan >= 0)
        {
          return scan;
        }
      }
    }

    // Try to parse the entire string as a number (simple case like "12345")
    try
    {
      return StringUtils::toInt32(native_id);
    }
    catch (const Exception::ConversionError&)
    {
      return std::nullopt;
    }
  }

  std::string USI::indexTypeToString(IndexType index_type)
  {
    switch (index_type)
    {
      case IndexType::SCAN:
        return "scan";
      case IndexType::INDEX:
        return "index";
      case IndexType::NATIVEID:
        return "nativeId";
      default:
        return "scan";
    }
  }

  USI::IndexType USI::indexTypeFromString(const std::string& type_string)
  {
    std::string lower = type_string;
    StringUtils::toLower(lower);

    if (lower == "scan")
    {
      return IndexType::SCAN;
    }
    else if (lower == "index")
    {
      return IndexType::INDEX;
    }
    else if (lower == "nativeid")
    {
      return IndexType::NATIVEID;
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Invalid USI index type", type_string);
    }
  }

  std::string USI::extractBasename(const std::string& filepath)
  {
    if (filepath.empty())
    {
      return "";
    }

    std::string path = filepath;

    // Handle file:// URIs
    if (StringUtils::hasPrefix(path, "file://"))
    {
      path = StringUtils::substr(path, 7);  // Remove "file://"
      // Handle Windows paths like file:///C:/path
      if (path.size() > 2 && path[0] == '/' && path[2] == ':')
      {
        path = StringUtils::substr(path, 1);  // Remove leading slash for Windows paths
      }
    }

    // Find the last path separator (handle both Unix and Windows)
    size_t last_sep = path.rfind('/');
    size_t last_win_sep = path.rfind('\\');
    
    if (last_win_sep != std::string::npos && (last_sep == std::string::npos || last_win_sep > last_sep))
    {
      last_sep = last_win_sep;
    }

    if (last_sep != std::string::npos)
    {
      return StringUtils::substr(path, last_sep + 1);
    }

    return path;
  }

  const std::string& USI::getCVAccession()
  {
    return CV_ACCESSION;
  }

  const std::string& USI::getCVName()
  {
    return CV_NAME;
  }

  std::ostream& operator<<(std::ostream& os, const USI& usi)
  {
    os << usi.toString();
    return os;
  }

} // namespace OpenMS
