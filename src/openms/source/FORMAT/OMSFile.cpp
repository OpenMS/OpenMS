// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/FORMAT/OMSFileLoad.h>
#include <OpenMS/FORMAT/OMSFileStore.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/SYSTEM/TempFiles.h>
#include <filesystem>
#ifdef _WIN32
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #include <windows.h>
#endif
#include <fstream>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS
{
namespace
{
  template<class Value>
  void storeAtomic(const std::string& filename, const Value& value, ProgressLogger::LogType log_type)
  {
    const auto utf8 = [](const std::filesystem::path& path) {
      const auto text = path.u8string();
      return std::string(reinterpret_cast<const char*>(text.data()), text.size());
    };
    const auto target = std::filesystem::absolute(std::filesystem::u8path(filename));
    TempDir staging(utf8(target.parent_path()));
    const auto temporary_path = std::filesystem::u8path(staging.getPath()) / "output.oms";
    const auto temporary = utf8(temporary_path);
    {
      Internal::OMSFileStore helper(temporary, log_type);
      helper.store(value);
    }
    // Rename only after validation, all writes and closing the SQLite connection.
    std::error_code error;
#ifdef _WIN32
    if (! MoveFileExW(temporary_path.c_str(), target.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
      error = std::error_code(static_cast<int>(GetLastError()), std::system_category());
#else
    std::filesystem::rename(temporary_path, target, error);
#endif
    if (error) throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename + ": " + error.message());
  }
} // namespace
void OMSFile::store(const std::string& filename, const IdentificationData& id_data)
{ storeAtomic(filename, id_data, log_type_); }

void OMSFile::store(const std::string& filename, const FeatureMap& features)
{
  if (features.getIdentificationData().empty())
  {
    auto converted = features;
    IdentificationDataConverter::importFeatureIDs(converted);
    storeAtomic(filename, converted, log_type_);
  }
  else
    storeAtomic(filename, features, log_type_);
}

  void OMSFile::store(const std::string& filename, const ConsensusMap& consensus)
  {
    if (consensus.getIdentificationData().empty())
    {
      auto converted = consensus;
      IdentificationDataConverter::importConsensusIDs(converted);
      storeAtomic(filename, converted, log_type_);
    }
    else
      storeAtomic(filename, consensus, log_type_);
  }

  void OMSFile::load(const std::string& filename, IdentificationData& id_data)
  {
    OpenMS::Internal::OMSFileLoad helper(filename, log_type_);
    IdentificationData loaded;
    helper.load(loaded);
    id_data = std::move(loaded);
  }

  void OMSFile::load(const std::string& filename, FeatureMap& features)
  {
    OpenMS::Internal::OMSFileLoad helper(filename, log_type_);
    FeatureMap loaded;
    helper.load(loaded);
    features = std::move(loaded);
  }

  void OMSFile::load(const std::string& filename, ConsensusMap& consensus)
  {
    OpenMS::Internal::OMSFileLoad helper(filename, log_type_);
    ConsensusMap loaded;
    helper.load(loaded);
    consensus = std::move(loaded);
  }

  void OMSFile::exportToJSON(const std::string& filename_in, const std::string& filename_out)
  {
    OpenMS::Internal::OMSFileLoad helper(filename_in, log_type_);
    ofstream output(filename_out.c_str());
    if (output.is_open())
    {
      helper.exportToJSON(output);
    }
    else
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_out);
    }
  }

}
