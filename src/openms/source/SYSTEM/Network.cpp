// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/Network.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/NetworkGetRequest.h>

#include <filesystem>
#include <fstream>

using namespace std;

namespace OpenMS
{

namespace
{
  // construct a filename from URL. Add number suffix if file already exists in dest_folder.
  std::string saveFileName_(const std::string& url, const std::string& dest_folder)
  {
    namespace fs = std::filesystem;
    // extract filename from URL path
    std::string path_part = url;
    auto query_pos = path_part.find('?');
    if (query_pos != std::string::npos) path_part = path_part.substr(0, query_pos);
    auto frag_pos = path_part.find('#');
    if (frag_pos != std::string::npos) path_part = path_part.substr(0, frag_pos);

    std::string basename = fs::path(path_part).filename().string();
    if (basename.empty()) basename = "download";

    fs::path dest(dest_folder);
    if (!fs::exists(dest / basename)) return basename;

    // already exists, don't overwrite
    int i = 0;
    while (fs::exists(dest / (basename + "." + std::to_string(i))))
      ++i;
    return basename + "." + std::to_string(i);
  }
} // anonymous namespace

void Network::downloadFile(const std::string& url, const std::string& download_folder)
{
  NetworkGetRequest query;
  query.setUrl(url);
  query.setTimeout(600); // 10 minutes timeout
  query.run();

  if (!query.hasError())
  {
    std::string folder = download_folder.empty() ? "./" : download_folder;
    std::string filename = folder + "/" + saveFileName_(url, folder);
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs)
    {
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Failed to open output file: " + filename);
    }
    const auto& data = query.getResponseBinary();
    ofs.write(data.data(), static_cast<std::streamsize>(data.size()));
    if (!ofs)
    {
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Failed to write downloaded data to: " + filename);
    }
    ofs.close();
    OPENMS_LOG_INFO << "Download of '" << url << "' successful." << endl;
    OPENMS_LOG_INFO << "Stored as '" << filename << "'." << endl;
  }
  else
  {
    String error = "Download of '" + url + "' failed!. Error: " + query.getErrorString() + '\n';
    throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
  }
}

} // namespace OpenMS
