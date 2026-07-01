// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/UniProtXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/UniProtXMLHandler.h>

namespace OpenMS
{
  UniProtXMLFile::UniProtXMLFile() = default;
  UniProtXMLFile::~UniProtXMLFile() = default;

  void UniProtXMLFile::load(const std::string& filename, std::vector<UniProtEntry>& entries)
  {
    entries.clear();
    loadStreaming(filename, [&entries](UniProtEntry&& e) { entries.push_back(std::move(e)); });
  }

  void UniProtXMLFile::loadStreaming(const std::string& filename, const std::function<void(UniProtEntry&&)>& callback)
  {
    Internal::UniProtXMLHandler handler(filename, callback);
    parse_(filename, &handler);
  }
} // namespace OpenMS
