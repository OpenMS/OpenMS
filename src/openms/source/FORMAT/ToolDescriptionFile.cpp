// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ToolDescriptionFile.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/HANDLERS/ToolDescriptionHandler.h>

namespace OpenMS
{

  ToolDescriptionFile::ToolDescriptionFile() :
    XMLFile("/SCHEMAS/ToolDescriptor_1_0.xsd", "1.0.0")
  {
  }

  ToolDescriptionFile::~ToolDescriptionFile() = default;

  void ToolDescriptionFile::load(const String & filename, std::vector<Internal::ToolDescription> & tds)
  {
    Internal::ToolDescriptionHandler handler(filename, schema_version_);
    parse_(filename, &handler);
    tds = handler.getToolDescriptions();
  }

  void ToolDescriptionFile::store(const String & filename, const std::vector<Internal::ToolDescription> & tds) const
  {
    Internal::ToolDescriptionHandler handler(filename, schema_version_);
    handler.setToolDescriptions(tds);
    save_(filename, &handler);
  }

} // namespace OpenMS
