// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ToolJSONFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace OpenMS
{
  bool ToolJSONFile::load(const String& filename, 
                          std::vector<Internal::ToolDescription>& tools,
                          std::map<String, String>& categories)
  {
    tools.clear();
    categories.clear();

    std::ifstream ifs(filename);
    if (!ifs.good())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    try
    {
      json doc = json::parse(ifs);

      // Load categories
      if (doc.contains("categories") && doc["categories"].is_object())
      {
        for (const auto& [key, value] : doc["categories"].items())
        {
          if (value.is_string())
          {
            categories[key] = value.get<std::string>();
          }
          else
          {
            OPENMS_LOG_WARN << "Skipping non-string category value for key: " << key << std::endl;
          }
        }
      }
      else
      {
        OPENMS_LOG_WARN << "No 'categories' section found in JSON file: " << filename << std::endl;
      }

      // Load tools
      if (doc.contains("tools") && doc["tools"].is_object())
      {
        for (const auto& [tool_name, tool_config] : doc["tools"].items())
        {
          if (!tool_config.is_object())
          {
            OPENMS_LOG_WARN << "Skipping non-object tool configuration for: " << tool_name << std::endl;
            continue;
          }

          // Check build-time requirements
          bool skip_tool = false;

          // Check GUI requirement
          if (tool_config.contains("requires_gui") && tool_config["requires_gui"].is_boolean() && 
              tool_config["requires_gui"].get<bool>())
          {
#ifndef WITH_GUI
            skip_tool = true;
#endif
          }

          // Check OpenSWATH requirement
          if (tool_config.contains("requires_openswath") && tool_config["requires_openswath"].is_boolean() && 
              tool_config["requires_openswath"].get<bool>())
          {
#ifdef DISABLE_OPENSWATH
            skip_tool = true;
#endif
          }

          // Check Parquet requirement
          if (tool_config.contains("requires_parquet") && tool_config["requires_parquet"].is_boolean() && 
              tool_config["requires_parquet"].get<bool>())
          {
#ifndef WITH_PARQUET
            skip_tool = true;
#endif
          }

          if (skip_tool)
          {
            continue;
          }

          // Get category
          String category;
          if (tool_config.contains("category") && tool_config["category"].is_string())
          {
            String category_id = tool_config["category"].get<std::string>();
            auto it = categories.find(category_id);
            if (it != categories.end())
            {
              category = it->second;
            }
            else
            {
              OPENMS_LOG_WARN << "Unknown category ID '" << category_id << "' for tool: " << tool_name << std::endl;
              category = category_id; // Use the ID as fallback
            }
          }
          else
          {
            OPENMS_LOG_WARN << "No category specified for tool: " << tool_name << std::endl;
            category = "Misc"; // Default category
          }

          // Create tool description
          Internal::ToolDescription tool(tool_name, category);
          tools.push_back(tool);
        }
      }
      else
      {
        OPENMS_LOG_WARN << "No 'tools' section found in JSON file: " << filename << std::endl;
      }
    }
    catch (const json::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", e.what());
    }
    catch (const std::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", e.what());
    }

    return true;
  }

  String ToolJSONFile::getDefaultConfigPath()
  {
    return File::getOpenMSDataPath() + "/TOOLS/tools.json";
  }

} // namespace OpenMS