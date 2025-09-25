// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <OpenMS/FORMAT/ToolDescriptionFile.h>
#include <OpenMS/FORMAT/ToolJSONFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <QStringList>
#include <QtCore/QDir>

namespace OpenMS
{
  ToolListType ToolHandler::getTOPPToolList(const bool includeGenericWrapper)
  {
    ToolListType tools_map;

    // Try to load from JSON first
    if (loadFromJSON(true) && tools_json_available_)
    {
      // Use JSON-based tool definitions
      for (const auto& tool : tools_json_)
      {
        tools_map[tool.name] = tool;
      }

      // Apply the same includeGenericWrapper logic as the hardcoded fallback
      if (!includeGenericWrapper)
      {
        tools_map.erase("GenericWrapper");
      }
    }

    // INTERNAL tools
    // this operation is expensive, as we need to parse configuration files (*.ttd)
    std::vector<Internal::ToolDescription> internal_tools = getInternalTools_();
    for (std::vector<Internal::ToolDescription>::const_iterator it = internal_tools.begin(); it != internal_tools.end(); ++it)
    {
      if (tools_map.find(it->name) == tools_map.end())
      {
        tools_map[it->name] = *it;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Duplicate tool name error: Trying to add internal tool '" + it->name, it->name);
      }
    }

    // EXTERNAL tools
    // this operation is expensive, as we need to parse configuration files (*.ttd)
    if (includeGenericWrapper)
    {
      tools_map["GenericWrapper"] = getExternalTools_();
    }

    return tools_map;
  }

  StringList ToolHandler::getTypes(const String& toolname)
  {
    Internal::ToolDescription ret;
    ToolListType tools;
    if (toolname == "GenericWrapper")
    {
      tools = getTOPPToolList(true);
    }
    else
    {
      tools = getTOPPToolList();
    }
    if (tools.find(toolname) != tools.end())
    {
      return tools[toolname].types;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Requested tool '" + toolname + "' does not exist!", toolname);
  }

  std::vector<Internal::ToolDescription> ToolHandler::getInternalTools_()
  {
    if (!tools_internal_loaded_)
    {
      loadInternalToolConfig_();
      tools_internal_loaded_ = true;
    }
    return tools_internal_;
  }

  String ToolHandler::getExternalToolsPath()
  {
    return File::getOpenMSDataPath() + "/TOOLS/EXTERNAL";
  }

  String ToolHandler::getInternalToolsPath()
  {
    return File::getOpenMSDataPath() + "/TOOLS/INTERNAL";
  }

  Internal::ToolDescription ToolHandler::getExternalTools_()
  {
    if (!tools_external_loaded_)
    {
      loadExternalToolConfig_();
      tools_external_loaded_ = true;
    }

    return tools_external_;
  }

  void ToolHandler::loadExternalToolConfig_()
  {
    QStringList files = getExternalToolConfigFiles_();
    for (int i = 0; i < files.size(); ++i)
    {
      ToolDescriptionFile tdf;
      std::vector<Internal::ToolDescription> tools;
      tdf.load(String(files[i]), tools);
      // add every tool from file to list
      for (Size i_t = 0; i_t < tools.size(); ++i_t)
      {
        if (i == 0 && i_t == 0)
          {
            tools_external_ = tools[i_t]; // init
          }
        else
          {
            tools_external_.append(tools[i_t]); // append
          }
      }
    }
    tools_external_.name = "GenericWrapper";
    tools_external_.category = "EXTERNAL";
  }

  void ToolHandler::loadInternalToolConfig_()
  {
    QStringList files = getInternalToolConfigFiles_();
    for (int i = 0; i < files.size(); ++i)
    {
      ToolDescriptionFile tdf;
      std::vector<Internal::ToolDescription> tools;
      tdf.load(String(files[i]), tools);
      // add every tool from file to list
      for (Size i_t = 0; i_t < tools.size(); ++i_t)
      {
        tools_internal_.push_back(tools[i_t]);
        tools_external_.category = "INTERNAL";
      }
    }
  }

  QStringList ToolHandler::getExternalToolConfigFiles_()
  {

    QStringList paths;
    // *.ttd default path
    paths << getExternalToolsPath().toQString();
    // OS-specific path
#ifdef OPENMS_WINDOWSPLATFORM
    paths << (getExternalToolsPath() + "/WINDOWS").toQString();
#else
    paths << (getExternalToolsPath() + "/LINUX").toQString();
#endif
    // additional environment
    if (getenv("OPENMS_TTD_PATH") != nullptr)
    {
      paths << String(getenv("OPENMS_TTD_PATH")).toQString();
    }

    QStringList all_files;
    for (int p = 0; p < paths.size(); ++p)
    {
      QDir dir(paths[p], "*.ttd");
      QStringList files = dir.entryList();
      for (int i = 0; i < files.size(); ++i)
      {
        files[i] = dir.absolutePath() + QDir::separator() + files[i];
      }
      all_files << files;
    }
    //StringList list = ListUtils::create<String>(getExternalToolsPath() + "/" + "msconvert.ttd");
    return all_files;
  }

  QStringList ToolHandler::getInternalToolConfigFiles_()
  {
    QStringList paths;
    // *.ttd default path
    paths << getInternalToolsPath().toQString();
    // OS-specific path
#ifdef OPENMS_WINDOWSPLATFORM
    paths << (getInternalToolsPath() + "/WINDOWS").toQString();
#else
    paths << (getInternalToolsPath() + "/LINUX").toQString();
#endif
    // additional environment
    if (getenv("OPENMS_TTD_INTERNAL_PATH") != nullptr)
    {
      paths << String(getenv("OPENMS_TTD_INTERNAL_PATH")).toQString();
    }

    QStringList all_files;
    for (int p = 0; p < paths.size(); ++p)
    {
      QDir dir(paths[p], "*.ttd");
      QStringList files = dir.entryList();
      for (int i = 0; i < files.size(); ++i)
      {
        files[i] = dir.absolutePath() + QDir::separator() + files[i];
      }
      all_files << files;
    }
    return all_files;
  }

  String ToolHandler::getCategory(const String& toolname)
  {
    ToolListType tools = getTOPPToolList(true);
    String s;
    if (tools.find(toolname) != tools.end())
    {
      s = tools[toolname].category;
    }

    return s;
  }

  bool ToolHandler::loadFromJSON(const bool use_json)
  {
    if (!use_json)
    {
      tools_json_available_ = false;
      return false;
    }

    if (tools_json_loaded_)
    {
      return tools_json_available_;
    }

    try
    {
      String json_path = ToolJSONFile::getDefaultConfigPath();
      if (File::exists(json_path))
      {
        tools_json_available_ = ToolJSONFile::load(json_path, tools_json_, categories_json_);
        tools_json_loaded_ = true;
        return tools_json_available_;
      }
      else
      {
        OPENMS_LOG_INFO << "JSON tool configuration file not found at: " << json_path << std::endl;
        tools_json_available_ = false;
        tools_json_loaded_ = true;
        return false;
      }
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "Failed to load JSON tool configuration: " << e.what() << std::endl;
      tools_json_available_ = false;
      tools_json_loaded_ = true;
      return false;
    }
  }

  // static
  Internal::ToolDescription ToolHandler::tools_external_ = Internal::ToolDescription();
  std::vector<Internal::ToolDescription> ToolHandler::tools_internal_;
  bool ToolHandler::tools_external_loaded_ = false;
  bool ToolHandler::tools_internal_loaded_ = false;

  // JSON-based static members
  std::vector<Internal::ToolDescription> ToolHandler::tools_json_;
  std::map<String, String> ToolHandler::categories_json_;
  bool ToolHandler::tools_json_loaded_ = false;
  bool ToolHandler::tools_json_available_ = false;

} // namespace
