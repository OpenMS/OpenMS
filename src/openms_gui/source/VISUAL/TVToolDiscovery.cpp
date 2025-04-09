// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: David Voigt $
// $Authors: David Voigt, Ruben Grünberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TVToolDiscovery.h>

#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <QCoreApplication>
#include <QDir>

#include <thread>

namespace OpenMS
{
  void TVToolDiscovery::loadToolParams()
  {
    // tool params are only loaded once by using a immediately evaluated lambda
    static bool _ [[maybe_unused]] = [&]() -> bool
    {
      // Get a map of all tools
      const auto &tools = ToolHandler::getTOPPToolList();
      // Launch threads for loading tool/util params.
      for (const auto& tool : tools)
      {
        tool_param_futures_.push_back(std::async(std::launch::async, getParamFromIni_, tool.first, false));
      }
      return true;
    }();
  }

  void TVToolDiscovery::loadPluginParams()
  {
    plugin_param_futures_.clear();
    plugins_.clear();
    const auto &plugins = getPlugins_();
    for (auto& plugin : plugins)
    {
      plugin_param_futures_.push_back(std::async(std::launch::async, getParamFromIni_, plugin, true));
    }
  }

  void TVToolDiscovery::waitForToolParams()
  {
    // Make sure that future results are only waited for and inserted in params_ once
    static bool _ [[maybe_unused]] = [&]() -> bool
    {
      // Make sure threads have been launched before waiting
      loadToolParams();
      // Wait for futures to finish
      for (auto& param_future : tool_param_futures_)
      {
        while (param_future.wait_for(std::chrono::milliseconds(10)) != std::future_status::ready)
        {
          // Keep GUI responsive while waiting
          QCoreApplication::processEvents();
        }
        // Make future results available in tool_params_
        tool_params_.insert("", param_future.get());
      }
      return true;
    }();
  }

  void TVToolDiscovery::waitForPluginParams()
  {
    // Make sure threads have been launched before waiting
    loadPluginParams();
    // Wait for futures to finish
    for (auto& param_future : plugin_param_futures_)
    {
      while (param_future.wait_for(std::chrono::milliseconds(10)) != std::future_status::ready)
      {
        // Keep GUI responsive while waiting
        QCoreApplication::processEvents();
      }
      // Make future results available in plugin_params_
      Param new_param = param_future.get();
      // Skip if the param is empty, that means something went wrong during execution
      if (new_param.empty()) continue;
      plugins_.push_back(new_param.begin().getTrace().begin()->name);
      plugin_params_.insert("", new_param);
    }
  }

  const Param& TVToolDiscovery::getToolParams()
  {
    // Make sure threads have been launched and waited for before accessing results
    waitForToolParams();
    return tool_params_;
  }

  const Param& TVToolDiscovery::getPluginParams()
  {
    plugin_params_.clear();
    waitForPluginParams();
    return plugin_params_;
  }

  Param TVToolDiscovery::getParamFromIni_(const String& tool_path, bool plugins)
  {
    static std::mutex io_mutex;
    FileHandler fh;
    // Temporary file path and arguments
    String path = File::getTemporaryFile();
    String working_dir = path.prefix(path.find_last_of('/'));
    QStringList args{"-write_ini", path.toQString()};
    Param tool_param;
    String executable;
    // Return empty param if tool executable cannot be found
    try
    {
      std::scoped_lock lock(io_mutex);
      // Is an executable already or has a sibling Executable
      executable = File::exists(tool_path) ? tool_path : File::findSiblingTOPPExecutable(tool_path);
    }
    catch (const Exception::FileNotFound& e)
    {
      std::scoped_lock lock(io_mutex);
      OPENMS_LOG_DEBUG << "TOPP tool: " << e << " not found during tool discovery. Skipping." << std::endl;
      return tool_param;
    }

    // Write tool ini to temporary file
    static std::atomic<int> running_processes{0}; // used to limit the number of parallel processes
    auto lam_out = [&](const String& out) { OPENMS_LOG_INFO << out; };
    auto lam_err = [&](const String& out) { OPENMS_LOG_INFO << out; };

    // Spawning a thread for all tools is no problem (if std::async decides to do so)
    // but spawning that many processes failed with not enough file handles on machines with large number of cores.
    // Restricting the number of running processes solves that issue.
    while (running_processes >= 6) 
    { 
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      QCoreApplication::processEvents();
    }

    ExternalProcess proc(lam_out, lam_err);
    // Write tool ini to temporary file
    ++running_processes;
    auto return_state = proc.run(executable.toQString(), args, working_dir.toQString(), true, ExternalProcess::IO_MODE::NO_IO);
    --running_processes;

    // Return empty param if writing the ini file failed
    if (return_state != ExternalProcess::RETURNSTATE::SUCCESS)
    {
      std::scoped_lock lock(io_mutex);
      OPENMS_LOG_DEBUG << "TOPP tool: " << executable << " error during execution: " << (uint32_t)return_state << "\n";
      return tool_param;
    }
    // Parse ini file to param object
    ParamXMLFile paramFile;
    try
    {
      paramFile.load((path).c_str(), tool_param);
    }
    catch(const Exception::FileNotFound& e)
    {
      std::scoped_lock lock(io_mutex);
      OPENMS_LOG_DEBUG << e << "\n" << "TOPP tool: " << executable <<
        " not able to write ini. Plugins must implement -write_ini parameter. Skipping." << std::endl;
      return tool_param;
    }
    
    if (plugins)
    {
      auto tool_name = tool_param.begin().getTrace().begin()->name;
      auto filename = File::basename(tool_path);
      tool_param.setValue(tool_name + ":filename", filename, "The filename of the plugin executable. This entry is automatically generated.");
    }

    return tool_param;
  }

  const std::vector<std::string>& TVToolDiscovery::getPlugins()
  {
    return plugins_;
  }

  const StringList TVToolDiscovery::getPlugins_()
  {
    StringList plugins;

    // here all supported file extensions can be added
    std::vector<std::string> valid_extensions {"", ".py"};
    const auto comparator = [valid_extensions](const std::string& plugin) -> bool
    {
        return !File::executable(plugin) ||
          (std::find(valid_extensions.begin(), valid_extensions.end(), plugin.substr(plugin.find_last_of('.'))) == valid_extensions.end());
    };

    if (File::fileList(plugin_path_, "*", plugins, true))
    {
      plugins.erase(std::remove_if(plugins.begin(), plugins.end(), comparator), plugins.end());
    }

    return plugins;
  }

  bool TVToolDiscovery::setPluginPath(const String& plugin_path, bool create)
  {
    if (!File::exists(plugin_path))
    {
      if (create)
      {
        QDir path = QDir(plugin_path.toQString());
        QString dir = path.dirName();
        path.cdUp();

        if (!path.mkdir(dir))
        {
          OPENMS_LOG_WARN << "Unable to create plugin directory " << plugin_path << std::endl;
          //plugin_path_ = plugin_path;
          return false;
        }
      }
      else
      {
        OPENMS_LOG_WARN << "Unable to set plugin directory: " << plugin_path << " does not exist." << std::endl;
        return false;
      }
    }

    plugin_path_ = plugin_path;
    return true;
  }

  const std::string TVToolDiscovery::getPluginPath()
  {
    return plugin_path_;
  }

  void TVToolDiscovery::setVerbose(int verbosity_level)
  {
    verbosity_level_ = verbosity_level;
  }

  std::string TVToolDiscovery::findPluginExecutable(const std::string& name)
  {
    //TODO: At the moment the usage of subdirectories in the plugin path are not possible
    //To change that, the tool scanner has to recursively search all directories in the plugin path
    if (!plugin_params_.exists(name + ":filename"))
    {
      return "";
    }
    return plugin_path_ + "/" + plugin_params_.getValue(name + ":filename").toString();
  }

}
