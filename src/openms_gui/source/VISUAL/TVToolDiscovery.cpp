// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: David Voigt $
// $Authors: David Voigt $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TVToolDiscovery.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <iostream>
#include <filesystem>

#include <QCoreApplication>


namespace OpenMS
{
  void TVToolDiscovery::loadToolParams()
  {
    // tool params are only loaded once by using a immediately evaluated lambda
    static bool _ [[maybe_unused]] = [&]() -> bool
    {
      // Get a map of all tools
      const auto &tools = ToolHandler::getTOPPToolList();
      const auto &utils = ToolHandler::getUtilList();
      // Launch threads for loading tool/util params.
      for (const auto& [name, _] : tools)
      {
        tool_param_futures_.push_back(std::async(std::launch::async, getParamFromIni_, name, nullptr));
      }
      for (const auto& [name, _] : utils)
      {
        tool_param_futures_.push_back(std::async(std::launch::async, getParamFromIni_, name, nullptr));
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
      std::cout << "starting work on " << plugin << std::endl;
      plugin_param_futures_.push_back(std::async(std::launch::async, getParamFromIni_, plugin, &plugins_));
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
        // Make future results available in params_
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
      // Make future results available in params_
      plugin_params_.insert("", param_future.get());
    }
  }

  //const std::vector<std::pair<String, Param>> &TVToolDiscovery::getToolParams()
  const Param& TVToolDiscovery::getToolParams()
  {
    std::cout << "STARTING" << std::endl;

    // Make sure threads have been launched and waited for before accessing results
    loadToolParams();
    waitForToolParams();
    return tool_params_;
  }

  //const std::vector<std::pair<String, Param>> &TVToolDiscovery::getPluginParams()
  const Param& TVToolDiscovery::getPluginParams()
  {
    plugin_params_.clear();
    waitForPluginParams();
    return plugin_params_;
  }

  Param TVToolDiscovery::getParamFromIni_(const String &tool_name, std::vector<std::pair<String, String>> *plugins)
  {
    FileHandler fh;
    // Temporary file path and arguments
    String path = File::getTemporaryFile();
    String working_dir = path.prefix(path.find_last_of('/'));
    QStringList args{"-write_ini", path.toQString()};
//    std::cout << "before finding exec " << tool_name << std::endl;
    Param tool_param;
    String executable;
    // Return empty param if tool executable cannot be found
    try
    {
      // Is an executable already or has a sibling Executable
      executable = File::exists(tool_name) ? tool_name : File::findSiblingTOPPExecutable(tool_name);
    }
    catch (const Exception::FileNotFound& e)
    {
      std::cerr << "TOPP tool: " << e << " not found during tool discovery. Skipping." << std::endl;
      return tool_param;
    }

    // Write tool ini to temporary file
    auto lam_out = [&](const String& out) { OPENMS_LOG_INFO << out; };
    auto lam_err = [&](const String& out) { OPENMS_LOG_INFO << out; };
    ExternalProcess proc(lam_out, lam_err);
    auto return_state = proc.run(executable.toQString(), args, working_dir.toQString(), true, ExternalProcess::IO_MODE::NO_IO);
    // Return empty param if writing the ini file failed
    if (return_state != ExternalProcess::RETURNSTATE::SUCCESS)
    {
      std::cerr << "TOPP tool: " << executable << " error during execution: " << (uint32_t)return_state << "\n";
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
      std::cerr << e << "\n" << "TOPP tool: " << executable << 
        " not able to write ini. Plugins must implement -write-ini flag. Skipping." << std::endl;
      return tool_param;
    }
    
    if (plugins) {
      auto param_line = tool_param.begin().getName();
      std::cout << "PLUGIN NAME: " << param_line.substr(0, param_line.find_first_of(':')) << std::endl;
     
      // is this a problem for thread safety that we use plugins?
      // fill list of plugins with the Plugin Name from the ini
      //plugins->emplace(plugins->end(), param_line.substr(0, param_line.find_first_of(':')));
      //We also save the actual filename of the tool, we need that to execute it later on
      //this way the tool name specified in the ini can be something different
      //otherwise it would be necessary to be EXACTLY the same, including file extension
      plugins->emplace_back(tool_name.suffix(tool_name.size() - tool_name.find_last_of('/') - 1),
                            param_line.substr(0, param_line.find_first_of(':')));
      std::cout << "PLUGINS LIST ENTRY 1: " << plugins->at(0).second << std::endl;
    }

    return tool_param;
  }

  // MAYBE USE THIS TO GET THE PLUGINS IN TOOLSDIALOG
  const std::vector<std::pair<String, String>> &TVToolDiscovery::getPlugins()
  {
    return plugins_;
  }

  const StringList TVToolDiscovery::getPlugins_()
  {
    StringList plugins;

    std::cout << "PLUGIN DETECTION" << std::endl; 
    // this is unused right now... change the return value in the comparator to use this
    std::vector<std::string> valid_extensions {".py"};

    if (File::fileList(plugin_path_, "*", plugins, true))
    {
      const auto comparator = [valid_extensions](const std::string& plugin) -> bool
      {
        return !File::executable(plugin) /*&& 
          (std::find(valid_extensions.begin(), valid_extensions.end(), std::filesystem::path(plugin).extension()) == valid_extensions.end())*/; 
      };
      plugins.erase(std::remove_if(plugins.begin(), plugins.end(), comparator), plugins.end());
    }

    // this is just for debugging
    for (auto& p : plugins)
    {
      std::cout << "plugin " << p << "\n";
    }

    std::cout << "END PLUGIN DETECTION" << std::endl;

    return plugins;
  }

  void TVToolDiscovery::setPluginPath(const std::string &path)
  {
    plugin_path_ = path;
  }

  const std::string &TVToolDiscovery::getPluginPath()
  {
    return plugin_path_;
  }

}
