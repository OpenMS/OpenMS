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
      for (const auto& tool : tools)
      {
        tool_param_futures_.push_back(std::async(std::launch::async, getParamFromIni_, tool.first, false));
      }
      for (const auto& util : utils)
      {
        tool_param_futures_.push_back(std::async(std::launch::async, getParamFromIni_, util.first, false));
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
      // Is an executable already or has a sibling Executable
      executable = File::exists(tool_path) ? tool_path : File::findSiblingTOPPExecutable(tool_path);
    }
    catch (const Exception::FileNotFound& e)
    {
      OPENMS_LOG_WARN << "TOPP tool: " << e << " not found during tool discovery. Skipping." << std::endl;
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
      OPENMS_LOG_WARN << "TOPP tool: " << executable << " error during execution: " << (uint32_t)return_state << "\n";
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
      OPENMS_LOG_WARN << e << "\n" << "TOPP tool: " << executable <<
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
