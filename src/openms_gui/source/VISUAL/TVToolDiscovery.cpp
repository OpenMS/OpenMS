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

#include <iostream>
#include <filesystem>

#include <QCoreApplication>


namespace OpenMS
{



  void TVToolDiscovery::loadParams()
  {
    // tool params are only loaded once by using a immediately evaluated lambda
    static bool _ [[maybe_unused]] = [&]() -> bool
    {
      // Get a map of all tools
      const auto &tools = ToolHandler::getTOPPToolList();
      const auto &utils = ToolHandler::getUtilList();
      const auto &plugins = getPlugins_();
      // Launch threads for loading tool/util params.
      for (auto& tool : tools)
      {
        const std::string name = tool.first;
        param_futures_[name] = std::async(std::launch::async, getParamFromIni_, name);
      }
      for (auto& util : utils)
      {
        const std::string name = util.first;
        param_futures_[name] = std::async(std::launch::async, getParamFromIni_, name);
      }
      for (auto& plugin : plugins)
      {
        std::cout << "starting work on " << plugin << std::endl;
        param_futures_[File::basename(plugin)] = std::async(std::launch::async, getParamFromIni_, plugin);
      }
      return true;
    }();
  }

  void TVToolDiscovery::waitForParams()
  {
    // Make sure that future results are only waited for and inserted in params_ once
    static bool _ [[maybe_unused]] = [&]() -> bool
    {
      // Make sure threads have been launched before waiting
      loadParams();
      // Wait for futures to finish
      for (auto&[name, param_future] : param_futures_)
      {
        while (param_future.wait_for(std::chrono::milliseconds(10)) != std::future_status::ready)
        {
          // Keep GUI responsive while waiting
          QCoreApplication::processEvents();
        }
        // Make future results available in params_
        params_.emplace(name, param_future.get());
      }
      return true;
    }();
  }

  const std::map<std::string, Param> &TVToolDiscovery::getToolParams()
  {
    std::cout << "STARTING" << std::endl;

    // Make sure threads have been launched and waited for before accessing results
    loadParams();
    waitForParams();
    return params_;
  }

  Param TVToolDiscovery::getParamFromIni_(const std::string &tool_name)
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
//    std::cout << "after finding exec " << executable << std::endl;
    // Write tool ini to temporary file
    ExternalProcess proc;
    auto return_state = proc.run(executable.toQString(), args, working_dir.toQString(), true, ExternalProcess::IO_MODE::NO_IO);
    // Return empty param if writing the ini file failed
    if (return_state != ExternalProcess::RETURNSTATE::SUCCESS)
    {
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
        " not able to write ini. Plugins must include -write-ini flag. Skipping." << std::endl;
      return tool_param;
    }


    if (executable.hasSuffix("test.py"))    
    {
      std::cout << "TOOL_PARAMS: " << tool_param << std::endl;
    }

    return tool_param;
  }

  const ToolListType &TVToolDiscovery::getPlugins()
  {
    ToolListType tool_map;
    std::map<std::string, Param> params;
    StringList plugins;

    params = TVToolDiscovery::getToolParams();
    plugins = TVToolDiscovery::getPlugins_();

    for (auto param : params)
    {
      std::cout << param.first << std::endl;
    }
  
    return tool_map;
  }

  const StringList TVToolDiscovery::getPlugins_()
  {
    StringList plugins;

    std::cout << "PLUGIN DETECTION" << std::endl; 
   
    std::vector<std::string> valid_extensions {".py"};
    const std::string plugin_path = File::absolutePath("./test/");

    if (File::fileList(plugin_path, "*", plugins, true))
    {
      const auto comparator = [plugin_path, valid_extensions](std::string plugin) -> bool
      {
        return !File::executable(plugin) /*&& 
          (std::find(valid_extensions.begin(), valid_extensions.end(), std::filesystem::path(plugin).extension()) == valid_extensions.end())*/; 
      };
      plugins.erase(std::remove_if(plugins.begin(), plugins.end(), comparator), plugins.end());
    }

    for (auto p : plugins) 
    {
      std::cout << "plugin " << p << std::endl; 
    }

    return plugins;
  }
}
