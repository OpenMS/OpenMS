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
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

#include <QCoreApplication>

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
        tool_param_futures_.push_back(std::async(std::launch::async, getParamFromIni_, tool.first));
      }
      return true;
    }();
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

  const Param& TVToolDiscovery::getToolParams()
  {
    // Make sure threads have been launched and waited for before accessing results
    waitForToolParams();
    return tool_params_;
  }

  Param TVToolDiscovery::getParamFromIni_(const std::string& tool_path)
  {
    static std::mutex io_mutex;
    FileHandler fh;
    // Temporary file path and arguments
    std::string path = File::getTemporaryFile();
    std::string working_dir = StringUtils::prefix(path, path.find_last_of('/'));
    std::vector<std::string> args{"-write_ini", path};
    Param tool_param;
    std::string executable;
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
    auto lam_out = [&](const std::string& out) { OPENMS_LOG_INFO << out; };
    auto lam_err = [&](const std::string& out) { OPENMS_LOG_INFO << out; };

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
    auto return_state = proc.run(executable, args, working_dir, true, ExternalProcess::IO_MODE::NO_IO);
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
        " not able to write ini (-write_ini failed). Skipping." << std::endl;
      return tool_param;
    }

    return tool_param;
  }

  void TVToolDiscovery::setVerbose(int verbosity_level)
  {
    verbosity_level_ = verbosity_level;
  }

}
