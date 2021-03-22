//
// Created by david on 18.03.21.
//

#include <OpenMS/VISUAL/TVToolDiscovery.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <QProcess>
#include <QCoreApplication>

#include <iostream>

namespace OpenMS {

    std::unordered_map<std::string, std::future<Param>> TVToolDiscovery::future_results_;
    std::unordered_map<std::string, Param> TVToolDiscovery::params_;
    bool TVToolDiscovery::params_ready_ = false;

    void TVToolDiscovery::loadParams() {
      // Get a map of all tools
      const auto& tools = ToolHandler::getTOPPToolList();
      const auto& utils = ToolHandler::getUtilList();
      // Get param for each tool/util
      for (const auto& pair : tools)
      {
        std::string tool_name = pair.first;
        future_results_.insert(
                std::make_pair(tool_name, std::async(std::launch::async, TVToolDiscovery::getParamFromIni_, tool_name))
                );
      }
      for (const auto& pair : utils)
      {
        std::string util_name = pair.first;
        future_results_.insert(
                std::make_pair(util_name, std::async(std::launch::async, TVToolDiscovery::getParamFromIni_, pair.first))
                );
      }
    }

    void TVToolDiscovery::waitForParams()
    {
      if (!params_ready_)
      {
        params_ready_ = true;
        for (auto &pair : future_results_)
        {
          while (pair.second.wait_for(std::chrono::milliseconds(10)) != std::future_status::ready)
          {
            QCoreApplication::processEvents();
          }
          params_.insert(std::make_pair(pair.first, pair.second.get()));
        }
      }
    }

    const std::unordered_map<std::string, Param>& TVToolDiscovery::getToolParams() {
    if (!params_ready_)
    {
      TVToolDiscovery::waitForParams();
    }
    return params_;
    }

    Param TVToolDiscovery::getParamFromIni_(const String &tool_name) {
      String path = File::getUniqueName() + ".ini";
      QStringList args{ "-write_ini", path.toQString()};
      QProcess qp;
      Param tool_param;
      String executable = File::findSiblingTOPPExecutable(tool_name);
      qp.start(executable.toQString(), args);
      const bool success = qp.waitForFinished(-1); // wait till job is finished
      if (qp.error() == QProcess::FailedToStart || !success || qp.exitStatus() != 0 || qp.exitCode() != 0 || !File::exists(path))
      {
        std::remove(path.c_str());
        return tool_param;
      }
      ParamXMLFile paramFile;
      paramFile.load((path).c_str(), tool_param);
      std::remove(path.c_str());
      return tool_param;
    }

}
