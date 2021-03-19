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

    std::vector<std::future<Param>> TVToolDiscovery::results_;

    void TVToolDiscovery::findTools() {
      // Get a map of all tools
      const auto& tools = ToolHandler::getTOPPToolList();
      const auto& utils = ToolHandler::getUtilList();
      // Get param for each tool/util
      for (const auto& pair : tools)
      {
        if (!pair.first.empty())
        results_.push_back(std::move(std::async(std::launch::async, [&]() { return TVToolDiscovery::getParamFromIni_(pair.first); })));
      }
      for (const auto& pair : utils)
      {
        if (!pair.first.empty())
        results_.push_back(std::move(std::async(std::launch::async, [&]() { return TVToolDiscovery::getParamFromIni_(pair.first); })));
      }
    }

    std::vector<Param> TVToolDiscovery::getToolParams() {
      std::vector<Param> res;
      for (auto& r : results_)
      {
        while(r.wait_for(std::chrono::milliseconds(10)) != std::future_status::ready)
        {
          QCoreApplication::processEvents();
        }
        res.push_back(r.get());
      }
      return res;
    }

    Param TVToolDiscovery::getParamFromIni_(const String &tool_name) {
      String path = File::getUniqueName() + ".ini";
      QStringList args{ "-write_ini", path.toQString()};
      QProcess qp;
      Param tool_param;
      String executable = File::findSiblingTOPPExecutable(tool_name);
      std::cout << executable << "\t" << tool_name << std::endl;
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
