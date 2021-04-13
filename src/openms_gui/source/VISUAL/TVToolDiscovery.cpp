// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer:  $
// $Authors: David Voigt $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TVToolDiscovery.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <QProcess>
#include <QCoreApplication>

#include <iostream>

namespace OpenMS {

    void TVToolDiscovery::loadParams()
    {
      // tool params are only loaded once by using a immediately evaluated lambda
      static bool _ [[maybe_unused]] = [&]() -> bool
      {
        // Get a map of all tools
        const auto& tools = ToolHandler::getTOPPToolList();
        const auto& utils = ToolHandler::getUtilList();
        // Launch threads for loading tool/util params.
        for (const auto& pair : tools)
        {
          std::string tool_name = pair.first;
          future_results_.insert(
                  std::make_pair(
                          tool_name,
                          std::async(getParamFromIni_, tool_name)
                          )
          );
        }
        for (const auto& pair : utils)
        {
          std::string util_name = pair.first;
          future_results_.insert(
                  std::make_pair(
                          util_name,
                          std::async(getParamFromIni_, util_name)
                          )
          );
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
        for (auto& pair : future_results_)
        {
          while (pair.second.wait_for(std::chrono::milliseconds(10)) != std::future_status::ready)
          {
            // Keep GUI responsive while waiting
            QCoreApplication::processEvents();
          }
          // Make future results available in params_
          params_.insert(std::make_pair(pair.first, pair.second.get()));
        }
        return true;
      }();
    }

    const std::unordered_map<std::string, Param>& TVToolDiscovery::getToolParams()
    {
      // Make sure threads have been launched and waited for before accessing results
      TVToolDiscovery::loadParams();
      TVToolDiscovery::waitForParams();
      return params_;
    }

    Param TVToolDiscovery::getParamFromIni_(const std::string& tool_name)
    {
      String path = File::getTempDirectory() + "/" + File::getUniqueName() + ".ini";
      QStringList args{ "-write_ini", path.toQString()};

      QProcess qp;
      Param tool_param;
      String executable = File::findSiblingTOPPExecutable(tool_name);
      qp.start(executable.toQString(), args, QProcess::NotOpen);

      const bool success = qp.waitForFinished(-1); // wait till job is finished
      if (qp.error() == QProcess::FailedToStart || !success || qp.exitStatus() != 0 || qp.exitCode() != 0 || !File::exists(path))
      {
        qp.close();
        std::remove(path.c_str());
        return tool_param;
      }
      ParamXMLFile paramFile;
      paramFile.load((path).c_str(), tool_param);

      qp.close();
      std::remove(path.c_str());
      return tool_param;
    }
}
