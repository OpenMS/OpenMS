// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// $Maintainer: Timo Sachsenberg $
// $Authors: Johannes Junker, Chris Bielow, Timo Sachsenberg $

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineExecutor.h>
#include <OpenMS/APPLICATIONS/PIPELINE/PipelineFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <OpenMS/SYSTEM/SystemSettings.h>
#include <atomic>
#include <chrono>
#include <csignal>
#include <filesystem>
#include <iostream>
#include <thread>

using namespace OpenMS;

namespace
{
std::atomic_flag interrupted = ATOMIC_FLAG_INIT;
void requestInterrupt(int signal)
{
  interrupted.test_and_set(std::memory_order_relaxed);
  // Cancellation stops the children; a second signal terminates at once if it hangs.
  std::signal(signal, SIG_DFL);
}

struct SignalGuard
{
  using Handler = void (*)(int);
  Handler interrupt_handler;
  Handler terminate_handler;
  SignalGuard()
  {
    interrupted.clear(std::memory_order_relaxed);
    interrupt_handler = install(SIGINT);
    terminate_handler = install(SIGTERM);
  }
  static Handler install(int signal)
  {
    const Handler previous = std::signal(signal, requestInterrupt);
    // Keep signals ignored that the caller ignores (nohup, background jobs of scripts).
    if (previous == SIG_IGN) { std::signal(signal, SIG_IGN); }
    return previous;
  }
  ~SignalGuard()
  {
    std::signal(SIGINT, interrupt_handler);
    std::signal(SIGTERM, terminate_handler);
  }
};

// Keep the monitor joinable on all supported standard libraries, including
// Xcode 16's libc++ where std::jthread is still an experimental feature.
struct InterruptMonitor
{
  std::atomic_bool stopped {false};
  std::thread worker;

  explicit InterruptMonitor(PipelineExecutor& executor):
      worker([this, &executor] {
        while (! stopped.load(std::memory_order_relaxed))
        {
          if (interrupted.test(std::memory_order_relaxed))
          {
            executor.cancel();
            return;
          }
          std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
      })
  {
  }

  ~InterruptMonitor()
  {
    stopped.store(true, std::memory_order_relaxed);
    worker.join();
  }
};
} // namespace

// clang-format off
/**
@page TOPP_ExecutePipeline ExecutePipeline

@brief Executes workflows created by TOPPAS.

This tool is the non-GUI, i.e. command line version for non-interactive execution of TOPPAS pipelines.
In order to really use this tool in batch-mode, you can provide a TOPPAS resource file (.trf) which specifies the
input files for the input nodes in your pipeline.

<B> *.trf files </B>

A TOPPAS resource file (<TT>*.trf</TT>) specifies the locations of input files for a pipeline.
It is an XML file following the normal TOPP INI file schema, i.e. it can be edited using the INIFileEditor or filled using a script (we do NOT provide
one - sorry). It can be exported from TOPPAS (<TT>File -> Save TOPPAS resource file</TT>). For two input nodes 1 and 2 with files
(<TT>dataA.mzML</TT>, <TT>dataB.mzML</TT>) and (<TT>dataC.mzML</TT>) respectively it has the following format.

\code
<?xml version="1.0" encoding="ISO-8859-1"?>
<PARAMETERS version="1.3" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/Param_1_3.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <NODE name="1" description="">
    <ITEMLIST name="url_list" type="string" description="">
      <LISTITEM value="file:///Users/jeff/dataA.mzML"/>
      <LISTITEM value="file:///Users/jeff/dataB.mzML"/>
    </ITEMLIST>
  </NODE>
  <NODE name="2" description="">
    <ITEMLIST name="url_list" type="string" description="">
      <LISTITEM value="file:///Users/jeff/dataC.mzML"/>
    </ITEMLIST>
  </NODE>
</PARAMETERS>
\endcode

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ExecutePipeline.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ExecutePipeline.html
*/
// clang-format on

/// @cond TOPPCLASSES
class TOPPExecutePipeline : public TOPPBase
{
public:
  TOPPExecutePipeline(): TOPPBase("ExecutePipeline", "Executes workflows created by TOPPAS.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "The workflow to be executed.");
    setValidFormats_("in", {"toppas"});
    registerStringOption_("out_dir", "<directory>", "", "Directory for output files (default: user's home directory)", false);
    registerStringOption_("resource_file", "<file>", "", "A TOPPAS resource file (*.trf) specifying the files this workflow is to be applied to",
                          false);
    registerIntOption_("num_jobs", "<integer>", 1, "Maximum number of jobs running in parallel", false, false);
    setMinInt_("num_jobs", 1);
  }

  ExitCodes main_(int, const char**) override
  {
    try
    {
      return execute_();
    }
    catch (const Exception::BaseException&)
    {
      throw; // reported by TOPPBase
    }
    catch (const std::filesystem::filesystem_error& error)
    {
      OPENMS_LOG_ERROR << "Error: " << error.what() << '\n';
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    catch (const std::exception& error)
    {
      OPENMS_LOG_ERROR << "Error: " << error.what() << '\n';
      return UNKNOWN_ERROR;
    }
  }

  ExitCodes execute_()
  {
    namespace fs = std::filesystem;
    const auto filename = getStringOption_("in");
    PipelineGraph graph;
    PipelineFile().load(filename, graph);
    const auto resources = getStringOption_("resource_file");
    if (! resources.empty()) PipelineFile().loadResources(resources, graph);

    auto output = getStringOption_("out_dir");
    if (output.empty())
    {
      const auto basename = File::basename(filename);
      // Join natively: the home directory may be in the Windows ANSI code page and already ends in '/'.
      const auto directory = to_path(SystemSettings::getUserDirectory()) / to_path(basename.substr(0, basename.find('.')));
      const auto utf8 = directory.generic_u8string();
      output.assign(utf8.begin(), utf8.end());
      std::cout << "No output directory specified. Using " << output << '\n';
    }
    else if (std::error_code error; ! fs::is_directory(to_path(output), error))
    {
      OPENMS_LOG_ERROR << "The specified output directory does not exist: " << output << '\n';
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    PipelineExecutor::Options options;
    options.output_directory = output;
    options.temp_directory = SystemSettings::getTempDirectory();
    options.num_jobs = static_cast<Size>(getIntOption_("num_jobs"));

    PipelineExecutor executor;
    SignalGuard signals;
    // atomic_flag is guaranteed lock-free and can safely communicate from a
    // signal handler to the monitor thread. Termination happens outside the handler.
    InterruptMonitor monitor(executor);
    const auto result = executor.run(graph, options, [](const PipelineExecutor::Event& event) {
      if (! event.text.empty())
      {
        auto& stream = event.type == PipelineExecutor::Event::Type::NODE_FAILED ? std::cerr : std::cout;
        stream << event.text;
        if (event.text.back() != '\n') stream << '\n';
      }
    });
    // Child exit codes are passed on; ExitCodes cannot represent other values.
    if (result.exit_code < EXECUTION_OK || result.exit_code > EXTERNAL_PROGRAM_NOTFOUND) { return EXTERNAL_PROGRAM_ERROR; }
    return static_cast<ExitCodes>(result.exit_code);
  }
};

int main(int argc, const char** argv)
{
  TOPPExecutePipeline tool;
  return tool.mainWithUtf8Arguments(argc, argv);
}
/// @endcond
