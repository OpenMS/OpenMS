// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <thread>

using namespace OpenMS;
namespace fs = std::filesystem;

/** A real TOPP subprocess used only by PipelineExecutor_test. */
class PipelineTestTool final : public TOPPBase
{
public:
  PipelineTestTool(): TOPPBase("PipelineTestTool", "Test subprocess for the workflow executor.", {}, false)
  {
  }

private:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", {}, "List input.", false);
    registerInputFile_("single_in", "<file>", "", "Single file input.", false);
    registerTOPPSubsection_("ETool", "Nested external tool parameters.");
    registerInputFileList_("ETool:in", "<files>", {}, "Nested input requiring INI binding.", false);
    registerOutputFile_("out", "<file>", "", "Single file output.", false);
    registerOutputFile_("untyped_out", "<file>", "", "Output whose format must be detected from content.", false);
    registerOutputFileList_("out_list", "<files>", {}, "List output.", false);
    registerOutputDir_("out_dir", "<directory>", "", "Directory output.", false);
    setValidFormats_("in", {"fasta"});
    setValidFormats_("single_in", {"fasta"});
    setValidFormats_("ETool:in", {"fasta"});
    setValidFormats_("out", {"fasta"});
    setValidFormats_("out_list", {"fasta"});
    registerStringOption_("mode", "<mode>", "copy", "Whether to copy inputs or fail.", false);
    setValidStrings_("mode", {"copy", "fail", "skip_output"});
    registerIntOption_("exit_code", "<code>", 9, "Exit code in failure mode.", false);
    setMinInt_("exit_code", 1);
    setMaxInt_("exit_code", 127);
    registerIntOption_("delay_ms", "<milliseconds>", 0, "Delay before writing outputs.", false);
    setMinInt_("delay_ms", 0);
    registerStringOption_("trace_directory", "<directory>", "", "Optional per-process timing traces.", false);
  }

  static long long now_()
  { return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count(); }

  static bool write_(const fs::path& path, const std::string& contents)
  {
    fs::create_directories(path.parent_path());
    std::ofstream stream(path, std::ios::binary);
    stream << contents;
    return stream.good();
  }

  ExitCodes main_(int, const char**) override
  {
    const auto started = now_();
    const auto trace_directory = getStringOption_("trace_directory");
    fs::path trace;
    if (! trace_directory.empty())
    {
      fs::create_directories(to_path(trace_directory));
      trace = to_path(trace_directory) / (File::getUniqueName() + ".trace");
      // Each child owns one file; no concurrent appends or platform-specific locks.
      if (! write_(fs::path(trace).concat(".started"), std::to_string(started))) return CANNOT_WRITE_OUTPUT_FILE;
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(getIntOption_("delay_ms")));
    if (getStringOption_("mode") == "fail") return static_cast<ExitCodes>(getIntOption_("exit_code"));
    if (getStringOption_("mode") == "skip_output") return EXECUTION_OK;

    auto inputs = getStringList_("in");
    const auto single = getStringOption_("single_in");
    if (! single.empty()) inputs.push_back(single);
    const auto nested = getStringList_("ETool:in");
    inputs.insert(inputs.end(), nested.begin(), nested.end());
    std::string content;
    for (const auto& input : inputs)
    {
      std::ifstream stream(to_path(input), std::ios::binary);
      if (! stream) return INPUT_FILE_NOT_READABLE;
      content.append(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
    }
    if (content.empty()) content = ">generated\nPEPTIDE\n";

    const auto output = getStringOption_("out");
    if (! output.empty() && ! write_(to_path(output), content)) return CANNOT_WRITE_OUTPUT_FILE;
    const auto untyped = getStringOption_("untyped_out");
    if (! untyped.empty() && ! write_(to_path(untyped), content)) return CANNOT_WRITE_OUTPUT_FILE;
    for (const auto& name : getStringList_("out_list"))
    {
      if (! write_(to_path(name), content)) return CANNOT_WRITE_OUTPUT_FILE;
    }
    const auto directory = getOutputDirOption("out_dir");
    if (! directory.empty())
    {
      const auto path = to_path(directory);
      if (! write_(path / "direct.fasta", content) || ! write_(path / "nested" / "nested.fasta", content)) return CANNOT_WRITE_OUTPUT_FILE;
    }
    if (! trace.empty() && ! write_(trace, std::to_string(started) + " " + std::to_string(now_()) + "\n")) { return CANNOT_WRITE_OUTPUT_FILE; }
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  PipelineTestTool tool;
  return tool.main(argc, argv);
}
