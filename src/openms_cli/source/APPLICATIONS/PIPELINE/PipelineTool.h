// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineExecutor.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <filesystem>

namespace OpenMS::Internal
{
/// Detect aliases using the destination filesystem's own filename rules.
/// The exclusively owned probe directory contains no workflow data.
class PipelinePathClaims
{
public:
  explicit PipelinePathClaims(const std::filesystem::path& root);
  ~PipelinePathClaims();
  PipelinePathClaims(const PipelinePathClaims&) = delete;
  PipelinePathClaims& operator=(const PipelinePathClaims&) = delete;
  bool claim(const std::filesystem::path& path);

private:
  std::filesystem::path root_;
  std::filesystem::path probe_;
};

/// Private preparation and filename handling shared by executor jobs.
struct PipelineTool
{
  using FileBundle = PipelineExecutor::FileBundle;
  using Rounds = PipelineExecutor::Rounds;

  struct Port
  {
    enum class Kind
    {
      FILE,
      LIST,
      DIRECTORY
    };
    std::string name;
    Kind kind {Kind::FILE};
    std::vector<std::string> valid_types;
  };

  struct Descriptor
  {
    std::string executable;
    Param parameters;
    std::vector<Port> inputs;
    std::vector<Port> outputs;
  };

  struct Invocation
  {
    std::string executable;
    std::vector<std::string> arguments;
    std::string working_directory;
  };

  /// Preserve numeric process outcomes while unwinding preparation failures.
  struct Failure : Exception::BaseException
  {
    Failure(int code, const std::string& message): Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, code, message)
    {
    }
    Failure(const char* file, int line, const char* function, int code, const std::string& message):
        Exception::BaseException(file, line, function, "PipelineToolFailure", message),
        exit_code(code)
    {
    }
    int exit_code;
  };

  static Descriptor discover(const PipelineGraph::Node& node,
                             const std::string& temp_directory,
                             const std::atomic_bool& cancelled,
                             std::function<void(const std::string&)> log);

  static Rounds planOutputs(const PipelineGraph& graph,
                            const PipelineGraph::Node& node,
                            const Rounds& inputs,
                            const Descriptor& descriptor,
                            const std::string& run_temp);

  static Invocation invocation(const PipelineGraph& graph,
                               const PipelineGraph::Node& node,
                               const Descriptor& descriptor,
                               const FileBundle& inputs,
                               const FileBundle& outputs,
                               Size round,
                               const std::string& run_temp);

  /// Detect content formats after all rounds finish and rename without overwriting files.
  static void finalizeOutputs(Rounds& outputs);
};
} // namespace OpenMS::Internal
