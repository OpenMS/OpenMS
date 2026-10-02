// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// $Maintainer: Timo Sachsenberg $

#pragma once

#include <OpenMS/APPLICATIONS/OpenMS_CLIConfig.h>
#include <OpenMS/APPLICATIONS/PIPELINE/PipelineGraph.h>
#include <atomic>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace OpenMS
{
/**
  @brief Executes a TOPPAS workflow without a GUI or an event loop.

  run() blocks its calling thread. cancel() may be called from another thread.
  Events are serialized on the calling thread; they must be marshalled onto a
  GUI thread by the caller. An instance runs one workflow at a time.
*/
class OPENMS_CLI_DLLAPI PipelineExecutor
{
public:
  using Files = std::vector<std::string>;
  using FileBundle = std::map<std::string, Files>;
  using Rounds = std::vector<FileBundle>;

  enum class State
  {
    PENDING,
    RUNNING,
    SUCCEEDED,
    FAILED,
    CANCELLED,
    BLOCKED
  };

  struct NodeResult
  {
    State state {State::PENDING};
    Rounds outputs;
    Size completed_rounds {0};
    Size total_rounds {0};
    std::string error_message;
  };

  struct Result
  {
    int exit_code {0};
    std::string error_message;
    std::map<Size, NodeResult> nodes;
  };

  struct Options
  {
    std::string output_directory;
    /// Parent of an exclusively owned temporary run directory.
    std::string temp_directory;
    Size num_jobs {1};
    /// Required when retaining upstream results for a subsequent rerun.
    bool keep_temporary_files {false};
  };

  struct Event
  {
    enum class Type
    {
      LOG,
      NODE_SCHEDULED,
      NODE_STARTED,
      ROUND_COMPLETED,
      NODE_FINISHED,
      OUTPUT_WRITTEN,
      NODE_FAILED
    };
    Type type {Type::LOG};
    Size node_id {0};
    Size round {0};
    Size total {0};
    std::string text;
    int exit_code {0};
    /// Snapshot included with node/round events for presentation adapters.
    NodeResult result {};
  };

  using EventCallback = std::function<void(const Event&)>;

  PipelineExecutor() = default;
  PipelineExecutor(const PipelineExecutor&) = delete;
  PipelineExecutor& operator=(const PipelineExecutor&) = delete;

  /**
    @brief Run an immutable snapshot, optionally rerunning a downstream subgraph.

    A rerun requires completed upstream results in @p previous_result. Their
    files must still exist; this is in-memory rerun, not checkpoint recovery.
    Errors are returned as TOPP exit codes with contextual diagnostics.
  */
  Result run(const PipelineGraph& graph,
             const Options& options,
             EventCallback callback = {},
             const Result* previous_result = nullptr,
             std::optional<Size> start_node = std::nullopt);

  /// Cancel the current (or next) run. The request is consumed when run() returns.
  void cancel() noexcept;

private:
  std::atomic_bool cancelled_ {false};
  std::atomic_bool running_ {false};
};
} // namespace OpenMS
