// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/TOPPASPipeline.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <map>
#include <string>
#include <vector>

namespace OpenMS
{
  class TOPPASWorkflow;

  /**
    @brief Qt-free executor for TOPPAS workflows.

    PipelineExecutor turns a parsed @em .toppas workflow (OpenMS::TOPPASWorkflow) into an executable
    pipeline and runs it without the GUI. It builds on the Qt-free data-flow routing of
    OpenMS::TOPPASPipeline and replaces the historical Qt machinery (QProcess + signal/slot reactor +
    event loop) with a straightforward synchronous topological walk that launches each TOPP tool via
    the (blocking) OpenMS::ExternalProcess.

    It is the engine behind a Qt-free @c ExecutePipeline. Input file lists, mergers, splitters and
    output recycling are routed by the base class; this class only adds the actual tool invocation
    (INI generation, command line, subprocess) and the copying of results to the output directory.

    @ingroup Applications
  */
  class OPENMS_DLLAPI PipelineExecutor : public TOPPASPipeline
  {
  public:
    /// Runtime options.
    struct Options
    {
      std::string out_dir;  ///< directory the output nodes write their files to
      std::string tmp_dir;  ///< scratch directory for per-tool INI files and intermediate results
      int num_jobs = 1;     ///< maximum parallel jobs (currently executed serially)
    };

    /// Result of execute(); the numeric values double as process exit codes for ExecutePipeline.
    enum Result
    {
      OK = 0,              ///< the whole pipeline ran successfully
      TOOL_FAILED = 1,     ///< a tool returned an error / crashed
      INVALID = 2,         ///< the workflow graph is invalid (cycle, inconsistent rounds, ...)
      TOOL_NOT_FOUND = 14  ///< a required TOPP tool executable was not found (matches EXTERNAL_PROGRAM_NOTFOUND)
    };

    /// Build the executable graph from a parsed workflow. @p toppas_file is used to resolve relative
    /// input paths and to name the per-tool output directories.
    void fromWorkflow(const TOPPASWorkflow& wf, const std::string& toppas_file);

    /// Run the pipeline. @p input_overrides optionally replaces the files of input nodes, keyed by the
    /// input node's topological number (this is the .trf resource-file convention).
    Result execute(const Options& opts, const std::map<std::string, std::vector<std::string>>& input_overrides,
                   std::string& error_msg);

    /// A tool's input or output file parameter (resolved from its Param), used for edge<->index mapping.
    struct IOParam
    {
      std::string name;     ///< parameter name (matches the edge's source_out_param / target_in_param)
      bool is_list = false; ///< parameter holds a list of files
      std::string ext;      ///< first valid file extension (without the leading dot), or "" if unrestricted
    };

    /// Collect a tool's input (@p input == true) or output file parameters from its @p param, ordered
    /// identically to the historical TOPPASToolVertex (so the index matches what .toppas edges reference).
    static std::vector<IOParam> collectIO(const Param& param, bool input);

    /// Parse a TOPPAS resource file (.trf) into input-file overrides keyed by the input node's
    /// topological number (the value of each '<node>:url_list', with the 'file:'/'file://' scheme stripped).
    static std::map<std::string, std::vector<std::string>> loadResourceFile(const std::string& filename);

  private:
    /// Per-node description (TOOL nodes have a populated tool entry; output nodes carry their folder name).
    struct ToolDesc
    {
      bool is_tool = false;
      std::string name;
      std::string type;
      Param param;                  ///< the tool's parameters (bare, as stored in the .toppas)
      std::vector<IOParam> in;      ///< input-file params, sorted -> index == target_in_param
      std::vector<IOParam> out;     ///< output-file params, sorted -> index == source_out_param
      std::string output_folder_name; ///< OUTPUT / OUTPUT_FOLDER nodes only: custom output folder name (or empty)
    };

    /// One scheduled subprocess: running @p exe with @p args produces one round of one tool node.
    struct Task
    {
      int node = -1;             ///< the tool node this task belongs to
      std::string exe;           ///< tool executable
      std::vector<std::string> args; ///< command line (-type ... -ini <inifile>)
    };

    /// Result of a finished Task (as reported back by a worker thread).
    struct TaskResult
    {
      int node = -1;
      Result status = OK;
      std::string msg;
    };

    /// Thrown internally to abort execution with a specific result code.
    struct ExecError
    {
      Result result;
      std::string msg;
    };

    /// Set @p name in @p p to @p files (as a list, or as a single value).
    static void setFileParam_(Param& p, const std::string& name, const std::vector<std::string>& files, bool is_list);

    /// Prepare a tool node for execution: build the per-round output file names (filling @c output_files)
    /// and return one Task per round. When @p dry_run is true, only the output names are computed (no
    /// working directory is created and no INI is written) and an empty task list is returned.
    std::vector<Task> prepareToolTasks_(int node_index, const RoundPackages& inputs, bool dry_run = false);

    /// Copy a finished output node's incoming files into its output directory.
    void runOutputNode_(int node_index);

    /// Validate the input nodes the way the GUI's TOPPASScene::sanityCheck_ does for non-interactive runs:
    /// there must be at least one input node, and every connected input node must have a non-empty list of
    /// existing, non-duplicate files. Returns OK if the pipeline passes, INVALID otherwise.
    Result sanityCheck_(std::string& error_msg) const;

    /// Single-threaded, topological dry-run pass that routes the whole data flow (input seeding, mergers,
    /// splitters, recycling, tool output-name generation) WITHOUT launching any tool or copying any file.
    /// Mirrors TOPPAS's "dry run" so graph/routing errors are caught before any real tool is started.
    /// Returns OK on success, or the failing result code (with @p error_msg set).
    Result validateDryRun_(std::string& error_msg);

    /// Reset the per-node data-flow state (output_files / round_total / finished) between the dry and the
    /// real run. The INPUT nodes' seed @c input_files are preserved.
    void resetDataFlow_();

    /// Display name of a node mirroring TOPPASVertex::getName() (used to knit output directory names).
    std::string vertexName_(int node_index) const;

    /// Output sub-directory of an OUTPUT / OUTPUT_FOLDER node relative to the output root, mirroring
    /// TOPPASOutputVertex::getOutputDir: "TOPPAS_out/<folder>" if a custom folder name is set, else
    /// "TOPPAS_out/<NNN>-<source-vertex>-<source-out-param-without-colons>".
    std::string outputSubDir_(int node_index) const;

    std::vector<ToolDesc> tools_; ///< indexed by node index
    std::string toppas_stem_;     ///< basename of the .toppas (for output-dir names)
    Options opts_;
  };

} // namespace OpenMS
