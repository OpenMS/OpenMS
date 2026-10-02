// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// $Maintainer: Timo Sachsenberg $

#include "PipelineTool.h"

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineExecutor.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileNameUtils.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <OpenMS/SYSTEM/SystemSettings.h>
#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <fstream>
#include <future>
#include <mutex>
#include <set>
#include <utility>

#ifdef _WIN32
  #include <Windows.h>
#endif

namespace OpenMS
{
namespace
{
  namespace fs = std::filesystem;
  using Engine = PipelineExecutor;
  using Graph = PipelineGraph;
  using Kind = Graph::Kind;
  using Tool = Internal::PipelineTool;
  using Failure = Tool::Failure;
  using Descriptors = std::map<Size, Tool::Descriptor>;

  // Keep an explicit cancellation outcome distinct from cancelling peer jobs
  // while unwinding a genuine tool failure.
  struct Cancelled : Failure
  {
    explicit Cancelled(const std::string& message = "Pipeline cancelled."): Failure(TOPPBase::EXTERNAL_PROGRAM_ERROR, message)
    {
    }
  };

  std::string pathString(const fs::path& path)
  {
    const auto text = path.u8string();
    return {text.begin(), text.end()};
  }

  std::string context(const Graph::Node& node)
  { return "Node #" + std::to_string(node.topo_number) + (node.tool_name.empty() ? std::string {} : " (" + node.tool_name + ")"); }

  std::vector<const Graph::Edge*> incoming(const Graph& graph, Size id)
  {
    std::vector<const Graph::Edge*> edges;
    for (const auto& edge : graph.edges)
    {
      if (edge.target == id) edges.push_back(&edge);
    }
    return edges;
  }

  const Tool::Port& port(const std::vector<Tool::Port>& ports, const std::string& name, const Graph::Node& node)
  {
    const auto it = std::find_if(ports.begin(), ports.end(), [&](const auto& p) { return p.name == name; });
    if (it == ports.end())
    {
      throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA,
                    context(node) + ": unknown file parameter '" + name + "'. Refresh the workflow's tool parameters.");
    }
    return *it;
  }

  void resolveLegacyPorts(Graph& graph, const Descriptors& descriptors)
  {
    if (! graph.legacy_port_indices) return;
    auto resolve = [&](Size id, std::string& name, bool output) {
      const auto& node = graph.node(id);
      if (node.kind != Kind::TOOL || ! descriptors.count(id)) return;
      const auto& ports = output ? descriptors.at(id).outputs : descriptors.at(id).inputs;
      Size consumed = 0;
      unsigned long long index = 0;
      try
      {
        index = std::stoull(name, &consumed);
      }
      catch (const std::exception&)
      {
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": invalid legacy port index '" + name + "'.");
      }
      if (name.empty() || name.front() == '-' || consumed != name.size() || index >= ports.size())
      {
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": legacy port index is out of bounds.");
      }
      name = ports[index].name;
    };
    for (auto& edge : graph.edges)
    {
      resolve(edge.source, edge.source_port, true);
      resolve(edge.target, edge.target_port, false);
    }
    graph.legacy_port_indices = false;
  }

  void checkSourceFormat(const Graph& graph, const Descriptors& descriptors, const Graph::Edge& edge, const Tool::Port& target)
  {
    const auto& source = graph.node(edge.source);
    if (source.kind == Kind::TOOL)
    {
      const auto& output = port(descriptors.at(source.id).outputs, edge.source_port, source);
      if (output.kind == Tool::Port::Kind::DIRECTORY)
      {
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(source) + ": a directory output needs an output-folder node.");
      }
      if (! output.valid_types.empty() && ! target.valid_types.empty())
      {
        for (const auto& a : output.valid_types)
        {
          for (const auto& b : target.valid_types)
          {
            if (FileTypes::sameFormat(a, b)) return;
          }
        }
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA,
                      context(source) + ": incompatible formats for output '" + edge.source_port + "' and input '" + target.name + "'.");
      }
    }
    else if (source.kind == Kind::INPUT && ! target.valid_types.empty())
    {
      for (const auto& filename : source.files)
      {
        const auto type = FileNameUtils::getTypeByFileName(filename);
        const auto compression = FileNameUtils::compressionType(filename);
        bool matches = false;
        if (type != FileTypes::UNKNOWN && (compression == FileTypes::UNKNOWN || FileTypes::supportsCompressedReading(type, compression)))
        {
          for (const auto& format : target.valid_types)
          {
            if (FileTypes::nameToType(format) == type) matches = true;
          }
        }
        else if (type == FileTypes::UNKNOWN)
        {
          const auto dot = filename.rfind('.');
          if (dot != std::string::npos)
          {
            for (const auto& format : target.valid_types)
            {
              if (FileTypes::sameFormat(format, filename.substr(dot + 1))) matches = true;
            }
          }
        }
        if (! matches)
          throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, "Input '" + filename + "' is incompatible with parameter '" + target.name + "'.");
      }
    }
    else if (source.kind == Kind::MERGER || source.kind == Kind::SPLITTER)
    {
      for (const auto* upstream : incoming(graph, source.id))
        checkSourceFormat(graph, descriptors, *upstream, target);
    }
  }

  void validatePorts(const Graph& graph, const Descriptors& descriptors, const std::set<Size>& reachable)
  {
    for (const auto& edge : graph.edges)
    {
      if (! reachable.count(edge.target)) continue;
      const auto& source = graph.node(edge.source);
      const auto& target = graph.node(edge.target);
      if (source.kind == Kind::TOOL)
      {
        const auto& output = port(descriptors.at(source.id).outputs, edge.source_port, source);
        if ((output.kind == Tool::Port::Kind::DIRECTORY) != (target.kind == Kind::OUTPUT_DIRECTORY)
            && (output.kind == Tool::Port::Kind::DIRECTORY || target.kind == Kind::OUTPUT_DIRECTORY))
        {
          throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(target) + ": file/directory output binding does not match.");
        }
      }
      else if (target.kind == Kind::OUTPUT_DIRECTORY)
      {
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(target) + ": output folders require a tool directory port.");
      }
      if (target.kind == Kind::TOOL)
      {
        const auto& input = port(descriptors.at(target.id).inputs, edge.target_port, target);
        checkSourceFormat(graph, descriptors, edge, input);
      }
    }
  }

  const Engine::Files& filesAt(const Engine::Rounds& rounds, Size round, const std::string& name)
  {
    const auto& bundle = rounds.at(round);
    const auto it = bundle.find(name);
    if (it == bundle.end()) throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, "Upstream node did not provide output port '" + name + "'.");
    return it->second;
  }

  /// Align rounds first, then preserve the historical reverse edge order of mergers.
  Engine::Rounds inputsFor(const Graph& graph, const Graph::Node& node, const Engine::Result& result)
  {
    auto edges = incoming(graph, node.id);
    if (edges.empty()) throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": no upstream inputs.");
    // Output nodes publish every source round. Recycling only aligns the
    // inputs of processing/structural nodes; it must not suppress publication.
    if (node.kind == Kind::OUTPUT || node.kind == Kind::OUTPUT_DIRECTORY)
    {
      const auto& edge = *edges.front();
      const auto& previous = result.nodes.at(edge.source);
      if (previous.state != Engine::State::SUCCEEDED)
      {
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": upstream outputs are not complete.");
      }
      Engine::Rounds inputs(previous.outputs.size());
      for (Size round = 0; round < previous.outputs.size(); ++round)
      {
        inputs[round][""] = filesAt(previous.outputs, round, edge.source_port);
      }
      return inputs;
    }
    std::optional<Size> count;
    for (const auto* edge : edges)
    {
      const auto& previous = result.nodes.at(edge->source);
      if (previous.state != Engine::State::SUCCEEDED)
      {
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": required upstream node has no completed results.");
      }
      if (! graph.node(edge->source).recycle_output)
      {
        if (count && *count != previous.outputs.size())
          throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": incoming round counts differ; enable recycling where intended.");
        count = previous.outputs.size();
      }
    }
    if (! count || *count == 0)
      throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": at least one non-recycled input with positive round count is required.");
    for (const auto* edge : edges)
    {
      const auto rounds = result.nodes.at(edge->source).outputs.size();
      if (rounds == 0 || *count % rounds != 0)
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": recycled round count must be a positive divisor of the total.");
    }
    if (node.kind == Kind::MERGER) std::reverse(edges.begin(), edges.end());
    Engine::Rounds inputs(*count);
    for (Size round = 0; round < *count; ++round)
    {
      for (const auto* edge : edges)
      {
        const auto& upstream = result.nodes.at(edge->source).outputs;
        const auto& files = filesAt(upstream, round % upstream.size(), edge->source_port);
        auto& target = inputs[round][node.kind == Kind::TOOL ? edge->target_port : std::string {}];
        target.insert(target.end(), files.begin(), files.end());
      }
    }
    return inputs;
  }

  Engine::Rounds transform(const Graph::Node& node, const Engine::Rounds& inputs)
  {
    if (node.kind == Kind::MERGER && node.round_based) return inputs;
    Engine::Rounds outputs;
    if (node.kind == Kind::MERGER) outputs.resize(1);
    for (const auto& bundle : inputs)
    {
      for (const auto& filename : bundle.at(""))
      {
        if (node.kind == Kind::SPLITTER) outputs.push_back({{"", {filename}}});
        else
          outputs[0][""].push_back(filename);
      }
    }
    return outputs;
  }

  std::string displayName(const Graph::Node& node)
  {
    if (node.kind == Kind::TOOL) return node.tool_name;
    if (node.kind == Kind::MERGER) return "MergerVertex";
    if (node.kind == Kind::SPLITTER) return "SplitterVertex";
    return "InputVertex";
  }

  fs::path outputDirectory(const Graph& graph, const Graph::Node& node, const std::string& root)
  {
    const auto base = fs::canonical(to_path(root)) / "TOPPAS_out";
    auto contained = [&](const fs::path& directory) {
      const auto resolved = fs::weakly_canonical(directory);
      const auto relative = resolved.lexically_relative(base);
      if (relative.empty() || relative.has_root_path() || *relative.begin() == "..")
      {
        throw Failure(TOPPBase::ILLEGAL_PARAMETERS, context(node) + ": output folder resolves outside TOPPAS_out.");
      }
      return resolved;
    };
    if (! node.output_folder.empty())
    {
      const fs::path custom = to_path(node.output_folder);
      if (custom.has_root_path() || std::any_of(custom.begin(), custom.end(), [](const auto& part) { return part == ".."; }))
      {
        throw Failure(TOPPBase::ILLEGAL_PARAMETERS, context(node) + ": output folder must stay within TOPPAS_out.");
      }
      return contained(base / custom);
    }
    const auto* edge = incoming(graph, node.id).front();
    std::string number = std::to_string(node.topo_number);
    if (number.size() < 3) number.insert(0, 3 - number.size(), '0');
    std::string name = edge->source_port;
    name.erase(std::remove(name.begin(), name.end(), ':'), name.end());
    return contained(base / to_path(number + "-" + displayName(graph.node(edge->source)) + "-" + name));
  }

  void copyOutput(const fs::path& source, const fs::path& target)
  {
    fs::create_directories(target.parent_path());
    if (fs::exists(target) && fs::equivalent(source, target)) return;
    const fs::path staged = target.parent_path() / (Tool::temporaryName(".cp-") + ".part");
    try
    {
      fs::copy_file(source, staged);
      // Complete the copy before replacing an existing destination.
#ifdef _WIN32
      if (! MoveFileExW(staged.c_str(), target.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
      {
        throw fs::filesystem_error("Could not replace workflow output", staged, target, std::error_code(GetLastError(), std::system_category()));
      }
#else
      fs::rename(staged, target);
#endif
    }
    catch (...)
    {
      std::error_code ec;
      fs::remove(staged, ec);
      throw;
    }
  }

  Engine::Rounds publish(const Graph& graph,
                         const Graph::Node& node,
                         const Engine::Rounds& inputs,
                         const Descriptors& descriptors,
                         const std::string& root,
                         bool dry,
                         const std::atomic_bool& cancelled,
                         const Engine::EventCallback& emit,
                         Internal::PipelinePathClaims& claims)
  {
    const auto directory = outputDirectory(graph, node, root);
    const auto* edge = incoming(graph, node.id).front();
    Engine::Rounds outputs(inputs.size());
    struct Copy
    {
      fs::path source;
      fs::path target;
      Size round;
    };
    std::vector<Copy> copies;
    for (Size round = 0; round < inputs.size(); ++round)
    {
      for (const auto& filename : inputs[round].at(""))
      {
        if (cancelled.load()) throw Cancelled();
        const fs::path source = to_path(filename);
        if (node.kind == Kind::OUTPUT_DIRECTORY)
        {
          fs::path destination = directory;
          if (inputs.size() > 1) destination /= source.filename();
          outputs[round][""].push_back(pathString(destination));
          if (dry) continue;
          if (! fs::is_directory(source)) throw Failure(TOPPBase::CANNOT_WRITE_OUTPUT_FILE, "Missing output directory '" + filename + "'.");
          std::vector<fs::path> files;
          for (const auto& entry : fs::directory_iterator(source))
          {
            if (entry.is_regular_file()) files.push_back(entry.path());
          }
          if (files.empty()) throw Failure(TOPPBase::CANNOT_WRITE_OUTPUT_FILE, "Output directory '" + filename + "' contains no files.");
          std::sort(files.begin(), files.end());
          for (const auto& file : files)
          {
            if (cancelled.load()) throw Cancelled();
            const auto target = destination / file.filename();
            copies.push_back({file, target, round});
          }
        }
        else
        {
          auto type = FileTypes::UNKNOWN;
          if (graph.node(edge->source).kind == Kind::TOOL)
          {
            const auto& descriptor = descriptors.at(edge->source);
            const auto& output = port(descriptor.outputs, edge->source_port, graph.node(edge->source));
            if (output.valid_types.size() == 1) type = FileTypes::nameToType(output.valid_types.front());
            else if (! dry)
              type = FileHandler::getTypeByContent(filename);
            if (type == FileTypes::UNKNOWN && descriptor.parameters.exists(output.name + "_type"))
            {
              type = FileTypes::nameToType(descriptor.parameters.getValue(output.name + "_type").toString());
            }
          }
          std::string basename = pathString(source.filename());
          const auto suffix = basename.rfind("_tmp");
          if (suffix != std::string::npos && suffix + 4 < basename.size()
              && std::all_of(basename.begin() + suffix + 4, basename.end(), [](char c) { return c >= '0' && c <= '9'; }))
            basename.erase(suffix);
          std::string output_name = basename;
          if (type != FileTypes::UNKNOWN)
          {
            output_name = FileHandler::swapExtension(basename, type);
            if (const auto compression = FileNameUtils::compressionType(basename); compression != FileTypes::UNKNOWN)
            {
              output_name += "." + FileTypes::typeToName(compression);
            }
          }
          const auto target = directory / to_path(output_name);
          outputs[round][""].push_back(pathString(target));
          // Untyped outputs can change suffix after execution; reserve those then.
          if (! dry || type != FileTypes::UNKNOWN) { copies.push_back({source, target, round}); }
          if (! dry)
          {
            if (! fs::is_regular_file(source)) throw Failure(TOPPBase::CANNOT_WRITE_OUTPUT_FILE, "Missing output file '" + filename + "'.");
          }
        }
      }
    }
    // Reserve all destinations before writing any file from this output node.
    for (const auto& copy : copies)
    {
      if (! claims.claim(copy.target))
      {
        throw Failure(TOPPBase::CANNOT_WRITE_OUTPUT_FILE, "Workflow outputs collide at '" + pathString(copy.target) + "'.");
      }
    }
    if (! dry)
    {
      for (const auto& copy : copies)
      {
        if (cancelled.load()) { throw Cancelled(); }
        copyOutput(copy.source, copy.target);
        emit({Engine::Event::Type::OUTPUT_WRITTEN, node.id, copy.round, inputs.size(), pathString(copy.target)});
      }
    }
    return outputs;
  }

  void checkInputCardinality(const Graph::Node& node, const Tool::Descriptor& descriptor, const Engine::Rounds& inputs)
  {
    for (const auto& bundle : inputs)
    {
      for (const auto& [name, files] : bundle)
      {
        const auto& input = port(descriptor.inputs, name, node);
        if (files.empty() || (input.kind == Tool::Port::Kind::FILE && files.size() != 1))
        {
          throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA,
                        context(node) + ": input '" + name + "' has invalid file cardinality; use a splitter for scalar inputs.");
        }
      }
    }
  }

  struct TemporaryDirectory
  {
    fs::path path;
    bool keep;
    ~TemporaryDirectory()
    {
      if (! keep && ! path.empty())
      {
        std::error_code ec;
        fs::remove_all(path, ec);
      }
    }
  };
} // namespace

void PipelineExecutor::cancel() noexcept
{ cancelled_.store(true); }

PipelineExecutor::Result PipelineExecutor::run(const PipelineGraph& definition,
                                               const Options& options,
                                               EventCallback callback,
                                               const Result* previous_result,
                                               std::optional<Size> start_node)
{
  Result result;
  std::ofstream logfile;
  if (running_.exchange(true)) return {TOPPBase::ILLEGAL_PARAMETERS, "This executor is already running.", {}};
  struct RunningGuard
  {
    std::atomic_bool& running;
    std::atomic_bool& cancelled;
    ~RunningGuard()
    {
      cancelled.store(false);
      running.store(false);
    }
  } running_guard {running_, cancelled_};
  Size current_node = 0;
  bool cancellation_outcome = false;
  try
  {
    if (cancelled_.load()) throw Cancelled();
    if (options.num_jobs == 0) throw Failure(TOPPBase::ILLEGAL_PARAMETERS, "num_jobs must be positive.");
    Graph graph = definition;
    graph.validate();
    graph.assignTopologicalNumbers();
    const auto order = graph.topologicalOrder();
    if (std::none_of(graph.nodes.begin(), graph.nodes.end(), [](const auto& n) { return n.kind == Kind::INPUT; }))
    {
      throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, "The pipeline contains no input nodes.");
    }
    const fs::path out = fs::absolute(to_path(options.output_directory.empty() ? "." : options.output_directory));
    fs::create_directories(out);
    const fs::path parent = options.temp_directory.empty() ? to_path(SystemSettings::getTempDirectory()) : to_path(options.temp_directory);
    fs::create_directories(parent);
    TemporaryDirectory temporary {{}, options.keep_temporary_files};
    for (Size attempt = 0; attempt < 10 && temporary.path.empty(); ++attempt)
    {
      const auto candidate = parent / Tool::temporaryName("omp-");
      if (fs::create_directory(candidate)) temporary.path = fs::absolute(candidate);
    }
    if (temporary.path.empty()) throw Failure(TOPPBase::CANNOT_WRITE_OUTPUT_FILE, "Could not create an exclusive temporary run directory.");
    logfile.open(out / "TOPPAS.log", std::ios::trunc);
    if (! logfile) throw Failure(TOPPBase::CANNOT_WRITE_OUTPUT_FILE, "Cannot write TOPPAS.log in '" + pathString(out) + "'.");
    auto emit = [&](const Event& event) {
      if (! event.text.empty())
      {
        logfile << event.text;
        if (event.text.back() != '\n') logfile << '\n';
        logfile.flush();
      }
      if (callback) callback(event);
    };

    std::set<Size> affected;
    if (start_node)
    {
      graph.node(*start_node);
      if (! previous_result) throw Failure(TOPPBase::ILLEGAL_PARAMETERS, "Rerun requires previous upstream results.");
      affected.insert(*start_node);
      for (const auto id : order)
      {
        if (! affected.count(id)) continue;
        for (const auto& edge : graph.edges)
          if (edge.source == id) affected.insert(edge.target);
      }
    }
    else
    {
      for (const auto id : order)
        affected.insert(id);
    }
    for (const auto& node : graph.nodes)
    {
      result.nodes[node.id] = NodeResult {};
      if (! affected.count(node.id) && previous_result)
      {
        const auto found = previous_result->nodes.find(node.id);
        if (found != previous_result->nodes.end()) result.nodes[node.id] = found->second;
      }
    }

    // Ignore wholly disconnected editing leftovers, but never report success
    // for a live join whose required sibling branch cannot produce results.
    std::set<Size> reachable;
    for (const auto id : order)
    {
      current_node = id;
      const auto& node = graph.node(id);
      const auto edges = incoming(graph, id);
      const auto available = std::count_if(edges.begin(), edges.end(), [&](const auto* edge) { return reachable.count(edge->source); });
      if (node.kind == Kind::INPUT || (! edges.empty() && static_cast<Size>(available) == edges.size())) reachable.insert(id);
      else if (available != 0)
      {
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": a required branch is disconnected from all inputs.");
      }
      else
        result.nodes[id].state = State::BLOCKED;
    }

    // Tool introspection is explicit and uses the same cancellable backend as jobs.
    Descriptors descriptors;
    for (const auto& node : graph.nodes)
    {
      current_node = node.id;
      if (node.kind == Kind::TOOL && reachable.count(node.id))
      {
        try
        {
          descriptors.emplace(node.id, Tool::discover(node, pathString(temporary.path), cancelled_,
                                                      [&](const auto& text) { emit({Event::Type::LOG, node.id, 0, 0, text}); }));
        }
        catch (const Failure& error)
        {
          // There are no peer jobs yet, so this flag can only come from cancel().
          if (cancelled_.load()) throw Cancelled(error.what());
          throw;
        }
      }
    }
    resolveLegacyPorts(graph, descriptors);
    graph.validate(); // Legacy indices may have resolved to duplicate named bindings.
    validatePorts(graph, descriptors, reachable);

    // Validate retained boundary data before starting any processing jobs.
    // Reruns must not partially execute when a cached upstream file was removed.
    for (const auto& edge : graph.edges)
    {
      if (affected.count(edge.source) || ! affected.count(edge.target) || ! reachable.count(edge.target)) continue;
      current_node = edge.source;
      const auto& source = graph.node(edge.source);
      const auto& retained = result.nodes.at(edge.source);
      if (retained.state != State::SUCCEEDED)
      {
        throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(source) + ": rerun requires completed upstream results.");
      }
      const bool directory
        = source.kind == Kind::TOOL && port(descriptors.at(source.id).outputs, edge.source_port, source).kind == Tool::Port::Kind::DIRECTORY;
      for (Size round = 0; round < retained.outputs.size(); ++round)
      {
        for (const auto& file : filesAt(retained.outputs, round, edge.source_port))
        {
          if (! fs::exists(to_path(file)))
          {
            throw Failure(TOPPBase::INPUT_FILE_NOT_FOUND, context(source) + ": retained output '" + file + "' does not exist; rerun its producer.");
          }
          if (directory != fs::is_directory(to_path(file)))
          {
            throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(source) + ": retained output '" + file + "' has the wrong file type.");
          }
        }
      }
    }

    auto prepareInput = [&](const Graph::Node& node) {
      NodeResult input;
      input.state = State::SUCCEEDED;
      const bool connected = std::any_of(graph.edges.begin(), graph.edges.end(), [&](const auto& e) { return e.source == node.id; });
      std::set<fs::path> unique;
      if (connected && node.files.empty()) throw Failure(TOPPBase::INPUT_FILE_EMPTY, context(node) + ": empty input file list.");
      for (const auto& file : node.files)
      {
        if (connected)
        {
          if (! fs::exists(to_path(file))) throw Failure(TOPPBase::INPUT_FILE_NOT_FOUND, context(node) + ": input '" + file + "' does not exist.");
          if (! unique.insert(fs::canonical(to_path(file))).second)
            throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, context(node) + ": duplicate input file '" + file + "'.");
        }
        // Relative resource URLs are caller-relative, but child arguments must not
        // be mistaken for options (e.g. '-sample.fasta') or depend on child cwd.
        input.outputs.push_back({{"", {pathString(fs::absolute(to_path(file)))}}});
      }
      input.completed_rounds = input.total_rounds = input.outputs.size();
      return input;
    };

    // Preflight transforms all reachable rounds without invoking processing jobs.
    Result planned = result;
    // Creating the output root also lets the OS validate its naming semantics.
    // outputDirectory() checks containment before any claims are recorded.
    const auto publication_root = fs::canonical(out) / "TOPPAS_out";
    if (fs::is_symlink(publication_root)) { throw Failure(TOPPBase::ILLEGAL_PARAMETERS, "TOPPAS_out must not be a symbolic link."); }
    Internal::PipelinePathClaims planned_publications(publication_root);
    Internal::PipelinePathClaims actual_publications(publication_root);
    for (const auto id : order)
    {
      current_node = id;
      if (! affected.count(id)) continue;
      const auto& node = graph.node(id);
      auto& state = planned.nodes[id];
      if (node.kind == Kind::INPUT)
      {
        state = prepareInput(node);
        continue;
      }
      const auto edges = incoming(graph, id);
      if (edges.empty() || std::any_of(edges.begin(), edges.end(), [&](const auto* e) { return planned.nodes[e->source].state == State::BLOCKED; }))
      {
        state.state = State::BLOCKED;
        emit({Event::Type::LOG, id, 0, 0, context(node) + " is disconnected and will not run."});
        continue;
      }
      const auto inputs = inputsFor(graph, node, planned);
      if (node.kind == Kind::TOOL)
      {
        checkInputCardinality(node, descriptors.at(id), inputs);
        state.outputs = Tool::planOutputs(graph, node, inputs, descriptors.at(id), pathString(temporary.path));
      }
      else if (node.kind == Kind::MERGER || node.kind == Kind::SPLITTER)
        state.outputs = transform(node, inputs);
      else
        state.outputs = publish(graph, node, inputs, descriptors, pathString(out), true, cancelled_, emit, planned_publications);
      state.state = State::SUCCEEDED;
      state.completed_rounds = state.total_rounds = state.outputs.size();
    }
    for (const auto id : order)
    {
      if (affected.count(id) && planned.nodes[id].state == State::BLOCKED) result.nodes[id].state = State::BLOCKED;
    }

    struct Job
    {
      Size node;
      Size round;
      Tool::Invocation invocation;
    };
    struct Completion
    {
      Size node;
      Size round;
      ExternalProcess::Result process;
    };
    std::deque<Job> jobs;
    std::vector<std::future<Completion>> active;
    struct ProcessGuard
    {
      std::atomic_bool& cancelled;
      std::vector<std::future<Completion>>& active;
      void stopAndWait()
      {
        if (! active.empty()) cancelled.store(true);
        for (auto& future : active)
          if (future.valid()) future.wait();
        active.clear();
      }
      ~ProcessGuard()
      { stopAndWait(); }
    };
    std::mutex messages_mutex;
    std::condition_variable messages_ready;
    std::deque<Event> messages;
    ProcessGuard process_guard {cancelled_, active};
    auto drainMessages = [&](bool preserve_failure) {
      std::deque<Event> pending_messages;
      {
        std::lock_guard<std::mutex> lock(messages_mutex);
        pending_messages.swap(messages);
      }
      for (const auto& message : pending_messages)
      {
        if (! preserve_failure) emit(message);
        else
        {
          // A diagnostic callback must not replace the original process failure.
          try
          {
            emit(message);
          }
          catch (...)
          {
          }
        }
      }
    };
    std::set<Size> scheduled;
    for (const auto& [id, state] : result.nodes)
    {
      if (! affected.count(id) || state.state == State::BLOCKED) scheduled.insert(id);
    }
    auto snapshot = [&](Event::Type type, Size id, Size round = 0, const std::string& text = std::string {}) {
      Event event {type, id, round, result.nodes[id].total_rounds, text};
      event.result = result.nodes[id];
      emit(event);
    };
    auto complete = [&](Size id) {
      auto& state = result.nodes[id];
      state.state = State::SUCCEEDED;
      state.completed_rounds = state.total_rounds = state.outputs.size();
      snapshot(Event::Type::NODE_FINISHED, id);
    };

    try
    {
      while (true)
      {
        if (cancelled_.load()) throw Cancelled();
        for (const auto id : order)
        {
          current_node = id;
          if (scheduled.count(id)) continue;
          const auto& node = graph.node(id);
          const auto edges = incoming(graph, id);
          if (std::any_of(edges.begin(), edges.end(), [&](const auto* e) { return result.nodes[e->source].state != State::SUCCEEDED; })) continue;
          auto& state = result.nodes[id];
          scheduled.insert(id);
          if (node.kind == Kind::INPUT)
          {
            state = prepareInput(node);
            complete(id);
            continue;
          }
          const auto inputs = inputsFor(graph, node, result);
          if (node.kind == Kind::TOOL)
          {
            checkInputCardinality(node, descriptors.at(id), inputs);
            state.outputs = Tool::planOutputs(graph, node, inputs, descriptors.at(id), pathString(temporary.path));
            state.total_rounds = inputs.size();
            for (Size round = 0; round < inputs.size(); ++round)
            {
              jobs.push_back(
                {id, round,
                 Tool::invocation(graph, node, descriptors.at(id), inputs[round], state.outputs[round], round, pathString(temporary.path))});
            }
            snapshot(Event::Type::NODE_SCHEDULED, id);
          }
          else
          {
            state.state = State::RUNNING;
            snapshot(Event::Type::NODE_STARTED, id);
            if (node.kind == Kind::MERGER || node.kind == Kind::SPLITTER) state.outputs = transform(node, inputs);
            else
              state.outputs = publish(graph, node, inputs, descriptors, pathString(out), false, cancelled_, emit, actual_publications);
            complete(id);
          }
        }

        while (! jobs.empty() && active.size() < options.num_jobs)
        {
          Job job = std::move(jobs.front());
          jobs.pop_front();
          current_node = job.node;
          if (result.nodes[job.node].state != State::RUNNING)
          {
            result.nodes[job.node].state = State::RUNNING;
            snapshot(Event::Type::NODE_STARTED, job.node, job.round, context(graph.node(job.node)) + " started.");
          }
          active.push_back(std::async(std::launch::async, [&, job = std::move(job)]() {
            auto log = [&](const std::string& text) {
              std::lock_guard<std::mutex> lock(messages_mutex);
              messages.push_back({Event::Type::LOG, job.node, job.round, 0, text});
              messages_ready.notify_one();
            };
            ExternalProcess process(log, log);
            auto process_result = process.runWithResult(job.invocation.executable, job.invocation.arguments, job.invocation.working_directory, false,
                                                        ExternalProcess::IO_MODE::READ_WRITE, {}, nullptr, &cancelled_);
            messages_ready.notify_one();
            return Completion {job.node, job.round, std::move(process_result)};
          }));
        }
        drainMessages(false);

        for (Size i = 0; i < active.size();)
        {
          if (active[i].wait_for(std::chrono::milliseconds(0)) != std::future_status::ready)
          {
            ++i;
            continue;
          }
          auto completion = active[i].get();
          active.erase(active.begin() + i);
          current_node = completion.node;
          auto& state = result.nodes[completion.node];
          if (completion.process.state != ExternalProcess::RETURNSTATE::SUCCESS)
          {
            if (completion.process.state == ExternalProcess::RETURNSTATE::CANCELLED)
            {
              state.state = State::CANCELLED;
              throw Cancelled(context(graph.node(completion.node)) + ": " + completion.process.error_message);
            }
            state.state = State::FAILED;
            const int code = completion.process.state == ExternalProcess::RETURNSTATE::NONZERO_EXIT
                               ? completion.process.exit_code
                               : (completion.process.state == ExternalProcess::RETURNSTATE::FAILED_TO_START ? TOPPBase::EXTERNAL_PROGRAM_NOTFOUND
                                                                                                            : TOPPBase::EXTERNAL_PROGRAM_ERROR);
            throw Failure(code, context(graph.node(completion.node)) + ": " + completion.process.error_message);
          }
          ++state.completed_rounds;
          snapshot(Event::Type::ROUND_COMPLETED, completion.node, completion.round);
          if (state.completed_rounds == state.total_rounds)
          {
            Tool::finalizeOutputs(state.outputs);
            complete(completion.node);
          }
        }
        if (active.empty() && jobs.empty())
        {
          if (scheduled.size() == graph.nodes.size()) break;
          // Ready children are scheduled on the next pass. If none can become
          // ready, give a terminal error instead of waiting for a nonexistent event.
          const bool ready = std::any_of(order.begin(), order.end(), [&](Size id) {
            if (scheduled.count(id)) return false;
            const auto edges = incoming(graph, id);
            return std::all_of(edges.begin(), edges.end(), [&](const auto* edge) { return result.nodes[edge->source].state == State::SUCCEEDED; });
          });
          if (! ready) throw Failure(TOPPBase::INCOMPATIBLE_INPUT_DATA, "Pipeline has unfinished nodes but no runnable jobs.");
        }
        else
        {
          std::unique_lock<std::mutex> lock(messages_mutex);
          messages_ready.wait_for(lock, std::chrono::milliseconds(10), [&]() { return ! messages.empty() || cancelled_.load(); });
        }
      }
      drainMessages(false);
    }
    catch (...)
    {
      process_guard.stopAndWait();
      try
      {
        drainMessages(true);
      }
      catch (...)
      {
      }
      throw;
    }
    return result;
  }
  catch (const Cancelled& error)
  {
    cancellation_outcome = true;
    result.exit_code = error.exit_code;
    result.error_message = error.what();
  }
  catch (const Failure& error)
  {
    result.exit_code = error.exit_code;
    result.error_message = error.what();
  }
  catch (const fs::filesystem_error& error)
  {
    result.exit_code = TOPPBase::CANNOT_WRITE_OUTPUT_FILE;
    result.error_message = error.what();
  }
  catch (const std::exception& error)
  {
    result.exit_code = TOPPBase::INCOMPATIBLE_INPUT_DATA;
    result.error_message = error.what();
  }
  if (logfile.is_open())
  {
    logfile << result.error_message << '\n';
    logfile.flush();
  }
  auto failed = result.nodes.find(current_node);
  if (! cancellation_outcome && failed != result.nodes.end() && failed->second.state != State::CANCELLED)
  {
    failed->second.state = State::FAILED;
    failed->second.error_message = result.error_message;
  }
  for (auto& [id, state] : result.nodes)
  {
    if (state.state == State::RUNNING || (cancellation_outcome && state.state == State::PENDING)) { state.state = State::CANCELLED; }
    else if (state.state == State::PENDING)
      state.state = State::BLOCKED;
    if (cancellation_outcome && state.state == State::CANCELLED) state.error_message = result.error_message;
  }
  if (cancellation_outcome && (failed == result.nodes.end() || failed->second.state != State::CANCELLED))
  {
    failed = std::find_if(result.nodes.begin(), result.nodes.end(), [](const auto& item) { return item.second.state == State::CANCELLED; });
    if (failed != result.nodes.end()) current_node = failed->first;
  }
  if (callback)
  {
    Event event {Event::Type::NODE_FAILED, current_node, 0, 0, result.error_message, result.exit_code};
    if (failed != result.nodes.end()) event.result = failed->second;
    // Preserve the already established outcome even if reporting it fails.
    try
    {
      callback(event);
    }
    catch (...)
    {
    }
  }
  return result;
}
} // namespace OpenMS
