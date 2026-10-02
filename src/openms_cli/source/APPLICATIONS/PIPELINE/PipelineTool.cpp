// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "PipelineTool.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileNameUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <algorithm>
#include <cctype>
#include <charconv>
#include <filesystem>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <utility>

#ifdef _WIN32
  #include <Windows.h>
#endif

namespace OpenMS::Internal
{
namespace
{
  namespace fs = std::filesystem;
  using Port = PipelineTool::Port;

  struct PathLess
  {
    bool operator()(const fs::path& lhs, const fs::path& rhs) const
    {
#ifdef _WIN32
      // Compare native UTF-16, including non-ASCII case pairs. Windows accepts
      // both separator styles and generally aliases paths which differ only in case.
      const auto left = lhs.lexically_normal().make_preferred().native();
      const auto right = rhs.lexically_normal().make_preferred().native();
      return CompareStringOrdinal(left.c_str(), -1, right.c_str(), -1, TRUE) == CSTR_LESS_THAN;
#else
      return lhs < rhs;
#endif
    }
  };

  std::string pathString(const fs::path& path)
  {
    const auto value = path.u8string();
    return std::string(reinterpret_cast<const char*>(value.data()), value.size());
  }

  std::string paddedNumber(Size number)
  {
    const std::string text = std::to_string(number);
    return std::string(text.size() < 3 ? 3 - text.size() : 0, '0') + text;
  }

  // Parameters and tool types are user-editable; never let them escape the run directory.
  std::string pathComponent(std::string text)
  {
    for (char& c : text)
    {
      if (c == '/' || c == '\\' || c == ':' || c == '\0') { c = '_'; }
    }
    if (text.empty() || text == "." || text == "..") { text = "unnamed"; }
    return text;
  }

  fs::path toolDirectory(const PipelineGraph& graph, const PipelineGraph::Node& node, const std::string& run_temp)
  {
    std::string workflow = File::stemName(graph.filename);
    if (workflow.empty()) { workflow = "Untitled_workflow"; }
    std::string tool = paddedNumber(node.topo_number) + "_" + node.tool_name;
    if (! node.tool_type.empty()) { tool += "_" + node.tool_type; }
    return to_path(run_temp) / to_path(pathComponent(workflow)) / to_path(pathComponent(tool));
  }

  void writeParameters(const PipelineGraph::Node& node, const Param& parameters, const fs::path& filename)
  {
    Param save;
    save.insert(node.tool_name + ":1:", parameters);
    save.setSectionDescription(node.tool_name + ":1", "Instance '1' section for '" + node.tool_name + "'");
    ParamXMLFile().store(pathString(filename), save);
  }

  std::vector<Port> ports(const Param& parameters, bool inputs)
  {
    std::vector<Port> result;
    const std::vector<std::string> tags
      = inputs ? std::vector<std::string> {TOPPBase::TAG_INPUT_FILE} : std::vector<std::string> {TOPPBase::TAG_OUTPUT_FILE, TOPPBase::TAG_OUTPUT_DIR};
    for (const auto& tag : tags)
    {
      for (auto it = parameters.begin(); it != parameters.end(); ++it)
      {
        if (! it->tags.count(tag)) { continue; }
        Port port;
        port.name = it.getName();
        if (it->value.valueType() == ParamValue::STRING_LIST) { port.kind = Port::Kind::LIST; }
        else if (it->value.valueType() == ParamValue::STRING_VALUE)
        {
          port.kind = tag == TOPPBase::TAG_OUTPUT_DIR ? Port::Kind::DIRECTORY : Port::Kind::FILE;
        }
        else
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Non-string workflow file parameter '" + port.name + "'.");
        }
        for (const auto& restriction : it->valid_strings)
        {
          port.valid_types.push_back(restriction.starts_with("*.") ? restriction.substr(2) : restriction);
        }
        result.push_back(std::move(port));
      }
    }
    // Keep the historical file-before-list order, and make mixed list/directory ordering strict.
    std::sort(result.begin(), result.end(), [](const Port& lhs, const Port& rhs) {
      if (lhs.kind != rhs.kind) { return lhs.kind < rhs.kind; }
      return lhs.name < rhs.name;
    });
    return result;
  }

  void smartFileNames(std::vector<std::vector<std::string>>& names)
  {
    if (names.size() < 2 || names.front().size() != 1) { return; }
    const auto basename = to_path(names.front().front()).filename();
    for (const auto& round : names)
    {
      if (round.size() != 1 || to_path(round.front()).filename() != basename) { return; }
    }
    // When repeated basenames come from different samples, name the outputs after sample folders.
    for (auto& round : names)
    {
      std::error_code ec;
      // Preflight names refer to output folders which have not been created yet.
      const auto parent = fs::weakly_canonical(to_path(round.front()).parent_path(), ec);
      if (ec) { continue; }
      const auto candidate = pathString(parent.filename());
      if (candidate.size() > 2 && candidate.find(':') == std::string::npos) { round.front() = candidate; }
    }
  }

  std::string outputSuffix(const PipelineGraph::Node& node, const PipelineTool::Descriptor& descriptor, const Port& port)
  {
    if (port.kind == Port::Kind::DIRECTORY) { return "_dir"; }
    if (port.valid_types.size() == 1)
    {
      const auto type = FileTypes::nameToType(port.valid_types.front());
      return "." + (type == FileTypes::UNKNOWN ? port.valid_types.front() : FileTypes::typeToName(type));
    }
    if (const auto type_parameter = port.name + "_type"; descriptor.parameters.exists(type_parameter))
    {
      const auto type = descriptor.parameters.getValue(type_parameter).toString();
      if (! type.empty()) { return "." + type; }
    }
    if (node.tool_name == "FileMerger" && descriptor.parameters.exists("in_type"))
    {
      const auto type = descriptor.parameters.getValue("in_type").toString();
      if (! type.empty()) { return "." + type; }
    }
    if (! port.valid_types.empty())
    {
      // A writer with several advertised formats may dispatch by extension and
      // reject .unknown before there is any content to inspect. Preserve the
      // schema's declared preference (e.g. idXML before idparquet in adapters).
      const auto type = FileTypes::nameToType(port.valid_types.front());
      return "." + (type == FileTypes::UNKNOWN ? port.valid_types.front() : FileTypes::typeToName(type));
    }
    return ".unknown";
  }

  bool endsWithIgnoringCase(const std::string& text, const std::string& suffix)
  {
    return text.size() >= suffix.size() && std::equal(suffix.rbegin(), suffix.rend(), text.rbegin(), [](unsigned char lhs, unsigned char rhs) {
             return std::tolower(lhs) == std::tolower(rhs);
           });
  }
} // namespace

std::string PipelineTool::temporaryName(const std::string& prefix)
{
  // Hostnames and timestamps can consume most of MAX_PATH once repeated in a
  // run directory and its filename probe. Hash their unique token instead;
  // exclusive creation remains responsible for handling any name collision.
  const auto value = std::hash<std::string> {}(File::getUniqueName(false));
  char buffer[2 * sizeof(value)];
  const auto converted = std::to_chars(std::begin(buffer), std::end(buffer), value, 16);
  return prefix + std::string(buffer, converted.ptr);
}

PipelinePathClaims::PipelinePathClaims(const fs::path& root): root_(fs::weakly_canonical(fs::absolute(root)))
{
  fs::create_directories(root_);
  for (Size attempt = 0; attempt < 10; ++attempt)
  {
    const auto candidate = root_ / PipelineTool::temporaryName(".pc-");
    if (fs::create_directory(candidate))
    {
      probe_ = candidate;
      return;
    }
  }
  throw PipelineTool::Failure(TOPPBase::CANNOT_WRITE_OUTPUT_FILE, "Could not create a private filename probe directory.");
}

PipelinePathClaims::~PipelinePathClaims()
{
  std::error_code ignored;
  if (! probe_.empty()) { fs::remove_all(probe_, ignored); }
}

bool PipelinePathClaims::claim(const fs::path& path)
{
  const auto relative = fs::weakly_canonical(fs::absolute(path)).lexically_relative(root_);
  if (relative.empty() || relative == "." || *relative.begin() == "..")
  {
    throw PipelineTool::Failure(TOPPBase::CANNOT_WRITE_OUTPUT_FILE, "Output filename escapes its destination directory.");
  }
  auto marker = probe_;
  for (auto component = relative.begin(); component != relative.end(); ++component)
  {
    marker /= *component;
    // Reserved leaves remain empty; a subsequent output cannot use one as a parent.
    if (std::next(component) != relative.end() && fs::exists(marker) && fs::is_empty(marker)) { return false; }
  }
  fs::create_directories(marker.parent_path());
  // The OS resolves case, Unicode normalization, and Windows trailing dots here.
  // An existing marker also catches a destination that is another output's parent.
  return fs::create_directory(marker);
}

PipelineTool::Descriptor PipelineTool::discover(const PipelineGraph::Node& node,
                                                const std::string& temp_directory,
                                                const std::atomic_bool& cancelled,
                                                std::function<void(const std::string&)> log)
{
  Descriptor descriptor;
  try
  {
    descriptor.executable = File::findSiblingTOPPExecutable(node.tool_name);
  }
  catch (const Exception::FileNotFound&)
  {
    descriptor.executable = node.tool_name;
    if (! File::findExecutable(descriptor.executable))
    {
      throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, TOPPBase::EXTERNAL_PROGRAM_NOTFOUND,
                    "Could not find workflow tool '" + node.tool_name + "'.");
    }
  }

  const auto directory = to_path(temp_directory) / ("parameters_" + std::to_string(node.id));
  fs::create_directories(directory);
  const auto filename = directory / "defaults.ini";
  std::vector<std::string> arguments {"-write_ini", pathString(filename)};
  if (! node.tool_type.empty()) { arguments.insert(arguments.end(), {"-type", node.tool_type}); }
  ExternalProcess process(log, log);
  const auto result
    = process.runWithResult(descriptor.executable, arguments, "", false, ExternalProcess::IO_MODE::READ_ONLY, {}, nullptr, &cancelled);
  if (result.state != ExternalProcess::RETURNSTATE::SUCCESS)
  {
    const int code = result.state == ExternalProcess::RETURNSTATE::FAILED_TO_START
                       ? TOPPBase::EXTERNAL_PROGRAM_NOTFOUND
                       : (result.exit_code > 0 ? result.exit_code : TOPPBase::EXTERNAL_PROGRAM_ERROR);
    throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, code,
                  "Parameter discovery failed for '" + node.tool_name + "': " + result.error_message);
  }
  if (! fs::is_regular_file(filename))
  {
    throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, TOPPBase::EXTERNAL_PROGRAM_ERROR,
                  "Tool '" + node.tool_name + "' did not create its parameter file.");
  }
  Param defaults;
  ParamXMLFile().load(pathString(filename), defaults);
  descriptor.parameters = defaults.copy(node.tool_name + ":1:", true);
  if (descriptor.parameters.empty())
  {
    throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, TOPPBase::ILLEGAL_PARAMETERS,
                  "Tool '" + node.tool_name + "' returned an empty parameter schema.");
  }
  if (descriptor.parameters.exists("no_progress"))
  {
    descriptor.parameters.setValue("no_progress", "true", descriptor.parameters.getDescription("no_progress"),
                                   descriptor.parameters.getTags("no_progress"));
  }
  // Refresh metadata from the executable while preserving saved values, including the user's no_progress choice.
  if (! descriptor.parameters.update(node.parameters, false, false, true, true, OPENMS_LOG_WARN))
  {
    throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, TOPPBase::ILLEGAL_PARAMETERS,
                  "Saved parameters for '" + node.tool_name + "' are incompatible with its current schema.");
  }
  descriptor.inputs = ports(descriptor.parameters, true);
  descriptor.outputs = ports(descriptor.parameters, false);
  return descriptor;
}

PipelineTool::Rounds PipelineTool::planOutputs(const PipelineGraph& graph,
                                               const PipelineGraph::Node& node,
                                               const Rounds& inputs,
                                               const Descriptor& descriptor,
                                               const std::string& run_temp)
{
  if (inputs.empty()) { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "A tool requires at least one input round."); }
  const bool single_outputs
    = std::none_of(descriptor.outputs.begin(), descriptor.outputs.end(), [](const Port& port) { return port.kind == Port::Kind::LIST; });
  const Port* naming_port = nullptr;
  Size max_size = 0;
  for (const auto& port : descriptor.inputs)
  {
    const auto input = inputs.front().find(port.name);
    if (input == inputs.front().end() || input->second.empty()) { continue; }
    const auto edge = std::find_if(graph.edges.begin(), graph.edges.end(), [&](const PipelineGraph::Edge& candidate) {
      return candidate.target == node.id && candidate.target_port == port.name;
    });
    if (edge == graph.edges.end() || graph.node(edge->source).recycle_output) { continue; }
    const Size size = input->second.size();
    if ((single_outputs && (! naming_port || port.name == "in" || size == 1))
        || (! single_outputs && (! naming_port || size > max_size || (size == max_size && port.name == "in"))))
    {
      naming_port = &port;
      max_size = size;
    }
  }
  if (! naming_port)
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No non-recycled input can supply output filenames for '" + node.tool_name + "'.");
  }

  std::vector<std::vector<std::string>> names;
  for (const auto& round : inputs)
  {
    const auto input = round.find(naming_port->name);
    if (input == round.end() || input->second.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "An input round is missing files for '" + naming_port->name + "'.");
    }
    names.emplace_back();
    for (const auto& filename : input->second)
    {
      names.back().push_back(FileNameUtils::stripExtension(filename));
    }
  }
  smartFileNames(names);

  Rounds outputs(inputs.size());
  PipelinePathClaims allocated(to_path(run_temp));
  const auto directory = toolDirectory(graph, node, run_temp);
  for (const auto& port : descriptor.outputs)
  {
    if (std::none_of(graph.edges.begin(), graph.edges.end(),
                     [&](const PipelineGraph::Edge& edge) { return edge.source == node.id && edge.source_port == port.name; }))
    {
      continue;
    }
    auto component = port.name;
    component.erase(std::remove(component.begin(), component.end(), ':'), component.end());
    const auto path = directory / to_path(pathComponent(component.substr(0, 50)));
    const auto suffix = outputSuffix(node, descriptor, port);
    for (Size round = 0; round < inputs.size(); ++round)
    {
      const bool list_to_single = names[round].size() > 1 && port.kind == Port::Kind::FILE;
      for (const auto& input_file : names[round])
      {
        std::string basename = pathString(to_path(input_file).filename());
        fs::path filename;
        if (port.kind == Port::Kind::DIRECTORY)
        {
          // QFileInfo::baseName(), historically used here, takes the part before the first dot.
          filename = inputs.size() > 1 ? path / to_path(basename.substr(0, basename.find('.'))) : path;
        }
        else
        {
          if (list_to_single)
          {
            if (const auto start = basename.find("_to_"); start != std::string::npos && basename.find("_mrgd", start + 4) != std::string::npos)
            {
              basename.resize(start);
            }
            std::string last = pathString(to_path(names[round].back()).filename());
            if (const auto start = last.find("_to_"); start != std::string::npos)
            {
              const auto end = last.find("_mrgd", start + 4);
              if (end != std::string::npos) { last = last.substr(start + 4, end - start - 4); }
            }
            basename += "_to_" + last + "_mrgd";
          }
          if (! basename.ends_with(suffix)) { basename += suffix; }
          filename = path / to_path(basename);
        }
        if (! allocated.claim(filename))
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Workflow output filenames collide at '" + pathString(filename) + "'.");
        }
        outputs[round][port.name].push_back(pathString(filename));
        if (list_to_single || port.kind == Port::Kind::DIRECTORY) { break; }
      }
    }
  }
  return outputs;
}

PipelineTool::Invocation PipelineTool::invocation(const PipelineGraph& graph,
                                                  const PipelineGraph::Node& node,
                                                  const Descriptor& descriptor,
                                                  const FileBundle& inputs,
                                                  const FileBundle& outputs,
                                                  Size round,
                                                  const std::string& run_temp)
{
  Invocation invocation;
  invocation.executable = descriptor.executable;
  if (! node.tool_type.empty()) { invocation.arguments.insert(invocation.arguments.end(), {"-type", node.tool_type}); }
  Param parameters = descriptor.parameters;
  const auto directory = toolDirectory(graph, node, run_temp);
  fs::create_directories(directory);

  auto bind = [&](const FileBundle& bundle, const std::vector<Port>& ports, bool output) {
    for (const auto& port : ports)
    {
      const auto bound = bundle.find(port.name);
      if (bound == bundle.end()) { continue; }
      const auto& files = bound->second;
      if (files.empty() || (port.kind != Port::Kind::LIST && files.size() != 1))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Invalid number of files for workflow parameter '" + port.name + "'.");
      }
      if (output)
      {
        for (const auto& filename : files)
        {
          fs::create_directories(to_path(filename).parent_path());
        }
      }
      // Nested GenericWrapper parameters and long lists must use INI binding (Windows command-line limits).
      if (port.name.starts_with("ETool:") || files.size() > 10)
      {
        const ParamValue value = port.kind == Port::Kind::LIST ? ParamValue(files) : ParamValue(files.front());
        parameters.setValue(port.name, value, parameters.getDescription(port.name), parameters.getTags(port.name));
      }
      else
      {
        invocation.arguments.push_back("-" + port.name);
        invocation.arguments.insert(invocation.arguments.end(), files.begin(), files.end());
      }
    }
  };
  bind(inputs, descriptor.inputs, false);
  bind(outputs, descriptor.outputs, true);
  // Always use a private INI per round: workers may prepare and execute different rounds concurrently.
  const auto ini = directory / to_path(pathComponent(node.tool_name) + "_" + std::to_string(round) + ".ini");
  writeParameters(node, parameters, ini);
  invocation.arguments.insert(invocation.arguments.end(), {"-ini", pathString(ini)});
  return invocation;
}

void PipelineTool::finalizeOutputs(Rounds& outputs)
{
  struct Rename
  {
    fs::path source;
    fs::path target;
    fs::path staged;
    bool committed {false};
  };
  std::map<fs::path, std::string, PathLess> replacements;
  std::vector<Rename> renames;
  std::map<fs::path, Size, PathLess> counts;
  std::set<fs::path, PathLess> sources;
  for (const auto& round : outputs)
  {
    for (const auto& [port, files] : round)
    {
      for (const auto& filename : files)
      {
        if (! sources.insert(to_path(filename)).second) { continue; }
        if (! fs::exists(to_path(filename)))
        {
          throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, TOPPBase::UNEXPECTED_RESULT,
                        "Workflow tool did not create expected output '" + filename + "'.");
        }
        if (fs::is_directory(to_path(filename))) { continue; }
        std::string target = filename;
        const auto detected = FileHandler::getTypeByContent(filename);
        // Content sniffing cannot identify every advertised or custom format.
        // Keep the planned extension unless it supplies positive evidence.
        if (detected != FileTypes::UNKNOWN)
        {
          std::string suffix = FileTypes::typeToName(detected);
          if (endsWithIgnoringCase(filename, suffix)) { suffix = filename.substr(filename.size() - suffix.size()); }
          // Keep compression visible to downstream tools; stripExtension() removes both format and compression.
          if (const auto compression = FileNameUtils::compressionType(filename); compression != FileTypes::UNKNOWN)
          {
            suffix += "." + FileTypes::typeToName(compression);
          }
          target = FileNameUtils::stripExtension(filename) + "." + suffix;
        }
        renames.push_back({to_path(filename), to_path(target), {}, false});
        ++counts[to_path(target)];
      }
    }
  }
  std::map<fs::path, Size, PathLess> counters;
  std::set<fs::path, PathLess> targets;
  std::map<fs::path, std::unique_ptr<PipelinePathClaims>, PathLess> target_claims;
  for (auto& rename : renames)
  {
    const auto unnumbered = pathString(rename.target);
    if (counts[rename.target] > 1)
    {
      const auto prefix = FileNameUtils::stripExtension(unnumbered);
      rename.target = to_path(prefix + "_" + paddedNumber(++counters[rename.target]) + unnumbered.substr(prefix.size()));
    }
    auto& claims = target_claims[rename.target.parent_path()];
    if (! claims) { claims = std::make_unique<PipelinePathClaims>(rename.target.parent_path()); }
    if (! claims->claim(rename.target))
    {
      throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, TOPPBase::CANNOT_WRITE_OUTPUT_FILE,
                    "Content-based output filenames collide at '" + pathString(rename.target) + "'.");
    }
    targets.insert(rename.target);
    std::error_code ec;
    if (rename.source == rename.target || fs::equivalent(rename.source, rename.target, ec)) { rename.target = rename.source; }
    else if (fs::exists(rename.target) && ! sources.count(rename.target))
    {
      throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, TOPPBase::CANNOT_WRITE_OUTPUT_FILE,
                    "Refusing to overwrite an existing file while finalizing '" + pathString(rename.target) + "'.");
    }
    replacements[rename.source] = pathString(rename.target);
  }

  // Two-phase renaming handles swaps and suffix collisions without deleting a source file.
  try
  {
    for (auto& rename : renames)
    {
      if (rename.source == rename.target) { continue; }
      fs::path candidate;
      do
      {
        candidate = rename.source.parent_path() / temporaryName(".rn-");
      } while (fs::exists(candidate) || targets.count(candidate));
      fs::rename(rename.source, candidate);
      rename.staged = std::move(candidate);
    }
    for (auto& rename : renames)
    {
      if (rename.staged.empty()) { continue; }
      fs::rename(rename.staged, rename.target);
      rename.committed = true;
    }
  }
  catch (const fs::filesystem_error& error)
  {
    for (auto& rename : renames)
    {
      std::error_code ignored;
      if (rename.committed) { fs::rename(rename.target, rename.staged, ignored); }
    }
    for (auto& rename : renames)
    {
      std::error_code ignored;
      if (! rename.staged.empty()) { fs::rename(rename.staged, rename.source, ignored); }
    }
    throw Failure(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, TOPPBase::CANNOT_WRITE_OUTPUT_FILE,
                  "Could not finalize workflow outputs: " + std::string(error.what()));
  }
  for (auto& round : outputs)
  {
    for (auto& [port, files] : round)
    {
      for (auto& filename : files)
      {
        if (const auto found = replacements.find(to_path(filename)); found != replacements.end()) { filename = found->second; }
      }
    }
  }
}
} // namespace OpenMS::Internal
