// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/PipelineExecutor.h>

#include <OpenMS/APPLICATIONS/TOPPASWorkflow.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <filesystem>
#include <set>

namespace OpenMS
{
  namespace
  {
    // on-disk parameter tags marking input/output files (see TOPPBase::TAG_*)
    const std::string TAG_INPUT_FILE = "input file";
    const std::string TAG_OUTPUT_FILE = "output file";
    const std::string TAG_OUTPUT_DIR = "output dir";

    /// 3-char zero-padded number, e.g. 7 -> "007" (matches TOPPASToolVertex output dirs).
    std::string threeChars(Int n)
    {
      std::string s = std::to_string(n);
      while (s.size() < 3) s = "0" + s;
      return s;
    }

    /// Extract the child exit code embedded in ExternalProcess's error message ("...(exit code: N)..."),
    /// or -1 if it cannot be found. Used to recognise a tool that is missing its own external program
    /// (TOPP tools exit with EXTERNAL_PROGRAM_NOTFOUND == 14 in that case).
    int childExitCode(const std::string& msg)
    {
      const std::string marker = "(exit code: ";
      std::string::size_type pos = msg.find(marker);
      if (pos == std::string::npos) return -1;
      pos += marker.size();
      std::string::size_type end = msg.find(')', pos);
      if (end == std::string::npos) return -1;
      try
      {
        return std::stoi(msg.substr(pos, end - pos));
      }
      catch (...)
      {
        return -1;
      }
    }
  } // namespace

  std::vector<PipelineExecutor::IOParam> PipelineExecutor::collectIO_(const Param& param, bool input)
  {
    // tuple of (sort_type, IOParam) where sort_type mirrors TOPPASToolVertex::IOInfo::IOType
    // (0 = file, 1 = list, 2 = dir); IOInfo::operator< sorts files first, then by name.
    std::vector<std::pair<int, IOParam>> tmp;

    auto add = [&](const std::string& tag) {
      for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
      {
        if (!it->tags.count(tag)) continue;

        IOParam e;
        e.name = it.getName();
        // first valid extension (restrictions look like "*.mzML")
        for (const std::string& v : it->valid_strings)
        {
          if (StringUtils::hasPrefix(v, "*."))
          {
            e.ext = v.substr(2);
            break;
          }
        }
        int sort_type;
        if (it->value.valueType() == ParamValue::STRING_LIST)
        {
          e.is_list = true;
          sort_type = 1; // IOT_LIST
        }
        else
        {
          e.is_list = false;
          sort_type = (tag == TAG_OUTPUT_DIR) ? 2 /*IOT_DIR*/ : 0 /*IOT_FILE*/;
        }
        tmp.emplace_back(sort_type, e);
      }
    };

    if (input)
    {
      add(TAG_INPUT_FILE);
    }
    else
    {
      add(TAG_OUTPUT_FILE);
      add(TAG_OUTPUT_DIR);
    }

    // reproduce IOInfo::operator<: differing types -> the IOT_FILE (0) one is "smaller"; else by name
    std::stable_sort(tmp.begin(), tmp.end(), [](const std::pair<int, IOParam>& a, const std::pair<int, IOParam>& b) {
      if (a.first != b.first) return a.first == 0;
      return a.second.name < b.second.name;
    });

    std::vector<IOParam> result;
    result.reserve(tmp.size());
    for (auto& p : tmp) result.push_back(p.second);
    return result;
  }

  void PipelineExecutor::setFileParam_(Param& p, const std::string& name, const std::vector<std::string>& files, bool is_list)
  {
    if (is_list)
    {
      p.setValue(name, files);
    }
    else
    {
      if (files.size() > 1)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Multiple files given to single-file parameter '" + name + "'");
      }
      p.setValue(name, files.empty() ? std::string() : files[0]);
    }
  }

  void PipelineExecutor::fromWorkflow(const TOPPASWorkflow& wf, const std::string& toppas_file)
  {
    nodes_.clear();
    edges_.clear();
    tools_.clear();

    toppas_stem_ = File::stemName(toppas_file);
    if (toppas_stem_.empty()) toppas_stem_ = "Untitled_workflow";
    const std::string toppas_dir = File::path(toppas_file);

    std::map<std::string, int> id2idx;

    // --- nodes ---
    for (const TOPPASWorkflow::Node& wn : wf.nodes)
    {
      Node n;
      n.id = wn.id;
      n.recycle_output = wn.recycle_output;
      ToolDesc td;

      switch (wn.type)
      {
        case TOPPASWorkflow::NodeType::INPUT_FILE_LIST:
        {
          n.kind = NodeKind::INPUT;
          for (const std::string& f : wn.file_names)
          {
            // resolve relative input paths against the .toppas location, then normalize ('../' etc.)
            std::filesystem::path p(f);
            if (p.is_relative() && !toppas_dir.empty())
            {
              p = std::filesystem::path(toppas_dir) / p;
            }
            n.input_files.push_back(p.lexically_normal().string());
          }
          break;
        }
        case TOPPASWorkflow::NodeType::OUTPUT_FILE_LIST: n.kind = NodeKind::OUTPUT; break;
        case TOPPASWorkflow::NodeType::OUTPUT_FOLDER: n.kind = NodeKind::OUTPUT_FOLDER; break;
        case TOPPASWorkflow::NodeType::MERGER:
          n.kind = NodeKind::MERGER;
          n.round_based_merge = wn.round_based;
          break;
        case TOPPASWorkflow::NodeType::SPLITTER: n.kind = NodeKind::SPLITTER; break;
        case TOPPASWorkflow::NodeType::TOOL:
          n.kind = NodeKind::TOOL;
          td.is_tool = true;
          td.name = wn.tool_name;
          td.type = wn.tool_type;
          td.param = wn.parameters;
          td.in = collectIO_(wn.parameters, true);
          td.out = collectIO_(wn.parameters, false);
          break;
        default:
          continue; // unknown node type -> skip (keeps nodes_/tools_ aligned)
      }

      int idx = addNode(n);
      tools_.push_back(td);
      id2idx[wn.id] = idx;
    }

    // --- edges (resolve parameter names to indices for tool endpoints) ---
    for (const TOPPASWorkflow::Edge& we : wf.edges)
    {
      std::map<std::string, int>::const_iterator sit = id2idx.find(we.source_id);
      std::map<std::string, int>::const_iterator tit = id2idx.find(we.target_id);
      if (sit == id2idx.end() || tit == id2idx.end()) continue;
      int from = sit->second, to = tit->second;

      Int src_index = -1, tgt_index = -1;
      if (nodes_[from].kind == NodeKind::TOOL && we.source_out_param != TOPPASWorkflow::NO_NAME)
      {
        const std::vector<IOParam>& outs = tools_[from].out;
        for (Size k = 0; k < outs.size(); ++k)
        {
          if (outs[k].name == we.source_out_param) { src_index = (Int)k; break; }
        }
      }
      if (nodes_[to].kind == NodeKind::TOOL && we.target_in_param != TOPPASWorkflow::NO_NAME)
      {
        const std::vector<IOParam>& ins = tools_[to].in;
        for (Size k = 0; k < ins.size(); ++k)
        {
          if (ins[k].name == we.target_in_param) { tgt_index = (Int)k; break; }
        }
      }
      addEdge(from, to, src_index, tgt_index);
    }
  }

  PipelineExecutor::Result PipelineExecutor::execute(const Options& opts,
                                                     const std::map<std::string, std::vector<std::string>>& input_overrides,
                                                     std::string& error_msg)
  {
    opts_ = opts;

    if (!topoSort(error_msg)) return INVALID;

    // replace input-node files by the resource overrides (keyed by topological number)
    if (!input_overrides.empty())
    {
      for (Node& n : nodes_)
      {
        if (n.kind != NodeKind::INPUT) continue;
        std::map<std::string, std::vector<std::string>>::const_iterator it = input_overrides.find(std::to_string(n.topo_nr));
        if (it != input_overrides.end()) n.input_files = it->second;
      }
    }

    try
    {
      if (!computeDataFlow(error_msg)) return INVALID;
    }
    catch (const ExecError& e)
    {
      error_msg = e.msg;
      return e.result;
    }
    catch (const Exception::BaseException& e)
    {
      error_msg = e.what();
      return TOOL_FAILED;
    }
    catch (const std::exception& e) // e.g. std::filesystem_error -- must not escape to std::terminate
    {
      error_msg = e.what();
      return TOOL_FAILED;
    }
    return OK;
  }

  void PipelineExecutor::processToolNode_(int node_index, const RoundPackages& inputs)
  {
    Node& node = nodes_[node_index];
    const ToolDesc& td = tools_[node_index];
    const int rounds = node.round_total;

    node.output_files.assign(rounds, RoundPackage());

    // per-tool output directory: <tmp>/<workflow>/<topo>_<tool>[_<type>]
    std::string out_dir = opts_.tmp_dir + "/" + toppas_stem_ + "/" + threeChars(node.topo_nr) + "_" + td.name;
    if (!td.type.empty()) out_dir += "_" + td.type;
    if (!File::makeDir(out_dir) && !File::isDirectory(out_dir))
    {
      throw ExecError{TOOL_FAILED, "Could not create working directory '" + out_dir + "'"};
    }

    // which output parameters actually have an outgoing edge (only those need files)
    std::set<int> wanted_outputs;
    for (int ei : node.out_edges)
    {
      int p = (int)edges_[ei].source_out_param;
      if (p >= 0) wanted_outputs.insert(p);
    }

    // resolve the tool executable once
    std::string exe;
    try
    {
      exe = File::findSiblingTOPPExecutable(td.name);
    }
    catch (Exception::FileNotFound&)
    {
      throw ExecError{TOOL_NOT_FOUND, "Tool executable not found: " + td.name};
    }

    for (int round = 0; round < rounds; ++round)
    {
      Param ptmp = td.param;

      // set the input file parameters from the incoming round package
      const RoundPackage& in_pkg = inputs[round];
      Size ref_count = 1; // how many output files a list-output produces (mirror the largest input)
      for (RoundPackage::const_iterator it = in_pkg.begin(); it != in_pkg.end(); ++it)
      {
        int pidx = (int)it->first;
        if (pidx < 0 || pidx >= (int)td.in.size()) continue; // negative keys only occur for parameter-less nodes
        const IOParam& ip = td.in[pidx];
        setFileParam_(ptmp, ip.name, it->second.filenames, ip.is_list);
        // size list outputs from NON-recycling inputs only (mirror TOPPASToolVertex::updateCurrentOutputFileNames)
        bool from_recycler = it->second.edge_index >= 0 && nodes_[edges_[it->second.edge_index].from].recycle_output;
        if (!from_recycler) ref_count = std::max(ref_count, it->second.filenames.size());
      }

      // generate the output file names for every output parameter that feeds a downstream node
      for (int pidx : wanted_outputs)
      {
        if (pidx >= (int)td.out.size()) continue;
        const IOParam& op = td.out[pidx];
        const Size count = op.is_list ? ref_count : 1;
        std::string safe_name = op.name; // parameter names may contain ':' (nested sections)
        std::replace(safe_name.begin(), safe_name.end(), ':', '_');
        std::vector<std::string> ofiles;
        ofiles.reserve(count);
        for (Size k = 0; k < count; ++k)
        {
          std::string fn = out_dir + "/" + safe_name + "_r" + std::to_string(round) + "_" + std::to_string(k);
          if (!op.ext.empty()) fn += "." + op.ext;
          ofiles.push_back(fn);
        }
        setFileParam_(ptmp, op.name, ofiles, op.is_list);
        node.output_files[round][pidx].filenames = ofiles;
      }

      // write the per-round INI (TOPP tools expect "<tool>:1:<params>")
      Param ini;
      ini.setValue(td.name + ":1:toppas_dummy", "x");
      ini.insert(td.name + ":1:", ptmp);
      ini.remove(td.name + ":1:toppas_dummy");
      std::string ini_file = out_dir + "/" + td.name + (td.type.empty() ? "" : "_" + td.type) + "_r" + std::to_string(round) + ".ini";
      ParamXMLFile().store(ini_file, ini);

      // run the tool
      std::vector<std::string> args;
      if (!td.type.empty())
      {
        args.push_back("-type");
        args.push_back(td.type);
      }
      args.push_back("-ini");
      args.push_back(ini_file);

      ExternalProcess proc;
      std::string proc_err;
      ExternalProcess::RETURNSTATE state = proc.run(exe, args, "", false, proc_err);
      if (state == ExternalProcess::RETURNSTATE::FAILED_TO_START)
      {
        throw ExecError{TOOL_NOT_FOUND, "Could not start tool '" + td.name + "': " + proc_err};
      }
      if (state != ExternalProcess::RETURNSTATE::SUCCESS)
      {
        // a TOPP tool that cannot find its OWN external program exits with EXTERNAL_PROGRAM_NOTFOUND (14);
        // propagate that as TOOL_NOT_FOUND so a pipeline using e.g. CometAdapter without comet SKIPS (matches
        // the historical behavior, where the child exit code was passed through)
        if (childExitCode(proc_err) == 14)
        {
          throw ExecError{TOOL_NOT_FOUND, "Tool '" + td.name + "' is missing a required external program: " + proc_err};
        }
        throw ExecError{TOOL_FAILED, "Tool '" + td.name + "' failed: " + proc_err};
      }
    }
  }

  void PipelineExecutor::processOutputNode_(int node_index, const RoundPackages& inputs)
  {
    Node& node = nodes_[node_index];
    node.output_files.assign(node.round_total, RoundPackage());

    // each output node writes into its own sub-directory so that files of different output nodes (which
    // commonly share a parameter name like 'out') never clobber each other in the shared output directory
    std::string node_out_dir = opts_.out_dir + "/" + threeChars(node.topo_nr);
    if (!File::makeDir(node_out_dir) && !File::isDirectory(node_out_dir))
    {
      throw ExecError{TOOL_FAILED, "Could not create output directory '" + node_out_dir + "'"};
    }

    for (int round = 0; round < node.round_total; ++round)
    {
      const RoundPackage& pkg = inputs[round];
      for (RoundPackage::const_iterator it = pkg.begin(); it != pkg.end(); ++it)
      {
        for (const std::string& from : it->second.filenames)
        {
          namespace fs = std::filesystem;
          std::error_code exist_ec;
          if (!fs::exists(from, exist_ec))
          {
            throw ExecError{TOOL_FAILED, "Expected output does not exist: '" + from + "'"};
          }
          std::string to = node_out_dir + "/" + File::basename(from);
          std::error_code ec;
          if (fs::is_directory(from, exist_ec)) // output-directory parameter: copy the whole tree
          {
            fs::remove_all(to, ec);
            fs::copy(from, to, fs::copy_options::recursive | fs::copy_options::overwrite_existing, ec);
            if (ec) throw ExecError{TOOL_FAILED, "Could not copy directory '" + from + "' to '" + to + "': " + ec.message()};
          }
          else
          {
            if (File::exists(to)) File::remove(to);
            if (!File::copy(from, to))
            {
              throw ExecError{TOOL_FAILED, "Could not copy '" + from + "' to '" + to + "'"};
            }
          }
          node.output_files[round][-1].filenames.push_back(to);
        }
      }
    }
  }

} // namespace OpenMS
