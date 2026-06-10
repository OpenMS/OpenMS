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
#include <atomic>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <functional>
#include <mutex>
#include <set>
#include <thread>

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

  std::vector<PipelineExecutor::IOParam> PipelineExecutor::collectIO(const Param& param, bool input)
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

  std::map<std::string, std::vector<std::string>> PipelineExecutor::loadResourceFile(const std::string& filename)
  {
    std::map<std::string, std::vector<std::string>> result;
    Param p;
    ParamXMLFile().load(filename, p);
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      std::string key = it.getName(); // "<node>:url_list"
      std::string::size_type pos = key.rfind(":url_list");
      if (pos == std::string::npos) continue;
      std::string node = key.substr(0, pos);
      StringList urls = ListUtils::toStringList<std::string>(it->value);
      std::vector<std::string> files;
      for (std::string u : urls)
      {
        // strip a 'file:' / 'file://' URL scheme to get a local path
        if (StringUtils::hasPrefix(u, "file://")) u = u.substr(7);
        else if (StringUtils::hasPrefix(u, "file:")) u = u.substr(5);
        files.push_back(u);
      }
      result[node] = files;
    }
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
          td.in = collectIO(wn.parameters, true);
          td.out = collectIO(wn.parameters, false);
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

    // ------------------------------------------------------------------------------------------------
    // Global, dynamic, process-level scheduler (reproduces TOPPASScene's allowed_threads_ queue):
    // up to num_jobs TOOL processes run concurrently across the whole DAG. All graph logic runs on this
    // (coordinator) thread; worker threads only launch the blocking ExternalProcess subprocesses. The
    // "instant" nodes (input/merger/splitter/output) are processed inline as soon as they become ready.
    // ------------------------------------------------------------------------------------------------
    const int n = (int)nodes_.size();
    const int jobs = std::max(1, opts_.num_jobs);

    std::mutex tq_m;
    std::condition_variable tq_cv;
    std::deque<Task> task_q;
    std::mutex rq_m;
    std::condition_variable rq_cv;
    std::deque<TaskResult> result_q;
    std::atomic<bool> abort_flag{false}; ///< after an error, queued tasks are skipped (not run)
    std::atomic<bool> stop_flag{false};  ///< tells workers to exit once the task queue is drained

    auto push_task = [&](const Task& t) {
      { std::lock_guard<std::mutex> lk(tq_m); task_q.push_back(t); }
      tq_cv.notify_one();
    };
    auto pop_result = [&]() -> TaskResult {
      std::unique_lock<std::mutex> lk(rq_m);
      rq_cv.wait(lk, [&] { return !result_q.empty(); });
      TaskResult r = result_q.front();
      result_q.pop_front();
      return r;
    };

    // worker: pull tasks and launch the (blocking) subprocess; report the outcome
    auto worker = [&]() {
      while (true)
      {
        Task t;
        {
          std::unique_lock<std::mutex> lk(tq_m);
          tq_cv.wait(lk, [&] { return !task_q.empty() || stop_flag.load(); });
          if (task_q.empty()) break; // shutdown requested and nothing left to run
          t = task_q.front();
          task_q.pop_front();
        }
        TaskResult r;
        r.node = t.node;
        // an exception must never escape a std::thread callable (-> std::terminate); turn it into a result
        try
        {
          if (!abort_flag.load()) // after an error, queued tasks are skipped (just drained)
          {
            ExternalProcess proc;
            std::string perr;
            ExternalProcess::RETURNSTATE st = proc.run(t.exe, t.args, "", false, perr);
            if (st == ExternalProcess::RETURNSTATE::SUCCESS) { /* OK */ }
            else if (st == ExternalProcess::RETURNSTATE::FAILED_TO_START) { r.status = TOOL_NOT_FOUND; r.msg = "Could not start tool: " + perr; }
            // a TOPP tool missing its OWN external program exits 14 -> propagate as a skip (TOOL_NOT_FOUND)
            else if (childExitCode(perr) == 14) { r.status = TOOL_NOT_FOUND; r.msg = "Tool is missing a required external program: " + perr; }
            else { r.status = TOOL_FAILED; r.msg = "Tool failed: " + perr; }
          }
        }
        catch (const std::exception& e) { r.status = TOOL_FAILED; r.msg = std::string("internal error running tool: ") + e.what(); }
        catch (...) { r.status = TOOL_FAILED; r.msg = "internal error running tool"; }
        { std::lock_guard<std::mutex> lk(rq_m); result_q.push_back(r); }
        rq_cv.notify_one();
      }
    };
    std::vector<std::thread> pool;
    for (int i = 0; i < jobs; ++i) pool.emplace_back(worker);

    // Stop and join all workers. Allocation-free (only atomic stores + notify_all) so it cannot itself
    // throw: abort_flag makes still-queued tasks skip their subprocess; stop_flag makes each worker exit
    // once the queue is empty. This MUST run on every exit path (normal, error, or exception) -- otherwise
    // ~vector<thread> would destroy joinable threads and call std::terminate.
    auto shutdown = [&]() noexcept {
      abort_flag.store(true);
      stop_flag.store(true);
      tq_cv.notify_all();
      for (std::thread& w : pool) { if (w.joinable()) w.join(); }
    };

    Result error = OK;
    std::string emsg;
    // Everything between pool creation and shutdown() runs under this try so that ANY exception (e.g. a
    // std::bad_alloc) still joins the workers instead of letting ~vector<thread> call std::terminate.
    try
    {
    std::vector<int> in_degree(n, 0);
    for (const Edge& e : edges_) if (e.to >= 0 && e.from >= 0) ++in_degree[e.to];
    std::vector<int> remaining_tasks(n, 0);
    int outstanding = 0;

    // process a node whose inputs are all ready; dispatches tool rounds as tasks, runs instant nodes inline;
    // returns false (and sets error/emsg) on failure
    std::function<bool(int)> processReady = [&](int ni) -> bool {
      Node& node = nodes_[ni];
      try
      {
        switch (node.kind)
        {
          case NodeKind::INPUT:
            runInput_(node);
            break;
          case NodeKind::MERGER:
          {
            RoundPackages pkg;
            if (!buildRoundPackages(ni, pkg, emsg)) { error = INVALID; return false; }
            runMerger_(node, pkg);
            break;
          }
          case NodeKind::SPLITTER:
          {
            RoundPackages pkg;
            if (!buildRoundPackages(ni, pkg, emsg)) { error = INVALID; return false; }
            runSplitter_(node, pkg);
            break;
          }
          case NodeKind::OUTPUT:
          case NodeKind::OUTPUT_FOLDER:
            runOutputNode_(ni);
            break;
          case NodeKind::TOOL:
          {
            RoundPackages pkg;
            if (!buildRoundPackages(ni, pkg, emsg)) { error = INVALID; return false; }
            node.round_total = (int)pkg.size();
            std::vector<Task> tasks = prepareToolTasks_(ni, pkg);
            if (!tasks.empty())
            {
              remaining_tasks[ni] = (int)tasks.size();
              for (const Task& t : tasks) { push_task(t); ++outstanding; }
              return true; // not finished yet -- wait for the subprocesses to complete
            }
            break; // tool with zero rounds -> nothing to run
          }
        }
      }
      catch (const ExecError& e) { error = e.result; emsg = e.msg; return false; }
      catch (const Exception::BaseException& e) { error = TOOL_FAILED; emsg = e.what(); return false; }
      catch (const std::exception& e) { error = TOOL_FAILED; emsg = e.what(); return false; }

      // an instant node (or a tool with no rounds) is done -> mark finished and release its successors
      node.finished = true;
      for (int ei : node.out_edges)
      {
        int v = edges_[ei].to;
        if (--in_degree[v] == 0 && !processReady(v)) return false;
      }
      return true;
    };

    // kick off the source nodes (collected up-front so propagation never re-processes a node)
    std::vector<int> sources;
    for (int ni = 0; ni < n; ++ni) if (in_degree[ni] == 0) sources.push_back(ni);
    for (int ni : sources)
    {
      if (error != OK) break;
      if (!processReady(ni)) break;
    }

    // drive the tool subprocesses to completion, releasing successors as each tool node finishes
    while (outstanding > 0 && error == OK)
    {
      TaskResult r = pop_result();
      --outstanding;
      if (r.status != OK) { error = r.status; emsg = r.msg; break; }
      if (--remaining_tasks[r.node] == 0)
      {
        nodes_[r.node].finished = true;
        bool ok = true;
        for (int ei : nodes_[r.node].out_edges)
        {
          int v = edges_[ei].to;
          if (--in_degree[v] == 0 && !processReady(v)) { ok = false; break; }
        }
        if (!ok) break;
      }
    }
    } // end of the post-pool try block
    catch (const std::exception& e) { if (error == OK) { error = TOOL_FAILED; emsg = e.what(); } }
    catch (...) { if (error == OK) { error = TOOL_FAILED; emsg = "internal scheduler error"; } }

    shutdown(); // always joins every worker before 'pool' is destroyed (even on exception)

    if (error != OK) { error_msg = emsg; return error; }
    return OK;
  }

  std::vector<PipelineExecutor::Task> PipelineExecutor::prepareToolTasks_(int node_index, const RoundPackages& inputs)
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

    std::vector<Task> tasks;
    tasks.reserve(rounds);
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

      // build the task (the subprocess is launched later by a worker thread)
      Task t;
      t.node = node_index;
      t.exe = exe;
      if (!td.type.empty())
      {
        t.args.push_back("-type");
        t.args.push_back(td.type);
      }
      t.args.push_back("-ini");
      t.args.push_back(ini_file);
      tasks.push_back(std::move(t));
    }
    return tasks;
  }

  void PipelineExecutor::runOutputNode_(int node_index)
  {
    Node& node = nodes_[node_index];
    RoundPackages inputs = collectOutputNodeInputs_(node_index);
    node.round_total = (int)inputs.size();
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
