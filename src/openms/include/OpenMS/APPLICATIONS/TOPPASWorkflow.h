// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Qt-free in-memory representation of a TOPPAS workflow (.toppas) file.

    TOPPASWorkflow holds the serialized content of a @em .toppas file -- the header
    (version + description), the list of nodes (vertices) and the list of edges -- and
    provides (de)serialization via ParamXMLFile (the on-disk @em .toppas format).

    It is deliberately independent of Qt so that the workflow document can be read,
    written and inspected without the GUI. The Qt-based @c TOPPASScene uses this class
    for loading and storing, mapping between these plain structs and its visual
    QGraphicsItem vertices/edges.

    Conventions kept identical to the historical TOPPASScene serialization:
    - Node ids are stored verbatim as strings (TOPPASScene uses the topological index).
    - Input file names are stored verbatim (typically relative to the @em .toppas file);
      absolutization/relativization relative to the file location is the caller's job.
    - Edge source/target parameters are stored as strings: a parameter name, the
      sentinel @c "__no_name__", or -- for very old (pre-1.9) files -- a numeric index.
    - store()/toParam() always write the @em current %OpenMS version into @c info:version
      (matching the previous behavior).

    @ingroup Applications
  */
  class OPENMS_DLLAPI TOPPASWorkflow
  {
  public:
    /// Sentinel used for an edge endpoint that is not bound to a named parameter.
    static const std::string NO_NAME;

    /// Type of a workflow node (mirrors the @c toppas_type field).
    enum class NodeType
    {
      INPUT_FILE_LIST, ///< a list of input files fed into the pipeline
      OUTPUT_FILE_LIST, ///< collects output files
      OUTPUT_FOLDER,   ///< collects output into a folder
      TOOL,            ///< a TOPP tool invocation
      MERGER,          ///< merges incoming file lists
      SPLITTER,        ///< splits an incoming file list
      UNKNOWN          ///< unrecognized type (forward compatibility)
    };

    /// A single workflow node (vertex).
    struct Node
    {
      std::string id;                ///< node id as stored in the file (TOPPASScene: topological index)
      NodeType type = NodeType::UNKNOWN;
      double x_pos = 0.0;            ///< canvas x position
      double y_pos = 0.0;           ///< canvas y position
      bool recycle_output = false;  ///< output-recycling flag (since TOPPAS 1.9)

      // NodeType::TOOL
      std::string tool_name;        ///< name of the TOPP tool
      std::string tool_type;        ///< tool sub-type (algorithm), may be empty
      Param parameters;             ///< the tool's parameters (without any key prefix)

      // NodeType::INPUT_FILE_LIST
      std::vector<std::string> file_names; ///< verbatim file names as stored (usually relative)

      // NodeType::OUTPUT_FILE_LIST / NodeType::OUTPUT_FOLDER
      std::string output_folder_name; ///< custom output folder name (may be empty)

      // NodeType::MERGER
      bool round_based = true;      ///< round-based vs. wait-for-all merge mode (default: round-based)
    };

    /// A directed edge connecting an output of @c source_id to an input of @c target_id.
    struct Edge
    {
      std::string source_id;        ///< id of the source node
      std::string target_id;        ///< id of the target node
      std::string source_out_param; ///< source output parameter name (or NO_NAME, or a numeric index for very old files)
      std::string target_in_param;  ///< target input parameter name (or NO_NAME, or a numeric index for very old files)
    };

    /// %OpenMS version string that wrote the file (empty for pre-1.9 files lacking the tag).
    std::string version;
    /// Free-text workflow description (CDATA wrapper removed).
    std::string description;
    /// All nodes, in file order.
    std::vector<Node> nodes;
    /// All edges, in file order.
    std::vector<Edge> edges;

    /// Convert a @c toppas_type string to a NodeType (UNKNOWN if unrecognized).
    static NodeType typeFromString(const std::string& s);
    /// Convert a NodeType to its @c toppas_type string (empty for UNKNOWN).
    static std::string typeToString(NodeType t);

    /// Fill this workflow from an already parsed @em .toppas Param tree.
    void fromParam(const Param& param);
    /// Serialize this workflow into a Param tree in the @em .toppas layout.
    Param toParam() const;

    /// Load a @em .toppas file (ParamXMLFile format). Throws on file/parse errors.
    void load(const std::string& filename);
    /// Store this workflow to a @em .toppas file (ParamXMLFile format).
    void store(const std::string& filename) const;
  };

} // namespace OpenMS
