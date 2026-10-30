// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineGraph.h>

namespace OpenMS
{
/**
  @brief Read and write TOPPAS workflow and resource files without Qt.

  Loading is transactional. Relative embedded inputs resolve against the workflow
  directory. Resource files replace all input lists, including missing keys.
  The Param overloads allow TOPPAS to exchange graph definitions without temporary files.
*/
class OPENMS_CLI_DLLAPI PipelineFile
{
public:
  void load(const std::string& filename, PipelineGraph& graph) const;
  void store(const std::string& filename, const PipelineGraph& graph) const;

  void loadParam(const Param& parameters, PipelineGraph& graph, const std::string& filename = "") const;
  Param storeParam(const PipelineGraph& graph, const std::string& filename = "") const;

  void loadResources(const std::string& filename, PipelineGraph& graph) const;
  void loadResourceParam(const Param& parameters, PipelineGraph& graph) const;
};
} // namespace OpenMS
