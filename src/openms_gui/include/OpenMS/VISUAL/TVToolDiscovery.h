// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: David Voigt $
// $Authors: David Voigt $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <future>
#include <map>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
  /**
     @brief Scans for tools and generates a param for each asynchronously.

     @details All tools listed in the ToolHandler class are considered.

     @code
     TVToolDiscovery scanner;
     scanner.loadToolParams();
     // Do something else before explicitly waiting for the threads to finish
     ...
     // Wait when convenient. Keeps the GUI responsive while waiting
     scanner.waitForToolParams();
     // Access the params. If no special timing for waiting or loading is needed this function can be safely called directly.
     scanner.getToolParams();
     @endcode
   */
  class OPENMS_GUI_DLLAPI TVToolDiscovery
  {
  public:
    TVToolDiscovery() = default;

    TVToolDiscovery(const TVToolDiscovery&) = delete;

    TVToolDiscovery& operator=(const TVToolDiscovery&) = delete;

    ~TVToolDiscovery() = default;

    /// Start creating params for each tool/util asynchronously
    void loadToolParams();

    /**
       @brief Wait for all future params to finish evaluating.
       @details
       While waiting the GUI remains responsive. After waiting it is safe to access the params without further waiting.
     */
    void waitForToolParams();

    /**
       @brief Returns a Param object containing the params for each tool/util.
       @details
       Note that it is possible that not all param futures have been finished (or loaded) yet if this function is called.
       In that case, the function starts param parsing (loadParam()) and waits for completion (waitForToolParams())
       before returning the result.
     */
    const Param& getToolParams();

    /// set the verbosity level of the tool discovery for debug purposes
    void setVerbose(int verbosity_level);

  private:
    /// Returns param for a given tool/util. This function is thread-safe.
    static Param getParamFromIni_(const std::string& tool_path);

    /// The futures for asyncronous loading of the tools
    std::vector<std::future<Param>> tool_param_futures_;

    /// Contains all the params of the tools/utils
    Param tool_params_;

    /// Set to value > 0 to output tool discovery debug information
    int verbosity_level_ = 0;
  };
}