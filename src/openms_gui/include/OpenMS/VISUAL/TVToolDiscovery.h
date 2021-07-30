// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

namespace OpenMS
{
  /**
     @brief Scans for tools/utils and generates a param for each asynchronously.

     @details All tools and utils listed in the ToolHandler class are considered.

     @code
     TVToolDiscovery scanner;
     scanner.loadParams();
     // Do something else before explicitly waiting for the threads to finish
     ...
     // Wait when convenient. Keeps the GUI responsive while waiting
     scanner.waitForParams();
     // Access the params. If no special timing for waiting or loading is needed this function can be safely called directly.
     scanner.getToolParams();
     @endcode
   */
  class OPENMS_GUI_DLLAPI TVToolDiscovery
  {
  public:
    TVToolDiscovery() {};

    TVToolDiscovery(const TVToolDiscovery &) = delete;

    TVToolDiscovery &operator=(const TVToolDiscovery &) = delete;

    ~TVToolDiscovery() {};

    /// Start creating params for each tool/util asynchronously
    void loadParams();

    /**
       @brief Wait for all future params to finish evaluating.
       @details
       While waiting the GUI remains responsive. After waiting it is safe to access the params without further waiting.
     */
    void waitForParams();

    /**
       @brief Returns a hash map containing a param for each tool/util.
       @details
       Note that it is possible that not all param futures have been finished (or loaded) yet if this function is called.
       In that case, the function starts param parsing (loadParam()) and waits for completion (waitForParams())
       before returning the result.
     */
    const std::map<std::string, Param> &getToolParams();

  private:
    /// Returns param for a given tool/util. This function is thread-safe
    static Param getParamFromIni_(const std::string &tool_name);

    /// Contains a param future for each tool/util name
    std::map<std::string, std::future<Param>> param_futures_;

    /// Contains a mapping of each tool/util name to its param.
    std::map<std::string, Param> params_;
  };
}