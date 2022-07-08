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
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
  /**
     @brief Scans for tools/utils and generates a param for each asynchronously.

     @details All tools and utils listed in the ToolHandler class are considered.

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
    TVToolDiscovery() :
      plugin_path_() {};

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
    void waitForPluginParams();

    /**
       @brief Returns a Param object containing the params for each tool/util.
       @details
       Note that it is possible that not all param futures have been finished (or loaded) yet if this function is called.
       In that case, the function starts param parsing (loadParam()) and waits for completion (waitForToolParams())
       before returning the result.
     */
    const Param& getToolParams();

    /**
       @brief Returns a param containing the params for each plugin.
       @details
       This function will ALWAYS trigger a reload of the plugins.
    */
    const Param& getPluginParams();

    /// Returns the list of read plugin names as saved in the ini.
    const std::vector<std::string>& getPlugins();

    /**
     * @brief Sets the path that will be searched for Plugins
     * @param path The new path to set
     * @param create Attempt to create the directory if it does not already exist
     * @returns False if setting/creating the path fails. True otherwise.
     */
    [[maybe_unused]] bool setPluginPath(const String& path, bool create=false);

    /// Returns the current set path to search plugins in
    const std::string getPluginPath();

    /// Returns the path to the plugin executable or an empty string if the plugin name is unknown
    std::string findPluginExecutable(const std::string& name);

  private:
    /** Returns param for a given tool/util. This function is thread-safe. Additionally inserts names of tools into 
        plugin list
     */
    static Param getParamFromIni_(const String& tool_path, bool plugins=false);

    /** Start creating params for each plugin in the set plugin path asynchronously
     *  This should only be called from waitForPluginParams() or the names in the plugins vector are not correct
     */
    void loadPluginParams();

    /// Returns a list of executables that are found at the plugin path
    const StringList getPlugins_();

    /// The filepath to search pugins in
    std::string plugin_path_;

    /// The futures for asyncronous loading of the tools/utils and plugins
    std::vector<std::future<Param>> tool_param_futures_;
    std::vector<std::future<Param>> plugin_param_futures_;

    /// Contains all the params of the tools/utils
    Param tool_params_;

    /// Contains all the params of the plugins
    Param plugin_params_;

    /// The names of all loaded plugins, this is used to add the plugins to the list in the ToolsDialog
    std::vector<std::string> plugins_;
  };
}