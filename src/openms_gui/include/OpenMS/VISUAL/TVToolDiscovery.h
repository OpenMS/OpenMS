// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>

#include <future>
#include <vector>
#include <unordered_map>

namespace OpenMS
{
    class Param;
    class String;

    class OPENMS_DLLAPI TVToolDiscovery
    {
    private:
      TVToolDiscovery() {}

      /// Return param for a given tool/util
      static Param getParamFromIni_(const String& tool_name);

      /// Contains a future param for each tool/util name
      static std::unordered_map<std::string, std::future<Param>> future_results_;

      /// Indicates whether creating all params has finished yet
      static bool params_ready_;

      /// Contains a mapping of each tool/util name to its param.
      static std::unordered_map<std::string, Param> params_;

    public:
      /// Start creating params for each tool/util by starting asynchronous threads.
      static void loadParams();

      /**
       * Returns a hash map (query is in O(1) on average) containing a param for each tool/util. Note that
       * it is possible that not all param futures have been finished yet if this function is called before waitForParams().
       * Therefore it is possible that this function has to wait.
       */
      static const std::unordered_map<std::string, Param>& getToolParams();

      /**
       * @brief Wait for all future params to finish evaluating.
       *
       * While waiting the GUI remains responsive. After waiting it is safe to access the params without further waiting.
       */
      static void waitForParams();
    };
}