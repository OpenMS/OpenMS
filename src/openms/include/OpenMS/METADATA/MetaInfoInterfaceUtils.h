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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <algorithm>
#include <map>
#include <vector>


namespace OpenMS
{

  /**
    @brief Utilities operating on containers inheriting from MetaInfoInterface

    @ingroup MetaData
  */
  class /*OPENMS_DLLAPI -- disabled since it's template code only */ MetaInfoInterfaceUtils
  {
public:
    /// hide c'tors to avoid instantiation of utils class
    MetaInfoInterfaceUtils() = delete;
    MetaInfoInterfaceUtils(const MetaInfoInterfaceUtils&) = delete;
    MetaInfoInterfaceUtils& operator=(MetaInfoInterfaceUtils&) = delete;
    // no Move semantics for utils class

    ///@name Methods to find key sets
    //@{
    /**
      @brief Find keys in a collection of MetaInfoInterface objects which reach a certain frequency threshold.

      Searches the given iterator range for the keys of each element's MetaInfoInterface keys and returns those keys, which
      reach a certain frequency threshold. Common use cases 
      are @p min_frequency = 0 (i.e. take any key which occurs)
      and @p min_frequency = 100 (i.e. take only keys which are common to all elements in the iterator range).

      @tparam T_In Input container (e.g. std::vector or alike), containing objects which implement the MetaInfoInterface (i.e. support 'getKeys()')
      @tparam T_Out Output container of type T<String> (e.g. std::set<String>)
      @param start Iterator pointing to the initial position to search. (note: this does not need to correspond to the beginning of the container)
      @param end Iterator pointing to the end final position to search.
      @param min_frequency Minimum required frequency (in percent). Must be between 0-100. Other values are corrected to the closest value allowed.
      @return Returns a vector/list/set of keys passing the frequency criterion.
    */
    template<typename T_In, typename T_Out>
    static T_Out findCommonMetaKeys(const typename T_In::const_iterator& it_start, const typename T_In::const_iterator& it_end, float min_frequency)
    {
      // make sure min_frequency is within [0,100]
      min_frequency = std::min((float)100.0, std::max((float)0.0, min_frequency));

      std::map<String, uint> counter;
      typedef std::vector<String> KeysType;
      KeysType keys;
      for (typename T_In::const_iterator it = it_start; it != it_end; ++it)
      {
        it->getKeys(keys);
        for (KeysType::const_iterator itk = keys.begin(); itk != keys.end(); ++itk)
        {
          ++counter[*itk];
        }
      }
      // pick the keys which occur often enough
      const uint required_counts = uint(min_frequency / 100.0 * std::distance(it_start, it_end));
      T_Out common_keys;
      for (std::map<String, uint>::const_iterator it = counter.begin(); it != counter.end(); ++it)
      {
        if (it->second >= required_counts) 
        {
          common_keys.insert(common_keys.end(), it->first);
        }
      }
      return common_keys;
    }
  
  }; // class

} // namespace OPENMS

