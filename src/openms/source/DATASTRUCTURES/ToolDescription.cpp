// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/config.h>

using namespace std;

namespace OpenMS
{

  namespace Internal
  {

    FileMapping & FileMapping::operator=(const FileMapping & rhs)
    {
      if (this == &rhs) return *this;

      location = rhs.location;
      target = rhs.target;
      return *this;
    }

    OPENMS_DLLAPI MappingParam & MappingParam::operator=(const MappingParam & rhs)
    {
      if (this == &rhs) return *this;

      mapping = rhs.mapping;
      pre_moves = rhs.pre_moves;
      post_moves = rhs.post_moves;
      return *this;
    }

    // default C'Tor
    ToolDescriptionInternal::ToolDescriptionInternal() :
      is_internal(false),
      name(),
      category(),
      types()
    {
    }

    // C'Tor with arguments
    ToolDescriptionInternal::ToolDescriptionInternal(const bool p_is_internal, const String & p_name, const String & p_category, const StringList & p_types) :
      is_internal(p_is_internal),
      name(p_name),
      category(p_category),
      types(p_types)
    {
    }

    ToolDescriptionInternal::ToolDescriptionInternal(const String & p_name, const StringList & p_types) :
      is_internal(false),
      name(p_name),
      category(),
      types(p_types)
    {
    }

    ToolDescriptionInternal & ToolDescriptionInternal::operator=(const ToolDescriptionInternal & rhs)
    {
      if (this == &rhs)
        return *this;

      is_internal = rhs.is_internal;
      name = rhs.name;
      category = rhs.category;
      types = rhs.types;
      return *this;
    }

    bool ToolDescriptionInternal::operator==(const ToolDescriptionInternal & rhs) const
    {
      if (this == &rhs)
        return true;

      return is_internal == rhs.is_internal
             && name == rhs.name
             && category == rhs.category
             && types == rhs.types;
    }

    bool ToolDescriptionInternal::operator<(const ToolDescriptionInternal & rhs) const
    {
      if (this == &rhs)
        return false;

      return name + "." + ListUtils::concatenate(types, ",") < rhs.name + "." + ListUtils::concatenate(rhs.types, ",");
    }

    // default CTor
    ToolDescription::ToolDescription() :
      external_details()
    {
    }

    // C'Tor for internal TOPP tools
    ToolDescription::ToolDescription(const String & p_name, const String & p_category, const StringList & p_types) :
      ToolDescriptionInternal(true, p_name, p_category, p_types)
    {
    }

    void ToolDescription::addExternalType(const String & type, const ToolExternalDetails & details)
    {
      types.push_back(type);
      external_details.push_back(details);
    }

    void ToolDescription::append(const ToolDescription & other)
    {
      // sanity check
      if (is_internal != other.is_internal
         || name != other.name
          //|| category != other.category
         || (is_internal && external_details.size() > 0)
         || (other.is_internal && other.external_details.size() > 0)
         || (!is_internal && external_details.size() != types.size())
         || (!other.is_internal && other.external_details.size() != other.types.size())
          )
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Extending (external) ToolDescription failed!", "");
      }

      // append types and external information
      types.insert(types.end(), other.types.begin(), other.types.end());
      external_details.insert(external_details.end(), other.external_details.begin(), other.external_details.end());

      // check that types are unique
      std::set<String> unique_check;
      unique_check.insert(types.begin(), types.end());
      if (unique_check.size() != types.size())
      {
        LOG_ERROR << "A type appears at least twice for the TOPP/UTIL '" << name << "'. Types given are '" << ListUtils::concatenate(types, ", ") << "'\n";
        if (name == "GenericWrapper")
        {
          LOG_ERROR << "Check the .ttd files in your share/ folder and remove duplicate types!\n";
        }
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "see above!", "");
      }
    }

    ToolDescription & ToolDescription::operator=(const ToolDescription & rhs)
    {
      if (this == &rhs)
        return *this;

      ToolDescriptionInternal::operator=(rhs);
      external_details = rhs.external_details;
      return *this;
    }

  }

  Internal::ToolDescription bla;

} // namespace OpenMS
