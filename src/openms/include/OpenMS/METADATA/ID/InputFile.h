// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include <set>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Information about input files that were processed
    struct InputFile
    {
      String name;

      String experimental_design_id;

      std::set<String> primary_files;

      explicit InputFile(const String& name,
                         const String& experimental_design_id = "",
                         const std::set<String>& primary_files =
                         std::set<String>()):
        name(name), experimental_design_id(experimental_design_id),
        primary_files(primary_files)
      {
      }

      InputFile(const InputFile& other) = default;

      /// Merge in data from another object
      InputFile& merge(const InputFile& other)
      {
        if (experimental_design_id.empty())
        {
          experimental_design_id = other.experimental_design_id;
        }
        else if (!other.experimental_design_id.empty() && experimental_design_id != other.experimental_design_id)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION, 
                                        "Trying to overwrite InputFile experimental design id with conflicting value.", 
                                        experimental_design_id);
        }
        primary_files.insert(other.primary_files.begin(),
                             other.primary_files.end());
        return *this;
      }
    };

    typedef boost::multi_index_container<
      InputFile,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          InputFile, String, &InputFile::name>>>
      > InputFiles;
    typedef IteratorWrapper<InputFiles::iterator> InputFileRef;

  }
}
