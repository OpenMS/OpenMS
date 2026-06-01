// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
