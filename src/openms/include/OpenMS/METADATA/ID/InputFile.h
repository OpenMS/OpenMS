// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/IDDataContainer.h>

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/METADATA/ID/MetaData.h>


#include <set>

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Information about input files that were processed
    struct InputFile
    {
      std::string name;

      std::string experimental_design_id;

      std::set<std::string> primary_files;

      explicit InputFile(const std::string& name,
                         const std::string& experimental_design_id = "",
                         const std::set<std::string>& primary_files =
                         std::set<std::string>()):
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

    using InputFiles = IDDataContainer<InputFile, std::string, std::string>;
    extern template class OPENMS_DLLAPI IDDataContainer<InputFile, std::string, std::string>;
    typedef IteratorWrapper<InputFiles::iterator> InputFileRef;

  }
}
