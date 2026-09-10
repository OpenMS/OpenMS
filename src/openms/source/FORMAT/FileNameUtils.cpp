// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileNameUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/SYSTEM/PathUtils.h>

namespace OpenMS
{
  FileTypes::Type FileNameUtils::getTypeByFileName(const std::string& filename)
  {
    std::string basename = PathUtils::basename(filename), tmp;
    // special rules for "double extensions":
    if (StringUtils::hasSuffix(basename, ".pep.xml"))
    {
      return FileTypes::PEPXML;
    }
    if (StringUtils::hasSuffix(basename, ".prot.xml"))
    {
      return FileTypes::PROTXML;
    }
    if (StringUtils::hasSuffix(basename, ".xquest.xml"))
    {
      return FileTypes::XQUESTXML;
    }
    if (StringUtils::hasSuffix(basename, ".spec.xml"))
    {
      return FileTypes::SPECXML;
    }
    if (!StringUtils::has(basename, '.')) // no '.' => unknown type
    {
      // last chance, Bruker fid file
      if (basename == "fid")
      {
        return FileTypes::XMASS;
      }
      return FileTypes::UNKNOWN;
    }
    tmp = StringUtils::suffix(basename, '.');
    StringUtils::toUpper(tmp);
    if (tmp == "BZ2" || tmp == "GZ" || tmp == "ZIP")
    {
      // do not use getTypeByContent() here, as this is deadly for output files!
      return getTypeByFileName(StringUtils::prefix(filename, filename.size() - tmp.size() - 1)); // check name without compression suffix (e.g. bla.mzML.gz --> bla.mzML)
    }

    return FileTypes::nameToType(tmp);
  }

  bool FileNameUtils::hasValidExtension(const std::string& filename, const FileTypes::Type type)
  {
    FileTypes::Type ft = FileNameUtils::getTypeByFileName(filename);
    return (ft == type || ft == FileTypes::UNKNOWN);
  }

  std::string FileNameUtils::stripExtension(const std::string& filename)
  {
    if (!StringUtils::has(filename, '.'))
    {
      return filename;
    }
    // we don't just search for the last '.' and remove the suffix, because this could be wrong, e.g. bla.mzML.gz would become bla.mzML
    auto type = getTypeByFileName(filename);
    auto s_type = FileTypes::typeToName(type);
    size_t pos = StringUtils::toLowered(filename).rfind(StringUtils::toLowered(s_type)); // search backwards in entire string, because we could search for 'mzML' and have 'mzML.gz'
    if (pos == std::string::npos) // file type was FileTypes::UNKNOWN and we did not find '.unknown' as ending
    {
      size_t ext_pos = filename.rfind('.');
      size_t dir_sep = filename.find_last_of("/\\"); // look for '/' or '\'
      if (dir_sep != std::string::npos && dir_sep > ext_pos) // we found a directory separator after the last '.', e.g. '/my.dotted.dir/filename'! Ouch!
      { // do not strip anything, because there is no extension to strip
        return filename;
      }
      return StringUtils::prefix(filename, ext_pos);
    }
    return StringUtils::prefix(filename, pos - 1); // strip the '.' as well
  }

  std::string FileNameUtils::swapExtension(const std::string& filename, const FileTypes::Type new_type)
  {
    return stripExtension(filename) + "." + FileTypes::typeToName(new_type);
  }
}
