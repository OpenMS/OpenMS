// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileNameUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

namespace OpenMS
{
  namespace
  {
    /// Compression suffixes we see through: 'bla.mzML.gz' is an mzML, not a distinct format.
    bool isCompressionSuffix(const std::string& suffix_upper)
    {
      return suffix_upper == "GZ" || suffix_upper == "BZ2" || suffix_upper == "ZIP";
    }
  }

  FileTypes::Type FileNameUtils::matchExtension_(const std::string& filename, size_t& ext_start)
  {
    ext_start = std::string::npos;
    // npos + 1 wraps to 0, i.e. the whole string is the basename when there is no separator
    const size_t name_start = filename.find_last_of("\\/") + 1;

    if (filename.find('.', name_start) == std::string::npos) // no '.' in the basename => no extension to match
    {
      // last chance, Bruker fid file
      return StringUtils::substr(filename, name_start) == "fid" ? FileTypes::XMASS : FileTypes::UNKNOWN;
    }

    const size_t last_dot = filename.rfind('.'); // >= name_start, since the basename has a dot
    if (isCompressionSuffix(StringUtils::toUppered(StringUtils::substr(filename, last_dot + 1))))
    {
      // check the name without the compression suffix (e.g. bla.mzML.gz --> bla.mzML) and report the
      // whole '.mzML.gz' span as the extension. Do not use getTypeByContent() here, as this is deadly for output files!
      const FileTypes::Type inner = matchExtension_(StringUtils::prefix(filename, last_dot), ext_start);
      if (inner != FileTypes::UNKNOWN)
      {
        // usually ext_start already points at the inner extension. The exception is the extensionless
        // Bruker 'fid': for 'fid.gz' the compression suffix is the only extension-shaped span there is.
        if (ext_start == std::string::npos) ext_start = last_dot;
        return inner;
      }
      ext_start = last_dot;                          // e.g. 'archive.gz' => only '.gz' is an extension
      return FileTypes::UNKNOWN;
    }

    // Try every dot-delimited suffix of the basename. find() walks left to right and earlier dots yield
    // longer suffixes, so the first hit is the longest match, e.g. '.pep.xml' wins over '.xml'.
    for (size_t dot = filename.find('.', name_start); dot != std::string::npos; dot = filename.find('.', dot + 1))
    {
      const FileTypes::Type type = FileTypes::nameToType(StringUtils::substr(filename, dot + 1));
      if (type != FileTypes::UNKNOWN)
      {
        ext_start = dot;
        return type;
      }
    }

    ext_start = last_dot; // unknown format, but the last dot still delimits something extension-shaped
    return FileTypes::UNKNOWN;
  }

  FileTypes::Type FileNameUtils::getTypeByFileName(const std::string& filename)
  {
    size_t ext_start;
    return matchExtension_(filename, ext_start);
  }

  FileTypes::Type FileNameUtils::compressionType(const std::string& filename)
  {
    const size_t last_dot = filename.rfind('.');
    if (last_dot == std::string::npos) return FileTypes::UNKNOWN;
    const size_t name_start = filename.find_last_of("\\/") + 1;
    if (last_dot < name_start) return FileTypes::UNKNOWN; // the dot belongs to a directory, not to the basename
    const std::string suffix = StringUtils::toUppered(StringUtils::substr(filename, last_dot + 1));
    if (suffix == "GZ") return FileTypes::GZ;
    if (suffix == "BZ2") return FileTypes::BZ2;
    if (suffix == "ZIP") return FileTypes::ZIP;
    return FileTypes::UNKNOWN;
  }

  bool FileNameUtils::hasCompressionSuffix(const std::string& filename)
  {
    return compressionType(filename) != FileTypes::UNKNOWN;
  }

  bool FileNameUtils::hasValidExtension(const std::string& filename, const FileTypes::Type type)
  {
    FileTypes::Type ft = FileNameUtils::getTypeByFileName(filename);
    return (ft == type || ft == FileTypes::UNKNOWN);
  }

  std::string FileNameUtils::stripExtension(const std::string& filename)
  {
    // we don't just search for the last '.' and remove the suffix, because this could be wrong:
    // 'bla.mzML.gz' must lose both suffixes and 'bla.pep.xml' must lose the compound extension as a whole
    size_t ext_start;
    matchExtension_(filename, ext_start);
    if (ext_start == std::string::npos) // nothing extension-shaped, e.g. '/my.dotted.dir/filename'
    {
      return filename;
    }
    return StringUtils::prefix(filename, ext_start); // strip the '.' as well
  }

  std::string FileNameUtils::swapExtension(const std::string& filename, const FileTypes::Type new_type)
  {
    return stripExtension(filename) + "." + FileTypes::typeToName(new_type);
  }
}
