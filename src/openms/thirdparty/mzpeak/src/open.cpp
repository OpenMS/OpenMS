/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include "mzpeak/open.h"

#include <filesystem>
#include <memory>
#include <stdexcept>

#include "mzpeak/exception.h"
#include "mzpeak/io/directory.h"
#include "mzpeak/io/zip.h"

namespace MzPeak {

/******************************************************************************/
MzPeak::Index open(const fs::path& path, Validate validate)
{
  std::unique_ptr<MzPeak::IO::Archive> archive;

  if (fs::exists(path)) {
    if (fs::is_directory(path)) {
      archive = std::make_unique<MzPeak::IO::Directory>(path);
    } else {
      archive = std::make_unique<MzPeak::IO::Zip>(path);
    }
  } else {
    // FIXME:
    throw std::runtime_error("network access not implemented");
  }

  MzPeak::Index index(std::move(archive));

  if (validate == Validate::Checksums) {
    const ChecksumReport report = index.verify_checksums();
    if (!report.ok()) {
      // Name one member and both digests: "checksum mismatch" alone leaves the
      // caller unable to tell a damaged transfer from an edited archive.
      const ChecksumMismatch& first = report.mismatches.front();
      std::string msg("checksum mismatch in " + path.string() + ": " +
                      first.file_name + " hashes to " +
                      (first.actual.empty() ? std::string("nothing (member absent)")
                                            : first.actual) +
                      ", index records " + first.expected);
      if (report.mismatches.size() > 1) {
        msg += " (and " + std::to_string(report.mismatches.size() - 1) + " more)";
      }
      throw ChecksumError(msg);
    }
  }

  return index;
}

} // namespace MzPeak
