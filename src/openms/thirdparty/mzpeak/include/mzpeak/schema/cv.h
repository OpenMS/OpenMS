/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <optional>
#include <string>

namespace MzPeak::Schema {

/**
 * A controlled vocabulary term
 */
class CV {
public:
  /// Constructor.
  CV(std::string_view code, std::string_view accession)
      : code_(code)
      , accession_(accession)
  {
  }

  /// Destructor.
  ~CV() = default;

  /// Parse a string like "MS:1000511"
  static std::optional<CV> from_string(std::string_view);

  /// Convert this CV term to a string like "MS:1000511"
  std::string to_string() const;

  /// Access the CV vendor code.
  const std::string& code() const { return code_; }

  /// Access the accession number.
  const std::string& accession() const { return accession_; }

  /// Test for equality.
  bool operator==(const CV& other) const
  {
    return code_ == other.code_ && accession_ == other.accession_;
  }

  /// Less than.
  bool operator<(const CV& other) const
  {
    return code_ < other.code_ && accession_ < other.accession_;
  }

private:
  std::string code_;
  std::string accession_;
};

} // namespace MzPeak::Schema
