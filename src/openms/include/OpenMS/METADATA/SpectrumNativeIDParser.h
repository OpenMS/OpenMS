// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <boost/regex.hpp>

namespace OpenMS
{
  /**
    @brief Parser for extracting scan numbers from spectrum native IDs

    This class provides static functions for parsing native ID strings from various
    mass spectrometry file formats. Native IDs are vendor-specific identifiers that
    encode information such as scan numbers, file indices, and experiment numbers.

    @par Supported Native ID Formats

    The parser supports the following native ID formats based on CV (Controlled Vocabulary)
    accessions from the PSI-MS ontology:

    | CV Accession | Vendor/Format | Native ID Pattern | Example |
    |--------------|---------------|-------------------|---------|
    | MS:1000768 | Thermo | scan=NUMBER | scan=42 |
    | MS:1000769 | Waters | function=X process=Y scan=NUMBER | function=2 process=1 scan=100 |
    | MS:1000770 | WIFF (AB Sciex) | sample=X period=Y cycle=Z experiment=W | sample=1 period=1 cycle=42 experiment=1 |
    | MS:1000771 | Bruker/Agilent | scan=NUMBER | scan=42 |
    | MS:1000772 | Bruker BAF | scan=NUMBER | scan=42 |
    | MS:1000773 | Bruker FID | file=NUMBER | file=42 |
    | MS:1000774 | Index-based | index=NUMBER | index=42 (returns 43) |
    | MS:1000775 | Single peak list | file=NUMBER | file=42 |
    | MS:1000776 | Thermo/Bruker TDF | scan=NUMBER | scan=42 |
    | MS:1000777 | Generic spectrum | spectrum=NUMBER | spectrum=42 |
    | MS:1001508 | Agilent MassHunter | scanId=NUMBER | scanId=42 |
    | MS:1001530 | Numeric-only | NUMBER | 42 |

    @par Usage Examples

    @code
    // Extract scan number using CV accession
    Int scan = SpectrumNativeIDParser::extractScanNumber("scan=42", "MS:1000768");  // returns 42

    // Get regex pattern from native ID format
    std::string regex = SpectrumNativeIDParser::getRegExFromNativeID("scan=123");  // returns "scan=(?<GROUP>\d+)"

    // Check if string is a native ID
    bool is_native = SpectrumNativeIDParser::isNativeID("scan=123");  // returns true
    @endcode

    @see SpectrumLookup
  */
  class OPENMS_DLLAPI SpectrumNativeIDParser
  {
  public:
    /**
       @brief Extract the scan number from the native ID of a spectrum using a regular expression

       @param[in] native_id Spectrum native ID string
       @param[in] scan_regexp Regular expression to use (must contain the named group "?<SCAN>")
       @param[in] no_error Suppress the exception on failure and return -1 instead

       @throw Exception::ParseError if the scan number could not be extracted (unless @p no_error is set)

       @return Scan number of the spectrum (or -1 on failure to extract)

       @note The regular expression must contain a capture group, and the last matching
             subgroup is used as the scan number.
    */
    static Int extractScanNumber(const std::string& native_id,
                                 const boost::regex& scan_regexp,
                                 bool no_error = false);

    /**
       @brief Extract the scan number from the native ID using a CV accession

       @param[in] native_id Spectrum native ID string
       @param[in] native_id_type_accession CV accession specifying the native ID format
             (e.g., "MS:1000768" for Thermo, "MS:1000770" for WIFF)

       @return Scan number of the spectrum (or -1 on failure to extract)

       @note For WIFF files (MS:1000770), the return value is computed as cycle * 1000 + experiment.
       @note For index-based native IDs (MS:1000774), the return value is index + 1 for pepXML compatibility.
    */
    static Int extractScanNumber(const std::string& native_id,
                                 const std::string& native_id_type_accession);

    /**
       @brief Determine the regular expression to extract scan/index numbers from native IDs

       @param[in] native_id A native ID string to analyze

       @return Regular expression string with named group `(?<GROUP>\d+)` that matches
               the scan or index number in the native ID

       This function examines the prefix of the native ID to determine the appropriate
       regular expression pattern:
       - `scan=`, `controllerType=`, `function=` → `scan=(?<GROUP>\d+)`
       - `index=` → `index=(?<GROUP>\d+)`
       - `scanId=`, `scanID=` → `scanId=(?<GROUP>\d+)` or `scanID=(?<GROUP>\d+)`
       - `spectrum=` → `spectrum=(?<GROUP>\d+)`
       - `file=` → `file=(?<GROUP>\d+)`
       - Plain number → `(?<GROUP>\d+)`
    */
    static std::string getRegExFromNativeID(const std::string& native_id);

    /**
       @brief Check if a spectrum identifier is a native ID from a vendor file

       @param[in] id Spectrum identifier string to check

       @return True if the string matches a known native ID prefix pattern

       Recognized prefixes: scan=, scanId=, scanID=, controllerType=, function=, sample=, index=, spectrum=, file=
    */
    static bool isNativeID(const std::string& id);

  };

} // namespace OpenMS
