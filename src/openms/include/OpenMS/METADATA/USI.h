// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/config.h>

#include <optional>

namespace OpenMS
{
  /**
    @brief Utility class for handling Universal Spectrum Identifiers (USI).

    The Universal Spectrum Identifier (USI) is a standardized identifier format defined by the
    HUPO-PSI (Human Proteome Organization - Proteomics Standards Initiative) for uniquely
    referencing mass spectrometry spectra across public repositories. The USI format enables
    unambiguous identification of spectra in ProteomeXchange partner repositories (PRIDE, 
    PeptideAtlas, MassIVE, jPOST, iProX) and spectral libraries.

    @par USI Format
    The general format of a USI is:
    @code
    mzspec:<collection>:<ms_run>:<index_type>:<index>[:interpretation]
    @endcode

    Where:
    - @b collection: ProteomeXchange dataset identifier (e.g., PXD000561) or spectral library name
    - @b ms_run: Name of the MS run file (e.g., sample1.mzML)
    - @b index_type: Type of spectrum index, typically "scan" or "index"
    - @b index: Numeric identifier of the spectrum
    - @b interpretation: Optional peptide interpretation in ProForma format (sequence/charge)

    @par Examples
    - Basic spectrum reference: @c mzspec:PXD000561:Adult_Frontalcortex.mzML:scan:12345
    - With peptide interpretation: @c mzspec:PXD000561:Adult_Frontalcortex.mzML:scan:12345:PEPTIDEK/2
    - With modification: @c mzspec:PXD000561:sample.mzML:scan:1234:PEPT[Phospho]IDEK/2

    @par PSI-MS CV Term
    USI is defined in the PSI-MS ontology as MS:1003063 (universal spectrum identifier).

    @see SpectrumLookup, SpectrumSettings, PeptideIdentification
    @see https://www.psidev.info/usi (USI specification)

    @ingroup Metadata
  */
  class OPENMS_DLLAPI USI
  {
  public:
    /// Supported index types in USI
    enum class IndexType
    {
      SCAN,    ///< scan number (most common)
      INDEX,   ///< zero-based spectrum index
      NATIVEID ///< native spectrum identifier from the instrument
    };

    /// @name Constructors and Destructor
    //@{
    /// Default constructor (creates an empty/invalid USI)
    USI();

    /**
      @brief Construct a USI from its components.

      @param collection ProteomeXchange dataset identifier or spectral library name
      @param ms_run Name of the MS run file
      @param index_type Type of spectrum index (scan, index, or nativeId)
      @param index Numeric identifier of the spectrum
      @param interpretation Optional peptide interpretation (sequence/charge)
    */
    USI(const std::string& collection,
        const std::string& ms_run,
        IndexType index_type,
        const std::string& index,
        const std::string& interpretation = "");

    /**
      @brief Construct a USI from its string representation.

      @param usi_string Complete USI string to parse
      @throw Exception::ParseError if the USI string is malformed
    */
    explicit USI(const std::string& usi_string);

    /// Copy constructor
    USI(const USI&) = default;

    /// Move constructor
    USI(USI&&) noexcept = default;

    /// Destructor
    ~USI() = default;
    //@}

    /// @name Assignment operators
    //@{
    /// Copy assignment operator
    USI& operator=(const USI&) = default;

    /// Move assignment operator
    USI& operator=(USI&&) noexcept = default;
    //@}

    /// @name Comparison operators
    //@{
    /// Equality operator
    bool operator==(const USI& rhs) const;

    /// Inequality operator
    bool operator!=(const USI& rhs) const;

    /// Less-than operator (for sorting/containers)
    bool operator<(const USI& rhs) const;
    //@}

    /// @name Validity checking
    //@{
    /**
      @brief Check if this USI is valid (has all required components).

      A valid USI must have non-empty collection, ms_run, and index fields.

      @return True if the USI is valid, false otherwise
    */
    bool isValid() const;

    /**
      @brief Check if a string is a valid USI format.

      Validates that the string follows the USI specification format.

      @param usi_string std::string to validate
      @return True if the string is a valid USI, false otherwise
    */
    static bool isValidUSI(const std::string& usi_string);

    /**
      @brief Try to parse a USI string (non-throwing).

      This is the preferred method when you want to validate and use a USI
      in a single operation, avoiding double-parsing.

      @param usi_string std::string to parse
      @return Parsed USI if valid, std::nullopt otherwise

      @code
      // Prefer this pattern:
      if (auto usi = USI::tryParse(str)) {
        use(*usi);
      }
      // Over this (which parses twice):
      if (USI::isValidUSI(str)) {
        USI usi(str);
        use(usi);
      }
      @endcode
    */
    static std::optional<USI> tryParse(const std::string& usi_string);
    //@}

    /// @name Accessors
    //@{
    /// Returns the collection (dataset identifier or library name)
    const std::string& getCollection() const;

    /// Sets the collection (dataset identifier or library name)
    void setCollection(const std::string& collection);

    /// Returns the MS run file name
    const std::string& getMSRun() const;

    /// Sets the MS run file name
    void setMSRun(const std::string& ms_run);

    /// Returns the index type
    IndexType getIndexType() const;

    /// Sets the index type
    void setIndexType(IndexType index_type);

    /// Returns the spectrum index as a string
    const std::string& getIndex() const;

    /// Sets the spectrum index
    void setIndex(const std::string& index);

    /// Returns the interpretation (peptide sequence/charge)
    const std::string& getInterpretation() const;

    /// Sets the interpretation (peptide sequence/charge)
    void setInterpretation(const std::string& interpretation);

    /// Returns true if an interpretation is present
    bool hasInterpretation() const;
    //@}

    /// @name std::string conversion
    //@{
    /**
      @brief Converts the USI to its string representation.

      @return Complete USI string
    */
    std::string toString() const;

    /**
      @brief Parse a USI string and populate this object.

      @param usi_string USI string to parse
      @return True if parsing succeeded, false otherwise
    */
    bool fromString(const std::string& usi_string);
    //@}

    /// @name Static utility methods
    //@{
    /**
      @brief Create a USI from spectrum metadata.

      This is a convenience method for creating a USI from commonly available
      spectrum information.

      @param dataset_id ProteomeXchange dataset identifier (e.g., "PXD000561")
      @param filename MS file name (e.g., "sample.mzML")
      @param scan_number Scan number
      @param interpretation Optional ProForma interpretation (e.g., "PEPTIDEK/2")
      @return Constructed USI object
    */
    static USI createFromScanNumber(const std::string& dataset_id,
                                    const std::string& filename,
                                    int scan_number,
                                    const std::string& interpretation = "");

    /**
      @brief Create a USI from a native spectrum identifier.

      @param dataset_id ProteomeXchange dataset identifier
      @param filename MS file name
      @param native_id Native spectrum identifier from the instrument
      @return Constructed USI object
    */
    static USI createFromNativeID(const std::string& dataset_id,
                                  const std::string& filename,
                                  const std::string& native_id);

    /**
      @brief Extract scan number from a native spectrum ID.

      Attempts to extract a numeric scan number from various native ID formats.
      This is useful when converting native IDs to USI format.

      @param native_id Native spectrum identifier
      @return Scan number if extraction succeeded, std::nullopt otherwise
    */
    static std::optional<int> extractScanNumberFromNativeID(const std::string& native_id);

    /**
      @brief Get the index type string representation.

      @param index_type Index type enum value
      @return std::string representation ("scan", "index", or "nativeId")
    */
    static std::string indexTypeToString(IndexType index_type);

    /**
      @brief Parse index type from string.

      @param type_string std::string representation of index type
      @return Index type enum value
      @throw Exception::InvalidValue if the string is not a valid index type
    */
    static IndexType indexTypeFromString(const std::string& type_string);

    /**
      @brief Extract basename from a file path (removes directory path).

      This is useful when converting full file paths to MS run names for USI.
      Example: "/path/to/sample.mzML" -> "sample.mzML"
      Example: "file:///C:/data/sample.mzML" -> "sample.mzML"

      @param filepath Full file path or URI
      @return Basename without path component
    */
    static std::string extractBasename(const std::string& filepath);
    //@}

    /// @name CV term information
    //@{
    /// Returns the PSI-MS CV accession for USI (MS:1003063)
    static const std::string& getCVAccession();

    /// Returns the PSI-MS CV name for USI
    static const std::string& getCVName();
    //@}

  protected:
    std::string collection_;      ///< Dataset identifier or library name
    std::string ms_run_;          ///< MS run file name
    IndexType index_type_;   ///< Type of spectrum index
    std::string index_;           ///< Spectrum index value
    std::string interpretation_;  ///< Optional peptide interpretation

    /// Prefix for all USIs
    static const std::string USI_PREFIX;

    /// CV accession for USI
    static const std::string CV_ACCESSION;

    /// CV name for USI
    static const std::string CV_NAME;
  };

  /// Output stream operator
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const USI& usi);

} // namespace OpenMS
