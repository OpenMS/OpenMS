// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/OpenMSConfig.h>

#include <iosfwd>
#include <optional>
#include <utility>
#include <vector>

namespace OpenMS
{

  //--------------------------------------------------------------------------
  // Data Structures
  //--------------------------------------------------------------------------

  /**
    @brief Ion series types for mzPAF peak annotations

    Identifies the type of fragment ion or special ion in an mzPAF annotation.

    @ingroup Chemistry
  */
  enum class MzPAFIonSeries
  {
    A,          ///< a-ion (N-terminal, loses CO)
    B,          ///< b-ion (N-terminal)
    C,          ///< c-ion (N-terminal, ETD)
    X,          ///< x-ion (C-terminal)
    Y,          ///< y-ion (C-terminal)
    Z,          ///< z-ion (C-terminal, ETD)
    PRECURSOR,  ///< Precursor ion (p)
    IMMONIUM,   ///< Immonium ion (I)
    INTERNAL,   ///< Internal fragment (m)
    REPORTER,   ///< Reporter ion (r)
    FORMULA,    ///< Chemical formula ion (f)
    NAMED,      ///< Named compound (_)
    UNKNOWN     ///< Unknown or unrecognized ion type
  };

  /**
    @brief Unit for mass delta values in mzPAF annotations

    @ingroup Chemistry
  */
  enum class MzPAFDeltaUnit
  {
    DALTON,  ///< Mass delta in Daltons (Da)
    PPM      ///< Mass delta in parts per million
  };

  /**
    @brief Neutral loss in an mzPAF annotation

    Represents a neutral loss like -H2O or -NH3.

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI MzPAFNeutralLoss
  {
    EmpiricalFormula formula;   ///< Parsed chemical formula of the loss

    bool operator==(const MzPAFNeutralLoss& other) const;
    bool operator!=(const MzPAFNeutralLoss& other) const { return !(*this == other); }
  };

  /**
    @brief Mass delta in an mzPAF annotation

    Represents a mass difference with optional unit (Da or ppm).
    Example: /0.001 (Da) or /-1.4ppm

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI MzPAFMassDelta
  {
    double value = 0.0;         ///< Mass delta value
    MzPAFDeltaUnit unit = MzPAFDeltaUnit::DALTON; ///< Unit (DALTON or PPM)

    bool operator==(const MzPAFMassDelta& other) const;
    bool operator!=(const MzPAFMassDelta& other) const { return !(*this == other); }
  };

  /**
    @brief A single mzPAF peak annotation

    Represents one annotation for a peak in mzPAF format.
    mzPAF is the HUPO-PSI standard string notation for fragment ion peak annotations.

    Examples:
    - y4                  - Simple y-ion at position 4
    - b2-H2O              - b-ion with neutral loss
    - y4^2                - Doubly charged y-ion
    - y2+2i               - Second isotope peak
    - y4/0.001            - With mass delta in Da
    - y4/-1.4ppm          - With mass delta in ppm
    - y4*0.75             - With confidence score
    - IY                  - Immonium ion (tyrosine)
    - m3:6                - Internal fragment (positions 3-6)
    - r[TMT127N]          - Reporter ion
    - f{C16H22O}          - Chemical formula
    - 1\@y12              - Multi-analyte spectrum (analyte 1)
    - 0\@b2{LC[Carbamidomethyl]} - Embedded ProForma sequence

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI MzPAFAnnotation
  {
    std::optional<int> analyte_index;               ///< Analyte index (0\@, 1\@, etc.) for multi-analyte spectra
    MzPAFIonSeries ion_series = MzPAFIonSeries::UNKNOWN; ///< Ion series type (a, b, c, x, y, z, etc.)
    std::optional<int> ordinal;                     ///< Position/ordinal (4 in y4)
    std::optional<char> immonium_residue;           ///< Residue for immonium ions (Y in IY)
    std::optional<std::pair<int, int>> internal_range; ///< Range for internal fragments (m3:6 -> {3,6})
    std::optional<String> reporter_name;            ///< Name for reporter ions (r[TMT127N])
    std::optional<EmpiricalFormula> formula;        ///< Chemical formula for formula ions (f{C16H22O})
    std::optional<String> named_compound;           ///< Name for named compound ions (_[name])
    std::vector<MzPAFNeutralLoss> neutral_losses;   ///< Neutral losses (-H2O, -NH3, etc.)
    std::optional<int> isotope_offset;              ///< Isotope offset (+1i, +2i for M+1, M+2)
    std::optional<EmpiricalFormula> adduct;         ///< Adduct ion (+Na, +K, etc.)
    std::optional<int> charge;                      ///< Charge state (^2, ^3)
    std::optional<MzPAFMassDelta> mass_delta;       ///< Mass delta (/0.001, /-1.4ppm)
    std::optional<double> confidence;               ///< Confidence score (*0.75)
    std::optional<String> embedded_sequence;        ///< Embedded ProForma sequence string ({LC[Carbamidomethyl]})

    /// Check if this annotation has minimal valid data
    bool isValid() const;

    /// Equality comparison
    bool operator==(const MzPAFAnnotation& other) const;
    bool operator!=(const MzPAFAnnotation& other) const { return !(*this == other); }
  };

  /// Output operator for MzPAFAnnotation (for testing/debugging)
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const MzPAFAnnotation& ann);

  /**
    @brief Multiple mzPAF annotations for a single peak

    Represents comma-separated alternative annotations for a peak.
    Example: b2,y4^2 has two annotations.

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI MzPAFPeakAnnotations
  {
    std::vector<MzPAFAnnotation> annotations;  ///< List of alternative annotations

    /// Check if there are any annotations
    bool empty() const { return annotations.empty(); }

    /// Get number of annotations
    size_t size() const { return annotations.size(); }

    /// Equality comparison
    bool operator==(const MzPAFPeakAnnotations& other) const;
    bool operator!=(const MzPAFPeakAnnotations& other) const { return !(*this == other); }
  };

  //--------------------------------------------------------------------------
  // Error Handling
  //--------------------------------------------------------------------------

  /**
    @brief Error codes for mzPAF parsing errors

    @ingroup Chemistry
  */
  enum class MzPAFErrorCode
  {
    UNEXPECTED_CHARACTER,       ///< Unexpected character encountered
    UNCLOSED_BRACKET,           ///< Opening bracket without matching close
    INVALID_ION_SERIES,         ///< Invalid ion series character
    INVALID_NUMBER,             ///< Invalid numeric value
    INVALID_FORMULA,            ///< Invalid chemical formula
    INVALID_CHARGE,             ///< Invalid charge specification
    INVALID_DELTA,              ///< Invalid mass delta specification
    INVALID_CONFIDENCE,         ///< Invalid confidence score
    EMPTY_INPUT,                ///< Empty input string
    UNEXPECTED_END_OF_INPUT,    ///< Unexpected end of input
    INTERNAL_ERROR              ///< Internal parser error
  };

  /**
    @brief Convert error code to human-readable string

    @param[in] code The error code to convert
    @return A string describing the error code
  */
  OPENMS_DLLAPI const char* mzPAFErrorCodeToString(MzPAFErrorCode code);

  /**
    @brief Parse error for mzPAF notation

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI MzPAFParseError : public Exception::ParseError
  {
  public:
    MzPAFParseError(
      const char* file,
      int line,
      const char* function,
      MzPAFErrorCode error_code,
      size_t error_position,
      const String& input,
      const String& message
    ) noexcept;

    MzPAFErrorCode getErrorCode() const noexcept { return code_; }
    size_t getPosition() const noexcept { return position_; }
    String getFormattedMessage() const;

  private:
    MzPAFErrorCode code_;
    size_t position_;
    String context_before_;
    String context_after_;

    void extractContext_(const String& input, size_t pos);
  };

  //--------------------------------------------------------------------------
  // Parser API
  //--------------------------------------------------------------------------

  /**
    @brief Parser and writer for mzPAF (Peak Annotation Format) notation

    mzPAF is the HUPO-PSI standard string notation for fragment ion peak annotations,
    used in mzSpecLib and other spectral library formats.

    This class provides static methods for:
    - Parsing mzPAF strings into MzPAFAnnotation structures
    - Converting MzPAFAnnotation structures back to mzPAF strings
    - Integration with OpenMS PeptideHit::PeakAnnotation
    - Validation and format detection

    Usage example:
    @code
    // Parse a single annotation
    MzPAFAnnotation ann = MzPAF::parse("y4^2-H2O/0.001*0.75");

    // Parse multiple annotations (comma-separated)
    MzPAFPeakAnnotations anns = MzPAF::parseMultiple("b2,y4^2");

    // Convert back to string
    String s = MzPAF::toString(ann);

    // Check if a string is mzPAF format
    if (MzPAF::isMzPAFFormat("y4^2")) { ... }
    @endcode

    @see MzPAFAnnotation, MzPAFPeakAnnotations

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI MzPAF
  {
  public:
    //--------------------------------------------------------------------------
    // Parsing
    //--------------------------------------------------------------------------

    /**
      @brief Parse an mzPAF string into a single annotation

      @param[in] input The mzPAF string to parse
      @return The parsed annotation
      @throws MzPAFParseError if parsing fails

      @note If the input contains multiple comma-separated annotations, only the first is returned.
            Use parseMultiple() for multi-annotation strings.
    */
    static MzPAFAnnotation parse(const String& input);

    /**
      @brief Parse an mzPAF string with potentially multiple annotations

      @param[in] input The mzPAF string to parse (may contain comma-separated annotations)
      @return All parsed annotations
      @throws MzPAFParseError if parsing fails
    */
    static MzPAFPeakAnnotations parseMultiple(const String& input);

    /**
      @brief Try to parse an mzPAF string (non-throwing)

      @param[in] input The mzPAF string to parse
      @return The parsed annotation, or std::nullopt on failure
    */
    static std::optional<MzPAFAnnotation> tryParse(const String& input);

    /**
      @brief Try to parse multiple annotations (non-throwing)

      @param[in] input The mzPAF string to parse
      @return The parsed annotations, or empty on failure
    */
    static std::optional<MzPAFPeakAnnotations> tryParseMultiple(const String& input);

    //--------------------------------------------------------------------------
    // Serialization
    //--------------------------------------------------------------------------

    /**
      @brief Convert an annotation to mzPAF string

      @param[in] ann The annotation to convert
      @return The mzPAF string representation
    */
    static String toString(const MzPAFAnnotation& ann);

    /**
      @brief Convert multiple annotations to mzPAF string

      @param[in] anns The annotations to convert
      @return The mzPAF string representation (comma-separated)
    */
    static String toString(const MzPAFPeakAnnotations& anns);

    //--------------------------------------------------------------------------
    // PeakAnnotation Integration
    //--------------------------------------------------------------------------

    /**
      @brief Create a PeptideHit::PeakAnnotation from mzPAF data

      Converts an MzPAFAnnotation to the OpenMS PeakAnnotation format used in PeptideHit.

      @param[in] mzpaf The mzPAF annotation
      @param[in] mz The observed m/z value
      @param[in] intensity The peak intensity
      @return A PeptideHit::PeakAnnotation
    */
    static PeptideHit::PeakAnnotation toPeakAnnotation(
      const MzPAFAnnotation& mzpaf, double mz, double intensity);

    /**
      @brief Parse mzPAF annotations from a PeptideHit::PeakAnnotation

      Attempts to parse the annotation string from a PeakAnnotation as mzPAF.

      @param[in] peak_annotation The PeakAnnotation to parse
      @return Parsed mzPAF annotations, or empty if not valid mzPAF
    */
    static MzPAFPeakAnnotations fromPeakAnnotation(
      const PeptideHit::PeakAnnotation& peak_annotation);

    //--------------------------------------------------------------------------
    // Utilities
    //--------------------------------------------------------------------------

    /**
      @brief Check if a string appears to be in mzPAF format

      Performs a quick heuristic check to determine if a string looks like mzPAF notation.

      @param[in] annotation The string to check
      @return True if the string appears to be mzPAF format
    */
    static bool isMzPAFFormat(const String& annotation);

    /**
      @brief Calculate theoretical m/z for an annotation

      Calculates the theoretical m/z value for the annotated ion based on the sequence.

      @param[in] ann The annotation
      @param[in] sequence The peptide sequence
      @return Theoretical m/z, or std::nullopt if calculation not possible
    */
    static std::optional<double> calculateTheoreticalMZ(
      const MzPAFAnnotation& ann, const AASequence& sequence);

    /**
      @brief Check if ion series is a standard fragment ion (a, b, c, x, y, z)

      @param[in] series The ion series to check
      @return True if it's a standard fragment ion type
    */
    static bool isStandardFragmentIon(MzPAFIonSeries series);

    /**
      @brief Get the ion series character for an annotation

      @param[in] series The ion series enum
      @return The character representation (a, b, c, x, y, z, p, I, m, r, f, _)
    */
    static char ionSeriesToChar(MzPAFIonSeries series);

    /**
      @brief Parse ion series from character

      @param[in] c The character to parse
      @param[out] series Output ion series if successful
      @return True if parsing successful
    */
    static bool charToIonSeries(char c, MzPAFIonSeries& series);

  private:
    MzPAF() = delete;  // Static class, no instantiation
  };

} // namespace OpenMS
