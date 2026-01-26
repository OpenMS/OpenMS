// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <map>
#include <vector>
#include <string>

namespace OpenMS
{
  /** @ingroup Chemistry

      @brief Singleton database of monosaccharides for glycan notation

      This class provides access to monosaccharide symbols used in ProForma glycan notation.
      The data is loaded from a JSON file (CHEMISTRY/monosaccharides.json) which contains
      monosaccharide symbols, their full names, monoisotopic masses, chemical formulas,
      and synonyms.

      The monosaccharide data originates from the HUPO-PSI ProForma specification:
      https://github.com/HUPO-PSI/ProForma

      Common monosaccharides include:
      - Hex (Hexose, 162.0528 Da)
      - HexNAc (N-Acetylhexosamine, 203.0794 Da)
      - Fuc (Fucose, 146.0579 Da)
      - NeuAc / Neu5Ac (N-Acetylneuraminic acid / Sialic Acid, 291.0954 Da)
      - NeuGc / Neu5Gc (N-Glycolylneuraminic acid, 307.0903 Da)

      @see https://github.com/HUPO-PSI/ProForma for the ProForma specification
  */
  class OPENMS_DLLAPI MonosaccharideDB
  {
  public:
    /// Structure representing a monosaccharide
    struct OPENMS_DLLAPI Monosaccharide
    {
      String symbol;                    ///< Primary symbol (e.g., "Hex", "HexNAc")
      String name;                      ///< Full name or description
      double mass;                      ///< Monoisotopic mass in Daltons
      String formula;                   ///< Chemical formula (e.g., "C6H10O5")
      std::vector<String> synonyms;     ///< Alternative symbols/names
    };

    /** @name Accessors */
    //@{

    /// Returns a pointer to the singleton instance of the monosaccharide database
    /// This is thread-safe upon first and subsequent calls (Meyers' singleton)
    static MonosaccharideDB* getInstance();

    /// Check if a symbol (or synonym) is a known monosaccharide
    bool hasSymbol(const String& symbol) const;

    /**
     * @brief Get monosaccharide by symbol
     * @param symbol The monosaccharide symbol (e.g., "Hex") or a known synonym
     * @return Pointer to the monosaccharide, or nullptr if not found
     */
    const Monosaccharide* getMonosaccharide(const String& symbol) const;

    /**
     * @brief Get monosaccharide by symbol (throws if not found)
     * @param symbol The monosaccharide symbol or synonym
     * @return Reference to the monosaccharide
     * @throw Exception::ElementNotFound if symbol is not known
     */
    const Monosaccharide& getMonosaccharideOrThrow(const String& symbol) const;

    /// Get all known primary symbols
    std::vector<String> getAllSymbols() const;

    /// Get the number of monosaccharides in the database
    Size getNumberOfMonosaccharides() const;

    //@}

  private:
    /// Private constructor (singleton pattern)
    MonosaccharideDB();

    /// Destructor
    ~MonosaccharideDB() = default;

    /// Deleted copy constructor
    MonosaccharideDB(const MonosaccharideDB&) = delete;

    /// Deleted move constructor
    MonosaccharideDB(MonosaccharideDB&&) = delete;

    /// Deleted copy assignment
    MonosaccharideDB& operator=(const MonosaccharideDB&) = delete;

    /// Deleted move assignment
    MonosaccharideDB& operator=(MonosaccharideDB&&) = delete;

    /// Load monosaccharide data from JSON file
    void loadFromJSON_();

    /// Map from primary symbol to monosaccharide data
    std::map<String, Monosaccharide> monosaccharides_;

    /// Map from synonyms to primary symbol (for lookup)
    std::map<String, String> synonym_to_symbol_;
  };

} // namespace OpenMS
