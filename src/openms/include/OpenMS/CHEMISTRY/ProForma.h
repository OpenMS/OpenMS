// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ayesha Feroz $
// $Authors: Ayesha Feroz, Tom Müller$
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <vector>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
    // Forward declarations
    class ResidueModification;

    /**
     * @brief Parser and representation for ProForma peptidoform notation
     *
     * This class implements parsing and generation of ProForma strings according to
     * the HUPO-PSI ProForma specification version 2.0. ProForma is a standardized
     * notation for representing modified peptide sequences with precise modification
     * information.
     *
     * @section proforma_coverage ProForma 2.0 Coverage
     *
     * This implementation supports the following ProForma 2.0 features:
     *
     * @subsection supported_features Supported Features
     *
     * - **Basic amino acid sequences**: All 20 standard amino acids plus unconventional (U/O) and ambiguous (X/J/B/Z) amino acids
     * - **Controlled vocabulary modifications**: Support for UNIMOD, PSI-MOD, RESID, XLMOD, and GNOme ontologies
     *   - Format: `[CV:accession]` (e.g., `[UNIMOD:35]`, `[MOD:00046]`, `[RESID:AA0037]`)
     * - **Named modifications**: User-defined modification names in brackets
     *   - Format: `[ModificationName]` (e.g., `[Phospho]`, `[Acetyl]`)
     * - **Delta mass modifications**: Mass shift notation with explicit sign
     *   - Format: `[+/-mass]` (e.g., `[+15.99]`, `[-17.026]`)
     * - **N-terminal modifications**: Modifications at the peptide N-terminus
     *   - Format: `[modification]-sequence` (e.g., `[Acetyl]-PEPTIDE`)
     * - **C-terminal modifications**: Modifications at the peptide C-terminus
     *   - Format: `sequence-[modification]` (e.g., `PEPTIDE-[Amidation]`)
     * - **Range modifications**: Modifications spanning multiple consecutive residues
     *   - Format: `prefix(range)[modification]suffix` (e.g., `A(CDE)[+12.5]FGH`)
     * - **Ambiguous modifications**: Basic support for ambiguous modification placement
     *   - Format: `[modification]?` (e.g., `[+15.99]?`)
     *
     * @subsection unsupported_features Unsupported ProForma 2.0 Features
     *
     * The following ProForma 2.0 features are **NOT** currently implemented:
     *
     * - **Formula-based modifications**: `Formula:` notation (e.g., `[Formula:C2H6:z+2]`)
     * - **Glycan modifications**: Glycan composition notation with monosaccharides
     * - **Cross-linking modifications**: `#XL` and `#BRANCH` notation for cross-linked peptides
     * - **Labile modifications**: Modifications that can dissociate during analysis
     * - **Fixed modifications**: Global modifications applied to all matching residues
     * - **Isotope replacements**: Global isotope labeling (e.g., heavy nitrogen)
     * - **Charge states**: Global charge state information
     * - **Compound peptidoform ions**: Multiple peptidoforms in a single notation
     * - **Names/labels**: `(>name)`, `(>>name)`, `(>>>name)` notation
     * - **Advanced ambiguous modifications**:
     *   - Position rules and scoring
     *   - `#label` notation for unlocalised modifications
     *   - `COMKP`/`COMUP` flags for co-localization rules
     *   - Occurrence limits for ambiguous modifications
     * - **INFO tags**: `INFO:` metadata annotations
     *
     * @section examples Usage Examples
     *
     * @code{.cpp}
     * ProForma proforma;
     *
     * // Parse a simple modified peptide
     * proforma.fromProFormaString("EM[UNIMOD:35]EVEES[Phospho]PEK");
     *
     * // Convert back to ProForma string
     * String proforma_str = proforma.toProFormaString();
     *
     * // Add modifications programmatically
     * proforma.fromProFormaString("PEPTIDE");
     * proforma.addModification(0, 0, "Acetyl", 42.010565);  // N-terminal
     * proforma.addModification(3, 3, "", 79.966331);        // Mass-only at position 3
     * proforma.addModification(8, 8, "Amidation", 0.0);     // C-terminal
     *
     * // Remove modifications
     * proforma.removeModification(3);
     * @endcode
     *
     * @section specification_reference ProForma Specification Reference
     *
     * This implementation is based on the HUPO-PSI ProForma specification version 2.0.
     * The complete specification schema is available in the proforma.json file.
     *
     * For more information, see:
     * - ProForma specification: https://www.psidev.info/proforma
     * - HUPO-PSI GitHub: https://github.com/HUPO-PSI/ProForma
     *
     * @note This is a partial implementation focusing on the most commonly used
     *       ProForma features in mass spectrometry workflows. Additional features
     *       may be added in future versions.
     */
    class OPENMS_DLLAPI ProForma
    {
    public:
        /**
         * @brief Parse a ProForma string and populate the internal representation
         *
         * Parses a ProForma notation string and extracts the amino acid sequence
         * and modification information. The method supports all implemented ProForma
         * features (see class documentation for coverage details).
         *
         * @param proforma_str The ProForma notation string to parse
         *
         * @throws Exception::ParseError If the ProForma string has invalid syntax
         * @throws Exception::IllegalArgument If unsupported CV/ontology is used
         *
         * @section parsing_examples Parsing Examples
         *
         * @code{.cpp}
         * ProForma pf;
         *
         * // Basic sequence with named modification
         * pf.fromProFormaString("PEPTIDE[Phospho]");
         *
         * // CV-based modifications
         * pf.fromProFormaString("EM[UNIMOD:35]EVEES[MOD:00046]PEK");
         *
         * // Mass shift modifications
         * pf.fromProFormaString("PEPT[+79.966]IDE");
         *
         * // Terminal modifications
         * pf.fromProFormaString("[Acetyl]-PEPTIDE-[Amidation]");
         *
         * // Range modifications
         * pf.fromProFormaString("PEP(TIDE)[+12.5]");
         *
         * // Ambiguous modifications
         * pf.fromProFormaString("PEPTIDE[+15.99]?");
         * @endcode
         */
        void fromProFormaString(const String& proforma_str);

        /**
         * @brief Convert the internal representation to a ProForma string
         *
         * Generates a ProForma notation string from the current amino acid sequence
         * and modifications. The output follows ProForma 2.0 syntax conventions.
         *
         * @return ProForma notation string representing the current peptidoform
         *
         * @note The output may not be identical to the input string due to
         *       normalization of modification representations, but it will be
         *       semantically equivalent.
         */
        String toProFormaString() const;

        /**
         * @brief Add a modification to the sequence at specific position(s)
         *
         * Programmatically adds a modification to the peptide sequence. Supports
         * both single-position and range modifications.
         *
         * @param start_pos Starting position for the modification (0 = N-terminus,
         *                  1 = first residue, sequence_length+1 = C-terminus)
         * @param end_pos Ending position for the modification (same as start_pos for
         *                single-position modifications)
         * @param mod_id Modification identifier or name (empty string for mass-only modifications)
         * @param mass_shift Mass shift in Daltons (0.0 if not applicable)
         *
         * @throws Exception::OutOfRange If position is beyond sequence boundaries
         *
         * @section modification_examples Modification Examples
         *
         * @code{.cpp}
         * ProForma pf;
         * pf.fromProFormaString("PEPTIDE");
         *
         * // Add N-terminal acetylation
         * pf.addModification(0, 0, "Acetyl", 42.010565);
         *
         * // Add phosphorylation at position 3
         * pf.addModification(3, 3, "Phospho", 79.966331);
         *
         * // Add mass-only modification
         * pf.addModification(5, 5, "", 15.994915);
         *
         * // Add range modification from position 2 to 4
         * pf.addModification(2, 4, "RangeModification", 12.5);
         *
         * // Add C-terminal amidation
         * pf.addModification(8, 8, "Amidation", -0.984016);  // 8 = length + 1
         * @endcode
         */
        void addModification(size_t start_pos, size_t end_pos, const String& mod_id, double mass_shift);

        /**
         * @brief Remove a modification at a specific position
         *
         * Removes an existing modification from the specified position in the sequence.
         *
         * @param position Position of the modification to remove (0 = N-terminus,
         *                 1 = first residue, etc.)
         *
         * @throws Exception::OutOfRange If position is beyond sequence boundaries
         * @throws Exception::InvalidValue If no modification exists at the specified position
         *
         * @code{.cpp}
         * ProForma pf;
         * pf.fromProFormaString("PEP[Phospho]TIDE");
         * pf.removeModification(3);  // Remove phosphorylation from position 3
         * @endcode
         */
        void removeModification(size_t position);

    private:
        /// The amino acid sequence of the peptide
        AASequence sequence_;
        
        /// Map storing modifications at specific positions (position -> modification)
        std::unordered_map<size_t, ResidueModification> modifications_;
        
        /// Set of supported controlled vocabularies/ontologies
        std::unordered_set<String> supported_cvs_{"UNIMOD", "MOD", "RESID", "XLMOD", "GNO"};

        /**
         * @brief Parse controlled vocabulary (CV) modifications
         *
         * Handles CV-based modifications like [UNIMOD:35], [MOD:00046], etc.
         * Validates the CV prefix and accession format.
         *
         * @param modString The complete ProForma string being parsed
         * @param pos Current parsing position (updated after parsing)
         * @param residue_pos Position in sequence where modification applies
         */
        void parseCVModificationNames_(const String& modString, size_t& pos, size_t residue_pos);
        
        /**
         * @brief Parse named modifications
         *
         * Handles user-defined modification names in square brackets.
         *
         * @param modString The complete ProForma string being parsed
         * @param pos Current parsing position (updated after parsing)
         * @param residue_pos Position in sequence where modification applies
         */
        void parseStandardModification_(const String& modString, size_t& pos, size_t residue_pos);

        /**
         * @brief Parse delta mass notation
         *
         * Handles mass shift modifications like [+15.99] or [-17.026].
         * Supports ambiguous modification markers (?).
         *
         * @param modString The complete ProForma string being parsed
         * @param pos Current parsing position (updated after parsing)
         * @param residue_pos Position in sequence where modification applies
         */
        void parseDeltaMassNotation_(const String& modString, size_t& pos, size_t residue_pos);
        
        /**
         * @brief Parse N-terminal modifications
         *
         * Handles modifications at the N-terminus using [modification]- syntax.
         *
         * @param modString The complete ProForma string being parsed
         * @param pos Current parsing position (updated after parsing)
         */
        void parseNTerminalModification_(const String& modString, size_t& pos);
        
        /**
         * @brief Parse C-terminal modifications
         *
         * Handles modifications at the C-terminus using -[modification] syntax.
         *
         * @param modString The complete ProForma string being parsed
         * @param pos Current parsing position (updated after parsing)
         */
        void parseCTerminalModification_(const String& modString, size_t& pos);
        
        /**
         * @brief Parse range modifications
         *
         * Handles modifications spanning multiple residues using (range)[modification] syntax.
         *
         * @param modString The complete ProForma string being parsed
         * @param pos Current parsing position (updated after parsing)
         */
        void parseRangeModification_(const String& modString, size_t& pos);
        
        /**
         * @brief Validate controlled vocabulary modification format
         *
         * Checks that CV modifications have valid syntax (CV:accession) and
         * use supported controlled vocabularies.
         *
         * @param modification The modification string to validate (without brackets)
         * @throws Exception::IllegalArgument If CV is unsupported or format is invalid
         */
        void validateCVModification_(const String& modification);
    };
}