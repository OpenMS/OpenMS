// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ayesha Feroz $
// $Authors: Ayesha Feroz, Tom Müller$
// --------------------------------------------------------------------------
#ifndef OPENMS_CHEMISTRY_PROFORMA_H
#define OPENMS_CHEMISTRY_PROFORMA_H

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

    class OPENMS_DLLAPI ProForma
    {
    public:
        // Parse the ProForma string and populate the hash map
        void fromProFormaString(const String& proforma_str);

        // Convert to ProForma string
        String toProFormaString() const;

        // Add a modification to the sequence at a specific position
        void addModification(size_t start_pos, size_t end_pos, const String& mod_id, double mass_shift);

        // Remove a modification at a specific position
        void removeModification(size_t position);

    private:
        AASequence sequence_;
        std::unordered_map<size_t, ResidueModification> modifications_;
        std::unordered_set<String> supported_cvs_{"UNIMOD", "MOD", "RESID", "XLMOD", "GNO"};

        // Parsing methods
        void parseCVModificationNames(const String& modString, size_t& pos, size_t residue_pos);
        void parseStandardModification(const String& modString, size_t& pos, size_t residue_pos);
        void parseDeltaMassNotation(const String& modString, size_t& pos, size_t residue_pos);
        void parseNTerminalModification(const String& modString, size_t& pos);
        void parseCTerminalModification(const String& modString, size_t& pos);
        // NEW: Parsing method for range modifications
        void parseRangeModification(const String& modString, size_t& pos);
        void throwParseError(const String& message) const;
        void validateCVModification(const String& modification);
    };
}

#endif // OPENMS_CHEMISTRY_PROFORMA_H
