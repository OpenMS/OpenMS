#ifndef OPENMS_CHEMISTRY_PROFORMA_H
#define OPENMS_CHEMISTRY_PROFORMA_H

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
    struct ModificationAttributes
    {
        double mass_shift = 0.0;
        bool ambiguous_start = false;
        bool stop_position = false;
        std::string modification_name;
    };

    class OPENMS_DLLAPI ProForma
    {
    public:
        // Constructor
        explicit ProForma(const AASequence& seq);

        // Parse the ProForma string and populate the hash map
        AASequence fromProFormaString(const std::string& proforma_str);

        // Convert to ProForma string
        std::string toProFormaString() const;

        // Add a modification to the sequence at a specific position
        void addModification(size_t position, const std::string& mod_id, double mass_shift);

        // Remove a modification at a specific position
        void removeModification(size_t position);

    private:
        AASequence sequence_;
        std::unordered_map<size_t, ModificationAttributes> modifications_;
        std::unordered_set<std::string> supported_cvs_{"UNIMOD", "MOD", "RESID", "XLMOD", "GNO"};

        // Parsing methods
        void parseCVModificationNames(const std::string& modString, size_t& pos, size_t residue_pos);
        void parseStandardModification(const std::string& modString, size_t& pos, size_t residue_pos);
        void parseDeltaMassNotation(const std::string& modString, size_t& pos, size_t residue_pos);
        void parseNTerminalModification(const std::string& modString, size_t& pos);
        void parseCTerminalModification(const std::string& modString, size_t& pos);
        void throwParseError(const std::string& message) const;
        void validateCVModification(const std::string& modification);
    };
}

#endif // OPENMS_CHEMISTRY_PROFORMA_H
