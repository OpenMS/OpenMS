// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ayesha Feroz $
// $Authors: Ayesha Feroz, Tom Müller$
// --------------------------------------------------------------------------
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace OpenMS
{
    ProForma::ProForma(const AASequence& seq) : sequence_(seq) {}

    void ProForma::validateCVModification(const std::string& modification)
    {
        size_t colon_pos = modification.find(':');
        if (colon_pos == std::string::npos)
        {
            std::cout << "No CV prefix found in modification: " << modification << std::endl;
            return;
        }

        std::string cv = modification.substr(0, colon_pos);
        std::string accession = modification.substr(colon_pos + 1);

        std::cout << "Validating CV: " << cv << " with accession: " << accession << std::endl;

        if (supported_cvs_.find(cv) == supported_cvs_.end())
        {
            std::cout << "Unsupported CV detected: " << cv << std::endl;
            throw std::invalid_argument("Unsupported CV/ontology: " + cv);
        }

        if (accession.empty())
        {
            throw std::invalid_argument("Accession number cannot be empty in modification: " + modification);
        }
    }

    void ProForma::parseCVModificationNames(const std::string& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == std::string::npos || modEnd == std::string::npos)
        {
            throwParseError("Invalid modification format: Missing brackets for CV modification.");
        }

        std::string modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        std::cout << "Parsing CV modification: " << modification << " at position " << residue_pos << std::endl;

        validateCVModification(modification);

        ModificationAttributes attributes;
        attributes.modification_name = modification;

        if (modString.size() > modEnd + 1 && modString[modEnd + 1] == '?')
        {
            attributes.ambiguous_start = true;
            modEnd++;
        }

        modifications_[residue_pos] = attributes;
        pos = modEnd + 1;
    }

    void ProForma::parseStandardModification(const std::string& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == std::string::npos || modEnd == std::string::npos)
        {
            throwParseError("Invalid modification format: Missing brackets for standard modification.");
        }

        std::string modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        std::cout << "Parsing standard modification: " << modification << " at position " << residue_pos << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;

        modifications_[residue_pos] = attributes;
        pos = modEnd + 1;
    }

    void ProForma::parseDeltaMassNotation(const std::string& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == std::string::npos || modEnd == std::string::npos)
        {
            throwParseError("Invalid mass shift notation: Missing brackets.");
        }

        std::string modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        std::cout << "Parsing mass shift: " << modification << " at position " << residue_pos << std::endl;

        if (modification[0] != '+' && modification[0] != '-')
        {
            throwParseError("Invalid mass shift format: Missing +/- sign.");
        }

        try
        {
            double mass_shift = std::stod(modification);
            ModificationAttributes attributes;
            attributes.mass_shift = mass_shift;

            if (modString.size() > modEnd + 1 && modString[modEnd + 1] == '?')
            {
                attributes.ambiguous_start = true;
                modEnd++;
            }

            modifications_[residue_pos] = attributes;
        }
        catch (const std::invalid_argument&)
        {
            throwParseError("Invalid mass shift format: Could not convert to double.");
        }

        pos = modEnd + 1;
    }

    void ProForma::parseNTerminalModification(const std::string& modString, size_t& pos)
    {
        size_t modEnd = modString.find("]-", pos);
        if (modEnd == std::string::npos)
        {
            throwParseError("Invalid N-terminal modification format: Missing brackets and '-' indicator.");
        }

        std::string modification = modString.substr(1, modEnd - 1);
        std::cout << "Parsing N-terminal modification: " << modification << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;
        modifications_[0] = attributes;

        pos = modEnd + 2;
    }

    void ProForma::parseCTerminalModification(const std::string& modString, size_t& pos)
    {
        size_t modStart = modString.find("-[", pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == std::string::npos || modEnd == std::string::npos)
        {
            throwParseError("Invalid C-terminal modification format: Missing brackets and '-' indicator.");
        }

        std::string modification = modString.substr(modStart + 2, modEnd - modStart - 2);
        std::cout << "Parsing C-terminal modification: " << modification << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;
        modifications_[sequence_.size() + 1] = attributes;

        pos = modEnd + 1;
    }

    void ProForma::removeModification(size_t position)
    {
        std::cout << "Attempting to remove modification at position: " << position << std::endl;

        auto it = modifications_.find(position);
        if (it != modifications_.end())
        {
            std::cout << "Removing modification: " << it->second.modification_name << " at position: " << position << std::endl;
            modifications_.erase(it);
        }
        else
        {
            std::cerr << "No modification found at position: " << position << std::endl;
        }
    }

    void ProForma::addModification(size_t position, const std::string& mod_id, double mass_shift)
    {
        modifications_[position] = {mass_shift, false, false, mod_id};
    }

    AASequence ProForma::fromProFormaString(const std::string& proforma_str)
    {
        AASequence seq;
        size_t pos = 0;
        size_t residue_pos = 0;

        while (pos < proforma_str.size())
        {
            try
            {
                if (std::isalpha(proforma_str[pos]))
                {
                    seq = seq + ResidueDB::getInstance()->getResidue(proforma_str.substr(pos, 1));
                    residue_pos = seq.size();
                    pos++;
                }
                else if (proforma_str[pos] == '[')
                {
                    if (proforma_str.find(':', pos) != std::string::npos)
                    {
                        parseCVModificationNames(proforma_str, pos, residue_pos);
                    }
                    else if (proforma_str[pos + 1] == '+' || proforma_str[pos + 1] == '-')
                    {
                        parseDeltaMassNotation(proforma_str, pos, residue_pos);
                    }
                    else
                    {
                        parseStandardModification(proforma_str, pos, residue_pos);
                    }
                }
                else if (proforma_str[pos] == '-' && proforma_str[pos + 1] == '[')
                {
                    parseCTerminalModification(proforma_str, pos);
                }
                else if (proforma_str[pos] == '[' && proforma_str.find("]-") != std::string::npos)
                {
                    parseNTerminalModification(proforma_str, pos);
                }
                else
                {
                    pos++;
                }
            }
            catch (const std::runtime_error& e)
            {
                throw std::runtime_error("Error parsing ProForma string: " + std::string(e.what()));
            }
        }
        return seq;
    }

    std::string ProForma::toProFormaString() const
    {
        std::stringstream ss;

        if (modifications_.find(0) != modifications_.end())
        {
            ss << "[" << modifications_.at(0).modification_name << "]-";
        }

        for (size_t i = 1; i <= sequence_.size(); ++i)
        {
            ss << sequence_.getResidue(i - 1).getOneLetterCode();

            auto it = modifications_.find(i);
            if (it != modifications_.end())
            {
                if (!it->second.modification_name.empty())
                {
                    ss << "[" << it->second.modification_name << "]";
                }
                if (it->second.mass_shift != 0.0)
                {
                    ss << "[" << (it->second.mass_shift > 0 ? "+" : "") << it->second.mass_shift << "]";
                }
                if (it->second.ambiguous_start)
                {
                    ss << "?";
                }
            }
        }

        if (modifications_.find(sequence_.size() + 1) != modifications_.end())
        {
            ss << "-[" << modifications_.at(sequence_.size() + 1).modification_name << "]";
        }

        return ss.str();
    }

    void ProForma::throwParseError(const std::string& message) const
    {
        throw std::runtime_error("ProForma parsing error: " + message);
    }
}


