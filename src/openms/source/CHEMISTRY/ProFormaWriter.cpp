// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProFormaWriter.h>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <type_traits>

namespace OpenMS
{

  String ProFormaWriter::toString(const Peptidoform& peptidoform)
  {
    std::ostringstream os;

    // Write unlocalised modifications: [Phospho]? or [Phospho]^2?
    writeUnlocalisedMods_(os, peptidoform.unlocalised_mods);

    // Write labile modifications: {Glycan:Hex}
    writeLabileModifications_(os, peptidoform.labile_mods);

    // Write N-terminal modifications: [Acetyl]-
    writeNTermMods_(os, peptidoform.n_term_mods);

    // Write the sequence with modifications
    writeSequence_(os, peptidoform.sequence);

    // Write C-terminal modifications: -[Amidated]
    writeCTermMods_(os, peptidoform.c_term_mods);

    return String(os.str());
  }

  String ProFormaWriter::toString(const PeptidoformIon& ion)
  {
    std::ostringstream os;

    // Write global modifications if present (from the first chain that might have them)
    // Note: In ProForma, global mods appear before all chains
    // The AST doesn't have a separate global_mods field in PeptidoformIon,
    // so global mods would need to be handled differently if needed.

    // Write chains separated by //
    bool first = true;
    for (const auto& chain : ion.chains)
    {
      if (!first)
      {
        os << "//";
      }
      first = false;
      os << toString(chain);
    }

    // Write charge state if present
    if (ion.charge.has_value())
    {
      writeChargeState_(os, ion.charge.value());
    }

    return String(os.str());
  }

  void ProFormaWriter::writeGlobalMods_(std::ostream& os, const std::vector<GlobalModEntry>& mods)
  {
    for (const auto& entry : mods)
    {
      std::visit([&os](auto&& arg)
      {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, IsotopeReplacement>)
        {
          ProFormaWriter::writeIsotopeReplacement_(os, arg);
        }
        else if constexpr (std::is_same_v<T, GlobalModification>)
        {
          ProFormaWriter::writeGlobalModification_(os, arg);
        }
      }, entry);
    }
  }

  void ProFormaWriter::writeIsotopeReplacement_(std::ostream& os, const IsotopeReplacement& isotope)
  {
    os << '<' << isotope.isotope << '>';
  }

  void ProFormaWriter::writeGlobalModification_(std::ostream& os, const GlobalModification& mod)
  {
    os << '<';
    writeModification_(os, mod.modification);
    if (!mod.locations.empty())
    {
      os << '@';
      bool first = true;
      for (const auto& loc : mod.locations)
      {
        if (!first)
        {
          os << ',';
        }
        first = false;
        os << loc;
      }
    }
    os << '>';
  }

  void ProFormaWriter::writeModification_(std::ostream& os, const Modification& mod)
  {
    os << '[';
    bool first = true;
    for (const auto& alt : mod.alternatives)
    {
      if (!first)
      {
        os << '|';
      }
      first = false;

      // Write the modification tag
      writeModificationTag_(os, alt.first);

      // Write the label if present
      if (alt.second.has_value())
      {
        writeLabel_(os, alt.second.value());
      }
    }
    os << ']';
  }

  void ProFormaWriter::writeModificationTag_(std::ostream& os, const ModificationTag& tag)
  {
    std::visit([&os](auto&& arg)
    {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, CvAccession>)
      {
        ProFormaWriter::writeCvAccession_(os, arg);
      }
      else if constexpr (std::is_same_v<T, NamedMod>)
      {
        ProFormaWriter::writeNamedMod_(os, arg);
      }
      else if constexpr (std::is_same_v<T, MassDelta>)
      {
        ProFormaWriter::writeMassDelta_(os, arg);
      }
      else if constexpr (std::is_same_v<T, FormulaTag>)
      {
        ProFormaWriter::writeFormulaTag_(os, arg);
      }
      else if constexpr (std::is_same_v<T, GlycanComposition>)
      {
        ProFormaWriter::writeGlycanComposition_(os, arg);
      }
      else if constexpr (std::is_same_v<T, InfoTag>)
      {
        ProFormaWriter::writeInfoTag_(os, arg);
      }
    }, tag);
  }

  void ProFormaWriter::writeCvAccession_(std::ostream& os, const CvAccession& cv)
  {
    os << cvDatabaseToString_(cv.database) << ':' << cv.accession;
  }

  void ProFormaWriter::writeNamedMod_(std::ostream& os, const NamedMod& named)
  {
    if (named.cv_hint.has_value())
    {
      os << cvDatabaseToHint_(named.cv_hint.value()) << ':';
    }
    os << named.name;
  }

  void ProFormaWriter::writeMassDelta_(std::ostream& os, const MassDelta& delta)
  {
    // Write source prefix if present
    if (delta.source != MassDelta::Source::NONE)
    {
      os << massSourceToString_(delta.source) << ':';
    }

    // Use original_text if available for lossless roundtrip
    if (!delta.original_text.empty())
    {
      os << delta.original_text;
    }
    else
    {
      // Format the mass value with sign
      if (delta.mass >= 0)
      {
        os << '+';
      }
      // Use enough precision to preserve the value
      os << std::setprecision(6) << delta.mass;
    }
  }

  void ProFormaWriter::writeFormulaTag_(std::ostream& os, const FormulaTag& formula)
  {
    os << "Formula:" << formula.formula_string;
    if (formula.charge.has_value())
    {
      os << ":z";
      int charge = formula.charge.value();
      if (charge >= 0)
      {
        os << '+';
      }
      os << charge;
    }
  }

  void ProFormaWriter::writeGlycanComposition_(std::ostream& os, const GlycanComposition& glycan)
  {
    os << "Glycan:";
    for (const auto& component : glycan.components)
    {
      std::visit([&os](auto&& mono)
      {
        using T = std::decay_t<decltype(mono)>;
        if constexpr (std::is_same_v<T, String>)
        {
          os << mono;
        }
        else if constexpr (std::is_same_v<T, FormulaTag>)
        {
          os << "Formula:" << mono.formula_string;
          if (mono.charge.has_value())
          {
            os << ":z";
            int charge = mono.charge.value();
            if (charge >= 0)
            {
              os << '+';
            }
            os << charge;
          }
        }
      }, component.first);
      // Only write count if > 1
      if (component.second != 1)
      {
        os << component.second;
      }
    }
  }

  void ProFormaWriter::writeInfoTag_(std::ostream& os, const InfoTag& info)
  {
    os << "INFO:" << info.text;
  }

  void ProFormaWriter::writeLabel_(std::ostream& os, const Label& label)
  {
    os << '#' << label.identifier;
    if (label.score.has_value())
    {
      os << '(' << label.score.value() << ')';
    }
  }

  void ProFormaWriter::writeUnlocalisedMods_(std::ostream& os, const std::vector<UnlocalisedMod>& mods)
  {
    for (const auto& unloc : mods)
    {
      for (const auto& mod : unloc.modifications)
      {
        writeModification_(os, mod);
      }
      if (unloc.occurrence.has_value())
      {
        os << '^' << unloc.occurrence.value();
      }
      os << '?';
    }
  }

  void ProFormaWriter::writeLabileModifications_(std::ostream& os, const std::vector<LabileModification>& mods)
  {
    for (const auto& labile : mods)
    {
      os << '{';
      // Write the modification content without outer brackets
      bool first = true;
      for (const auto& alt : labile.modification.alternatives)
      {
        if (!first)
        {
          os << '|';
        }
        first = false;
        writeModificationTag_(os, alt.first);
        if (alt.second.has_value())
        {
          writeLabel_(os, alt.second.value());
        }
      }
      os << '}';
    }
  }

  void ProFormaWriter::writeNTermMods_(std::ostream& os, const std::vector<Modification>& mods)
  {
    if (!mods.empty())
    {
      for (const auto& mod : mods)
      {
        writeModification_(os, mod);
      }
      os << '-';
    }
  }

  void ProFormaWriter::writeCTermMods_(std::ostream& os, const std::vector<Modification>& mods)
  {
    if (!mods.empty())
    {
      os << '-';
      for (const auto& mod : mods)
      {
        writeModification_(os, mod);
      }
    }
  }

  void ProFormaWriter::writeSequence_(std::ostream& os, const std::vector<SequenceSection>& seq)
  {
    for (const auto& section : seq)
    {
      std::visit([&os](auto&& arg)
      {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, SequenceElement>)
        {
          ProFormaWriter::writeSequenceElement_(os, arg);
        }
        else if constexpr (std::is_same_v<T, AmbiguousRegion>)
        {
          ProFormaWriter::writeAmbiguousRegion_(os, arg);
        }
        else if constexpr (std::is_same_v<T, ModifiedRange>)
        {
          ProFormaWriter::writeModifiedRange_(os, arg);
        }
      }, section);
    }
  }

  void ProFormaWriter::writeSequenceElement_(std::ostream& os, const SequenceElement& elem)
  {
    os << elem.amino_acid;
    for (const auto& mod : elem.modifications)
    {
      writeModification_(os, mod);
    }
  }

  void ProFormaWriter::writeAmbiguousRegion_(std::ostream& os, const AmbiguousRegion& region)
  {
    os << "(?";
    for (const auto& elem : region.elements)
    {
      writeSequenceElement_(os, elem);
    }
    os << ')';
  }

  void ProFormaWriter::writeModifiedRange_(std::ostream& os, const ModifiedRange& range)
  {
    os << '(';
    for (const auto& elem : range.elements)
    {
      writeSequenceElement_(os, elem);
    }
    os << ')';
    for (const auto& mod : range.modifications)
    {
      writeModification_(os, mod);
    }
  }

  void ProFormaWriter::writeChargeState_(std::ostream& os, const ChargeState& charge)
  {
    os << '/';
    std::visit([&os](auto&& arg)
    {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, int>)
      {
        // Simple integer charge
        if (arg >= 0)
        {
          os << arg;
        }
        else
        {
          os << arg;
        }
      }
      else if constexpr (std::is_same_v<T, std::vector<AdductIon>>)
      {
        // Adduct list format: [M+2H]2+ or [Na:z+1,H:z+1]
        os << '[';

        // Calculate total charge for display
        int total_charge = 0;
        bool first = true;
        for (const auto& adduct : arg)
        {
          if (!first)
          {
            os << ',';
          }
          first = false;

          os << adduct.formula << ":z";
          if (adduct.charge >= 0)
          {
            os << '+';
          }
          os << adduct.charge;

          if (adduct.occurrence.has_value())
          {
            os << '^' << adduct.occurrence.value();
          }

          // Add contribution to total charge (considering occurrence)
          int count = adduct.occurrence.value_or(1);
          total_charge += adduct.charge * count;
        }

        os << ']';

        // Write the total charge with sign
        os << std::abs(total_charge);
        if (total_charge >= 0)
        {
          os << '+';
        }
        else
        {
          os << '-';
        }
      }
    }, charge);
  }

  const char* ProFormaWriter::cvDatabaseToString_(CvDatabase db)
  {
    switch (db)
    {
      case CvDatabase::UNIMOD: return "UNIMOD";
      case CvDatabase::MOD:    return "MOD";
      case CvDatabase::RESID:  return "RESID";
      case CvDatabase::XLMOD:  return "XLMOD";
      case CvDatabase::GNO:    return "GNO";
      default:                 return "UNIMOD";
    }
  }

  const char* ProFormaWriter::cvDatabaseToHint_(CvDatabase db)
  {
    switch (db)
    {
      case CvDatabase::UNIMOD: return "U";
      case CvDatabase::MOD:    return "M";
      case CvDatabase::RESID:  return "R";
      case CvDatabase::XLMOD:  return "X";
      case CvDatabase::GNO:    return "G";
      default:                 return "U";
    }
  }

  const char* ProFormaWriter::massSourceToString_(MassDelta::Source source)
  {
    switch (source)
    {
      case MassDelta::Source::OBS: return "Obs";
      case MassDelta::Source::U:   return "U";
      case MassDelta::Source::M:   return "M";
      case MassDelta::Source::R:   return "R";
      case MassDelta::Source::X:   return "X";
      case MassDelta::Source::G:   return "G";
      case MassDelta::Source::NONE:
      default: return "";
    }
  }

} // namespace OpenMS
