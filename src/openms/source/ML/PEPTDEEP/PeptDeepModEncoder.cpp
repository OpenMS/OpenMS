// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PEPTDEEP/PeptDeepModEncoder.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <unordered_map>

namespace OpenMS::ML
{
  namespace
  {
    /// AlphaPeptDeep's mod_elements, in model order. The index of a symbol here IS its
    /// position in the per-site feature vector, so this list must not be reordered.
    const std::vector<std::string>& modElements()
    {
      static const std::vector<std::string> elements = {
        "C", "H", "N", "O", "P", "S", "B", "F",
        "I", "K", "U", "V", "W", "X", "Y", "Ac",
        "Ag", "Al", "Am", "Ar", "As", "At", "Au", "Ba",
        "Be", "Bi", "Bk", "Br", "Ca", "Cd", "Ce", "Cf",
        "Cl", "Cm", "Co", "Cr", "Cs", "Cu", "Dy", "Er",
        "Es", "Eu", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge",
        "He", "Hf", "Hg", "Ho", "In", "Ir", "Kr", "La",
        "Li", "Lr", "Lu", "Md", "Mg", "Mn", "Mo", "Na",
        "Nb", "Nd", "Ne", "Ni", "No", "Np", "Os", "Pa",
        "Pb", "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra",
        "Rb", "Re", "Rh", "Rn", "Ru", "Sb", "Sc", "Se",
        "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te",
        "Th", "Ti", "Tl", "Tm", "Xe", "Yb", "Zn", "Zr",
        "2H", "13C", "15N", "18O", "?"};
      return elements;
    }

    /// Symbol -> index, including the OpenMS isotope spellings. OpenMS writes isotopes as
    /// "(13)C" where PeptDeep writes "13C", so both spellings resolve to the same slot.
    const std::unordered_map<std::string, size_t>& symbolIndex()
    {
      static const std::unordered_map<std::string, size_t> index = []()
      {
        std::unordered_map<std::string, size_t> m;
        const auto& els = modElements();
        for (size_t i = 0; i < els.size(); ++i) { m[els[i]] = i; }
        m["(2)H"] = m.at("2H");
        m["(13)C"] = m.at("13C");
        m["(15)N"] = m.at("15N");
        m["(18)O"] = m.at("18O");
        return m;
      }();
      return index;
    }

    /// Adds the atom counts of one modification's diff formula into a single mod_x row.
    void addDiffFormula_(const ResidueModification& mod, float* row)
    {
      const EmpiricalFormula& diff = mod.getDiffFormula();
      for (auto it = diff.begin(); it != diff.end(); ++it)
      {
        // A negative count is normal: a modification that removes atoms (e.g. a loss)
        // has a negative contribution, exactly as in AlphaPeptDeep.
        row[PeptDeepModEncoder::elementIndex(it->first->getSymbol())] += static_cast<float>(it->second);
      }
    }
  } // namespace

  size_t PeptDeepModEncoder::elementIndex(const std::string& symbol)
  {
    const auto& index = symbolIndex();
    const auto it = index.find(symbol);
    // Anything PeptDeep does not model explicitly is folded into the "?" slot rather
    // than dropped, so its mass delta still registers as "something is modified here".
    return it == index.end() ? unknownElementIndex() : it->second;
  }

  size_t PeptDeepModEncoder::unknownElementIndex()
  {
    return modElements().size() - 1;
  }

  const std::vector<std::string>& PeptDeepModEncoder::elementOrder()
  {
    return modElements();
  }

  void PeptDeepModEncoder::encode(const AASequence& seq, float* mod_x, size_t rows)
  {
    const size_t required = seq.size() + 2;
    if (rows < required)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "mod_x has fewer rows than the peptide needs (length+2).", std::to_string(rows));
    }
    if (!seq.isModified()) { return; }

    // Row 0 is the N-terminus, rows 1..nAA the residues, row nAA+1 the C-terminus.
    if (seq.hasNTerminalModification())
    {
      addDiffFormula_(*seq.getNTerminalModification(), mod_x);
    }
    for (size_t i = 0; i < seq.size(); ++i)
    {
      const ResidueModification* mod = seq[i].getModification();
      if (mod != nullptr)
      {
        addDiffFormula_(*mod, mod_x + (i + 1) * MOD_FEATURE_SIZE);
      }
    }
    if (seq.hasCTerminalModification())
    {
      addDiffFormula_(*seq.getCTerminalModification(), mod_x + (seq.size() + 1) * MOD_FEATURE_SIZE);
    }
  }
} // namespace OpenMS::ML
