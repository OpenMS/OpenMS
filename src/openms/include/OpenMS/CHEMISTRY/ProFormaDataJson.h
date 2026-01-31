// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ProForma.h>
#include <nlohmann/json.hpp>

namespace OpenMS
{

  //--------------------------------------------------------------------------
  // JSON Serialization (nlohmann/json ADL pattern)
  //--------------------------------------------------------------------------

  /// @name CvDatabase enum JSON serialization
  /// @{

  /// Convert CvDatabase enum to JSON string
  inline void to_json(nlohmann::json& j, const ProForma::CvDatabase& db)
  {
    switch (db)
    {
      case ProForma::CvDatabase::UNIMOD: j = "UNIMOD"; break;
      case ProForma::CvDatabase::MOD:    j = "MOD"; break;
      case ProForma::CvDatabase::RESID:  j = "RESID"; break;
      case ProForma::CvDatabase::XLMOD:  j = "XLMOD"; break;
      case ProForma::CvDatabase::GNO:    j = "GNO"; break;
    }
  }

  /// Construct CvDatabase enum from JSON string
  inline void from_json(const nlohmann::json& j, ProForma::CvDatabase& db)
  {
    std::string s = j.get<std::string>();
    if (s == "UNIMOD") db = ProForma::CvDatabase::UNIMOD;
    else if (s == "MOD") db = ProForma::CvDatabase::MOD;
    else if (s == "RESID") db = ProForma::CvDatabase::RESID;
    else if (s == "XLMOD") db = ProForma::CvDatabase::XLMOD;
    else if (s == "GNO") db = ProForma::CvDatabase::GNO;
    else throw std::invalid_argument("Unknown ProForma::CvDatabase: " + s);
  }

  /// @}


  /// @name CvAccession JSON serialization
  /// @{

  /// Convert CvAccession to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::CvAccession& cv)
  {
    j = nlohmann::json{{"database", cv.database}, {"accession", static_cast<std::string>(cv.accession)}};
  }

  /// Construct CvAccession from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::CvAccession& cv)
  {
    j.at("database").get_to(cv.database);
    cv.accession = j.at("accession").get<std::string>();
  }

  /// @}


  /// @name NamedMod JSON serialization
  /// @{

  /// Convert NamedMod to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::NamedMod& nm)
  {
    j = nlohmann::json{{"name", static_cast<std::string>(nm.name)}};
    if (nm.cv_hint.has_value())
    {
      j["cv_hint"] = nm.cv_hint.value();
    }
  }

  /// Construct NamedMod from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::NamedMod& nm)
  {
    nm.name = j.at("name").get<std::string>();
    if (j.contains("cv_hint") && !j.at("cv_hint").is_null())
    {
      nm.cv_hint = j.at("cv_hint").get<ProForma::CvDatabase>();
    }
    else
    {
      nm.cv_hint = std::nullopt;
    }
  }

  /// @}


  /// @name ProForma::MassDelta::Source enum JSON serialization
  /// @{

  /// Convert ProForma::MassDelta::Source enum to JSON string
  inline void to_json(nlohmann::json& j, const ProForma::MassDelta::Source& src)
  {
    switch (src)
    {
      case ProForma::MassDelta::Source::NONE: j = "NONE"; break;
      case ProForma::MassDelta::Source::OBS:  j = "OBS"; break;
      case ProForma::MassDelta::Source::U:    j = "U"; break;
      case ProForma::MassDelta::Source::M:    j = "M"; break;
      case ProForma::MassDelta::Source::R:    j = "R"; break;
      case ProForma::MassDelta::Source::X:    j = "X"; break;
      case ProForma::MassDelta::Source::G:    j = "G"; break;
    }
  }

  /// Construct ProForma::MassDelta::Source enum from JSON string
  inline void from_json(const nlohmann::json& j, ProForma::MassDelta::Source& src)
  {
    std::string s = j.get<std::string>();
    if (s == "NONE") src = ProForma::MassDelta::Source::NONE;
    else if (s == "OBS") src = ProForma::MassDelta::Source::OBS;
    else if (s == "U") src = ProForma::MassDelta::Source::U;
    else if (s == "M") src = ProForma::MassDelta::Source::M;
    else if (s == "R") src = ProForma::MassDelta::Source::R;
    else if (s == "X") src = ProForma::MassDelta::Source::X;
    else if (s == "G") src = ProForma::MassDelta::Source::G;
    else throw std::invalid_argument("Unknown ProForma::MassDelta::Source: " + s);
  }

  /// @}


  /// @name MassDelta JSON serialization
  /// @{

  /// Convert MassDelta to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::MassDelta& md)
  {
    j = nlohmann::json{
      {"source", md.source},
      {"mass", md.mass},
      {"original_text", static_cast<std::string>(md.original_text)}
    };
  }

  /// Construct MassDelta from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::MassDelta& md)
  {
    j.at("source").get_to(md.source);
    j.at("mass").get_to(md.mass);
    md.original_text = j.at("original_text").get<std::string>();
  }

  /// @}


  /// @name FormulaTag JSON serialization
  /// @{

  /// Convert FormulaTag to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::FormulaTag& ft)
  {
    j = nlohmann::json{{"formula_string", static_cast<std::string>(ft.formula_string)}};
    if (ft.charge.has_value())
    {
      j["charge"] = ft.charge.value();
    }
  }

  /// Construct FormulaTag from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::FormulaTag& ft)
  {
    ft.formula_string = j.at("formula_string").get<std::string>();
    if (j.contains("charge") && !j.at("charge").is_null())
    {
      ft.charge = j.at("charge").get<int>();
    }
    else
    {
      ft.charge = std::nullopt;
    }
  }

  /// @}


  /// @name ProForma::GlycanComposition::Monosaccharide JSON serialization
  /// @{

  /// Convert Monosaccharide variant to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::GlycanComposition::Monosaccharide& mono)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, OpenMS::String>)
      {
        j = nlohmann::json{{"type", "name"}, {"value", static_cast<std::string>(arg)}};
      }
      else if constexpr (std::is_same_v<T, ProForma::FormulaTag>)
      {
        j = nlohmann::json{{"type", "formula"}, {"value", arg}};
      }
    }, mono);
  }

  /// Construct Monosaccharide variant from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::GlycanComposition::Monosaccharide& mono)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "name")
    {
      mono = String(j.at("value").get<std::string>());
    }
    else if (type == "formula")
    {
      mono = j.at("value").get<ProForma::FormulaTag>();
    }
    else
    {
      throw std::invalid_argument("Unknown ProForma::GlycanComposition::Monosaccharide type: " + type);
    }
  }

  /// @}


  /// @name GlycanComposition JSON serialization
  /// @{

  /// Convert GlycanComposition to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::GlycanComposition& gc)
  {
    j = nlohmann::json::array();
    for (const auto& [mono, count] : gc.components)
    {
      nlohmann::json component;
      component["monosaccharide"] = mono;
      component["count"] = count;
      j.push_back(component);
    }
  }

  /// Construct GlycanComposition from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::GlycanComposition& gc)
  {
    gc.components.clear();
    for (const auto& item : j)
    {
      ProForma::GlycanComposition::Monosaccharide mono;
      from_json(item.at("monosaccharide"), mono);
      int count = item.at("count").get<int>();
      gc.components.emplace_back(mono, count);
    }
  }

  /// @}


  /// @name InfoTag JSON serialization
  /// @{

  /// Convert InfoTag to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::InfoTag& it)
  {
    j = nlohmann::json{{"text", static_cast<std::string>(it.text)}};
  }

  /// Construct InfoTag from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::InfoTag& it)
  {
    it.text = j.at("text").get<std::string>();
  }

  /// @}


  /// @name PositionConstraint JSON serialization
  /// @{

  /// Convert PositionConstraint to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::PositionConstraint& pc)
  {
    j = nlohmann::json{
      {"residues", std::string(pc.residues.begin(), pc.residues.end())},
      {"n_term", pc.n_term},
      {"c_term", pc.c_term}
    };
  }

  /// Construct PositionConstraint from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::PositionConstraint& pc)
  {
    std::string residue_str = j.at("residues").get<std::string>();
    pc.residues.assign(residue_str.begin(), residue_str.end());
    if (j.contains("n_term")) pc.n_term = j.at("n_term").get<bool>();
    if (j.contains("c_term")) pc.c_term = j.at("c_term").get<bool>();
  }

  /// @}


  /// @name ModificationTag variant JSON serialization
  /// @{

  /// Convert ModificationTag variant to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::ModificationTag& tag)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, ProForma::CvAccession>)
      {
        j = nlohmann::json{{"type", "cv_accession"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::NamedMod>)
      {
        j = nlohmann::json{{"type", "named_mod"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::MassDelta>)
      {
        j = nlohmann::json{{"type", "mass_delta"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::FormulaTag>)
      {
        j = nlohmann::json{{"type", "formula"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::GlycanComposition>)
      {
        j = nlohmann::json{{"type", "glycan"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::InfoTag>)
      {
        j = nlohmann::json{{"type", "info"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::PositionConstraint>)
      {
        j = nlohmann::json{{"type", "position"}, {"value", arg}};
      }
    }, tag);
  }

  /// Construct ModificationTag variant from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::ModificationTag& tag)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "cv_accession")
    {
      tag = j.at("value").get<ProForma::CvAccession>();
    }
    else if (type == "named_mod")
    {
      tag = j.at("value").get<ProForma::NamedMod>();
    }
    else if (type == "mass_delta")
    {
      tag = j.at("value").get<ProForma::MassDelta>();
    }
    else if (type == "formula")
    {
      tag = j.at("value").get<ProForma::FormulaTag>();
    }
    else if (type == "glycan")
    {
      tag = j.at("value").get<ProForma::GlycanComposition>();
    }
    else if (type == "info")
    {
      tag = j.at("value").get<ProForma::InfoTag>();
    }
    else if (type == "position")
    {
      tag = j.at("value").get<ProForma::PositionConstraint>();
    }
    else
    {
      throw std::invalid_argument("Unknown ProForma::ModificationTag type: " + type);
    }
  }

  /// @}


  /// @name ProForma::Label::Type enum JSON serialization
  /// @{

  /// Convert ProForma::Label::Type enum to JSON string
  inline void to_json(nlohmann::json& j, const ProForma::Label::Type& lt)
  {
    switch (lt)
    {
      case ProForma::Label::Type::CROSSLINK: j = "CROSSLINK"; break;
      case ProForma::Label::Type::BRANCH:    j = "BRANCH"; break;
      case ProForma::Label::Type::AMBIGUOUS: j = "AMBIGUOUS"; break;
    }
  }

  /// Construct ProForma::Label::Type enum from JSON string
  inline void from_json(const nlohmann::json& j, ProForma::Label::Type& lt)
  {
    std::string s = j.get<std::string>();
    if (s == "CROSSLINK") lt = ProForma::Label::Type::CROSSLINK;
    else if (s == "BRANCH") lt = ProForma::Label::Type::BRANCH;
    else if (s == "AMBIGUOUS") lt = ProForma::Label::Type::AMBIGUOUS;
    else throw std::invalid_argument("Unknown ProForma::Label::Type: " + s);
  }

  /// @}


  /// @name Label JSON serialization
  /// @{

  /// Convert Label to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::Label& lbl)
  {
    j = nlohmann::json{
      {"type", lbl.type},
      {"identifier", static_cast<std::string>(lbl.identifier)}
    };
    if (lbl.score.has_value())
    {
      j["score"] = lbl.score.value();
    }
  }

  /// Construct Label from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::Label& lbl)
  {
    j.at("type").get_to(lbl.type);
    lbl.identifier = j.at("identifier").get<std::string>();
    if (j.contains("score") && !j.at("score").is_null())
    {
      lbl.score = j.at("score").get<double>();
    }
    else
    {
      lbl.score = std::nullopt;
    }
  }

  /// @}


  /// @name Modification JSON serialization
  /// @{

  /// Convert Modification to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::Modification& mod)
  {
    j = nlohmann::json::array();
    for (const auto& [tag, label] : mod.alternatives)
    {
      nlohmann::json alt;
      alt["tag"] = tag;
      if (label.has_value())
      {
        alt["label"] = label.value();
      }
      j.push_back(alt);
    }
  }

  /// Construct Modification from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::Modification& mod)
  {
    mod.alternatives.clear();
    mod.resolved_mod = nullptr;  // Reset to avoid stale pointers
    for (const auto& item : j)
    {
      ProForma::ModificationTag tag = item.at("tag").get<ProForma::ModificationTag>();
      std::optional<ProForma::Label> label;
      if (item.contains("label") && !item.at("label").is_null())
      {
        label = item.at("label").get<ProForma::Label>();
      }
      mod.alternatives.emplace_back(tag, label);
    }
  }

  /// @}


  /// @name SequenceElement JSON serialization
  /// @{

  /// Convert SequenceElement to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::SequenceElement& se)
  {
    j = nlohmann::json{
      {"amino_acid", std::string(1, se.amino_acid)},
      {"modifications", se.modifications}
    };
  }

  /// Construct SequenceElement from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::SequenceElement& se)
  {
    std::string aa = j.at("amino_acid").get<std::string>();
    if (aa.size() != 1)
    {
      throw std::invalid_argument("amino_acid must be exactly 1 character, got length " + std::to_string(aa.size()));
    }
    se.amino_acid = aa[0];
    j.at("modifications").get_to(se.modifications);
  }

  /// @}


  /// @name AmbiguousRegion JSON serialization
  /// @{

  /// Convert AmbiguousRegion to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::AmbiguousRegion& ar)
  {
    j = nlohmann::json{{"elements", ar.elements}};
  }

  /// Construct AmbiguousRegion from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::AmbiguousRegion& ar)
  {
    j.at("elements").get_to(ar.elements);
  }

  /// @}


  /// @name ModifiedRange JSON serialization
  /// @{

  /// Convert ModifiedRange to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::ModifiedRange& mr)
  {
    j = nlohmann::json{
      {"elements", mr.elements},
      {"modifications", mr.modifications}
    };
  }

  /// Construct ModifiedRange from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::ModifiedRange& mr)
  {
    j.at("elements").get_to(mr.elements);
    j.at("modifications").get_to(mr.modifications);
  }

  /// @}


  /// @name SequenceSection variant JSON serialization
  /// @{

  /// Convert SequenceSection variant to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::SequenceSection& ss)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, ProForma::SequenceElement>)
      {
        j = nlohmann::json{{"type", "element"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::AmbiguousRegion>)
      {
        j = nlohmann::json{{"type", "ambiguous_region"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::ModifiedRange>)
      {
        j = nlohmann::json{{"type", "modified_range"}, {"value", arg}};
      }
    }, ss);
  }

  /// Construct SequenceSection variant from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::SequenceSection& ss)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "element")
    {
      ss = j.at("value").get<ProForma::SequenceElement>();
    }
    else if (type == "ambiguous_region")
    {
      ss = j.at("value").get<ProForma::AmbiguousRegion>();
    }
    else if (type == "modified_range")
    {
      ss = j.at("value").get<ProForma::ModifiedRange>();
    }
    else
    {
      throw std::invalid_argument("Unknown ProForma::SequenceSection type: " + type);
    }
  }

  /// @}


  /// @name UnlocalisedMod JSON serialization
  /// @{

  /// Convert UnlocalisedMod to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::UnlocalisedMod& um)
  {
    j = nlohmann::json{{"modifications", um.modifications}};
    if (um.occurrence.has_value())
    {
      j["occurrence"] = um.occurrence.value();
    }
  }

  /// Construct UnlocalisedMod from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::UnlocalisedMod& um)
  {
    j.at("modifications").get_to(um.modifications);
    if (j.contains("occurrence") && !j.at("occurrence").is_null())
    {
      um.occurrence = j.at("occurrence").get<int>();
    }
    else
    {
      um.occurrence = std::nullopt;
    }
  }

  /// @}


  /// @name LabileModification JSON serialization
  /// @{

  /// Convert LabileModification to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::LabileModification& lm)
  {
    j = nlohmann::json{{"modification", lm.modification}};
  }

  /// Construct LabileModification from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::LabileModification& lm)
  {
    j.at("modification").get_to(lm.modification);
  }

  /// @}


  /// @name GlobalModification JSON serialization
  /// @{

  /// Convert GlobalModification to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::GlobalModification& gm)
  {
    std::vector<std::string> locs;
    for (const auto& loc : gm.locations)
    {
      locs.push_back(static_cast<std::string>(loc));
    }
    j = nlohmann::json{
      {"modification", gm.modification},
      {"locations", locs}
    };
  }

  /// Construct GlobalModification from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::GlobalModification& gm)
  {
    j.at("modification").get_to(gm.modification);
    gm.locations.clear();
    for (const auto& loc : j.at("locations"))
    {
      gm.locations.push_back(String(loc.get<std::string>()));
    }
  }

  /// @}


  /// @name IsotopeReplacement JSON serialization
  /// @{

  /// Convert IsotopeReplacement to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::IsotopeReplacement& ir)
  {
    j = nlohmann::json{{"isotope", static_cast<std::string>(ir.isotope)}};
  }

  /// Construct IsotopeReplacement from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::IsotopeReplacement& ir)
  {
    ir.isotope = j.at("isotope").get<std::string>();
  }

  /// @}


  /// @name GlobalModEntry variant JSON serialization
  /// @{

  /// Convert GlobalModEntry variant to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::GlobalModEntry& gme)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, ProForma::IsotopeReplacement>)
      {
        j = nlohmann::json{{"type", "isotope_replacement"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ProForma::GlobalModification>)
      {
        j = nlohmann::json{{"type", "global_modification"}, {"value", arg}};
      }
    }, gme);
  }

  /// Construct GlobalModEntry variant from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::GlobalModEntry& gme)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "isotope_replacement")
    {
      gme = j.at("value").get<ProForma::IsotopeReplacement>();
    }
    else if (type == "global_modification")
    {
      gme = j.at("value").get<ProForma::GlobalModification>();
    }
    else
    {
      throw std::invalid_argument("Unknown ProForma::GlobalModEntry type: " + type);
    }
  }

  /// @}


  /// @name AdductIon JSON serialization
  /// @{

  /// Convert AdductIon to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::AdductIon& ai)
  {
    j = nlohmann::json{
      {"formula", static_cast<std::string>(ai.formula)},
      {"charge", ai.charge}
    };
    if (ai.occurrence.has_value())
    {
      j["occurrence"] = ai.occurrence.value();
    }
  }

  /// Construct AdductIon from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::AdductIon& ai)
  {
    ai.formula = j.at("formula").get<std::string>();
    j.at("charge").get_to(ai.charge);
    if (j.contains("occurrence") && !j.at("occurrence").is_null())
    {
      ai.occurrence = j.at("occurrence").get<int>();
    }
    else
    {
      ai.occurrence = std::nullopt;
    }
  }

  /// @}


  /// @name ChargeState variant JSON serialization
  /// @{

  /// Convert ChargeState variant to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::ChargeState& cs)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, int>)
      {
        j = nlohmann::json{{"type", "simple"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, std::vector<ProForma::AdductIon>>)
      {
        j = nlohmann::json{{"type", "adducts"}, {"value", arg}};
      }
    }, cs);
  }

  /// Construct ChargeState variant from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::ChargeState& cs)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "simple")
    {
      cs = j.at("value").get<int>();
    }
    else if (type == "adducts")
    {
      cs = j.at("value").get<std::vector<ProForma::AdductIon>>();
    }
    else
    {
      throw std::invalid_argument("Unknown ProForma::ChargeState type: " + type);
    }
  }

  /// @}


  /// @name Peptidoform JSON serialization
  /// @{

  /// Convert Peptidoform to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::Peptidoform& pf)
  {
    j = nlohmann::json{
      {"global_mods", pf.global_mods},
      {"unlocalised_mods", pf.unlocalised_mods},
      {"labile_mods", pf.labile_mods},
      {"n_term_mods", pf.n_term_mods},
      {"sequence", pf.sequence},
      {"c_term_mods", pf.c_term_mods}
    };
    if (pf.name.has_value())
    {
      j["name"] = static_cast<std::string>(pf.name.value());
    }
    if (pf.charge.has_value())
    {
      j["charge"] = pf.charge.value();
    }
  }

  /// Construct Peptidoform from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::Peptidoform& pf)
  {
    if (j.contains("global_mods"))
    {
      j.at("global_mods").get_to(pf.global_mods);
    }
    else
    {
      pf.global_mods.clear();  // Clear when field is omitted
    }
    j.at("unlocalised_mods").get_to(pf.unlocalised_mods);
    j.at("labile_mods").get_to(pf.labile_mods);
    j.at("n_term_mods").get_to(pf.n_term_mods);
    j.at("sequence").get_to(pf.sequence);
    j.at("c_term_mods").get_to(pf.c_term_mods);
    if (j.contains("name") && !j.at("name").is_null())
    {
      pf.name = String(j.at("name").get<std::string>());
    }
    else
    {
      pf.name = std::nullopt;
    }
    if (j.contains("charge") && !j.at("charge").is_null())
    {
      pf.charge = j.at("charge").get<ProForma::ChargeState>();
    }
    else
    {
      pf.charge = std::nullopt;
    }
  }

  /// @}


  /// @name PeptidoformIon JSON serialization
  /// @{

  /// Convert PeptidoformIon to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::PeptidoformIon& pfi)
  {
    j = nlohmann::json{{"chains", pfi.chains}, {"is_chimeric", pfi.is_chimeric}};
    if (pfi.name.has_value())
    {
      j["name"] = static_cast<std::string>(pfi.name.value());
    }
    if (pfi.charge.has_value())
    {
      j["charge"] = pfi.charge.value();
    }
  }

  /// Construct PeptidoformIon from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::PeptidoformIon& pfi)
  {
    j.at("chains").get_to(pfi.chains);
    if (j.contains("name") && !j.at("name").is_null())
    {
      pfi.name = String(j.at("name").get<std::string>());
    }
    else
    {
      pfi.name = std::nullopt;
    }
    if (j.contains("charge") && !j.at("charge").is_null())
    {
      pfi.charge = j.at("charge").get<ProForma::ChargeState>();
    }
    else
    {
      pfi.charge = std::nullopt;
    }
    if (j.contains("is_chimeric"))
    {
      pfi.is_chimeric = j.at("is_chimeric").get<bool>();
    }
    else
    {
      pfi.is_chimeric = false;
    }
  }

  /// @}


  /// @name CrossLinkGroup JSON serialization
  /// @{

  /// Convert CrossLinkGroup to JSON object
  inline void to_json(nlohmann::json& j, const ProForma::CrossLinkGroup& clg)
  {
    nlohmann::json sites_json = nlohmann::json::array();
    for (const auto& [chain_idx, site_idx] : clg.sites)
    {
      sites_json.push_back(nlohmann::json{{"chain_index", chain_idx}, {"site_index", site_idx}});
    }
    j = nlohmann::json{
      {"label", static_cast<std::string>(clg.label)},
      {"sites", sites_json}
    };
  }

  /// Construct CrossLinkGroup from JSON object
  inline void from_json(const nlohmann::json& j, ProForma::CrossLinkGroup& clg)
  {
    clg.label = j.at("label").get<std::string>();
    clg.sites.clear();
    for (const auto& site : j.at("sites"))
    {
      size_t chain_idx = site.at("chain_index").get<size_t>();
      size_t site_idx = site.at("site_index").get<size_t>();
      clg.sites.emplace_back(chain_idx, site_idx);
    }
  }

  /// @}


  // Note: The convenience functions toJSON(), peptidoformFromJSON(), peptidoformIonFromJSON()
  // are declared in ProFormaData.h and implemented in ProFormaDataJson.cpp.
  // This header provides the nlohmann::json ADL overloads needed by those implementations.

} // namespace OpenMS
