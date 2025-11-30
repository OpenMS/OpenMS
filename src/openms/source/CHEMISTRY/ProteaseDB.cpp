// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <fstream>
using namespace std;

namespace OpenMS
{
  ProteaseDB::ProteaseDB():
    DigestionEnzymeDB<DigestionEnzymeProtein, ProteaseDB>()  // no file - we'll add built-in enzymes first
  {
    // Add built-in enzymes
    addBuiltInEnzymes_();
    
    // Try to load from file if present (allows user customization/extension)
    readEnzymesFromFileIfPresent_("CHEMISTRY/Enzymes.xml");
  }

  void ProteaseDB::addBuiltInEnzymes_()
  {
    // Built-in enzyme definitions from share/OpenMS/CHEMISTRY/Enzymes.xml
    // Empty EmpiricalFormula and empty string parameters match the original XML
    // where certain values were not specified for some enzymes.

    // Trypsin
    addEnzyme_(new DigestionEnzymeProtein(
      "Trypsin",
      "(?<=[KRX])(?!P)",
      set<String>(),
      "Trypsin cleaves following a K or R residue unless the next residue is P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001251",
      "[KR]|{P}",
      1,   // CometID
      -1,  // MSGFID
      0    // OMSSAID
    ));

    // Arg-C
    addEnzyme_(new DigestionEnzymeProtein(
      "Arg-C",
      "(?<=[RX])(?!P)",
      set<String>{"Clostripain", "argc", "arg_c"},
      "Arg-C cleaves following R residue unless the next residue is P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001303",
      "[R]|{P}",
      5,   // CometID
      -1,  // MSGFID
      1    // OMSSAID
    ));

    // Arg-C/P
    addEnzyme_(new DigestionEnzymeProtein(
      "Arg-C/P",
      "(?<=[RX])",
      set<String>(),
      "Arg-C/P cleaves after R residues.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "",
      "[R]|[X]",
      12,  // CometID
      6,   // MSGFID
      -1   // OMSSAID
    ));

    // Asp-N
    addEnzyme_(new DigestionEnzymeProtein(
      "Asp-N",
      "(?=[DBX])",
      set<String>{"asp_n"},
      "Asp-N cleaves before D(or B).",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001304",
      "[X]|[BD]",
      6,   // CometID
      -1,  // MSGFID
      12   // OMSSAID
    ));

    // Asp-N/B
    addEnzyme_(new DigestionEnzymeProtein(
      "Asp-N/B",
      "(?=[DX])",
      set<String>(),
      "Asp-N/B cleaves before D(while B is ignored).",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "",
      "[X]|[D]",
      16,  // CometID
      7,   // MSGFID
      -1   // OMSSAID
    ));

    // Asp-N_ambic
    addEnzyme_(new DigestionEnzymeProtein(
      "Asp-N_ambic",
      "(?=[DBEZX])",
      set<String>(),
      "Asp-N Ammonium bicarbonate cleaves before D(or B) or E(or Z).",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001305",
      "[X]|[DE]",
      17,  // CometID
      -1,  // MSGFID
      19   // OMSSAID
    ));

    // Chymotrypsin
    addEnzyme_(new DigestionEnzymeProtein(
      "Chymotrypsin",
      "(?<=[FYWLJX])(?!P)",
      set<String>(),
      "Chymotrypsin cleaves following F, Y, W or L(or J) residue unless the next residue is P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001306",
      "[FYWL]|{P}",
      10,  // CometID
      -1,  // MSGFID
      3    // OMSSAID
    ));

    // Chymotrypsin/P
    addEnzyme_(new DigestionEnzymeProtein(
      "Chymotrypsin/P",
      "(?<=[FYWLJX])",
      set<String>(),
      "Chymotrypsin cleaves following F, Y, W or L(or J) residue.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "",
      "[FYWL]|[X]",
      15,  // CometID
      2,   // MSGFID
      -1   // OMSSAID
    ));

    // CNBr
    addEnzyme_(new DigestionEnzymeProtein(
      "CNBr",
      "(?<=[MX])",
      set<String>(),
      "CNBr cleaves following M.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001307",
      "[M]|[X]",
      7,   // CometID
      -1,  // MSGFID
      2    // OMSSAID
    ));

    // Formic_acid
    addEnzyme_(new DigestionEnzymeProtein(
      "Formic_acid",
      "(?<=[DBX])(?=[DBX])",
      set<String>(),
      "Formic_acid cuts after D(or B) and next residue is D (or B).",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001308",
      "[D]|[D]",
      18,  // CometID
      -1,  // MSGFID
      4    // OMSSAID
    ));

    // Lys-C
    addEnzyme_(new DigestionEnzymeProtein(
      "Lys-C",
      "(?<=[KX])(?!P)",
      set<String>{"lys_c"},
      "Lys-C cuts after K if not followed by P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001309",
      "[K]|{P}",
      3,   // CometID
      -1,  // MSGFID
      5    // OMSSAID
    ));

    // Lys-N
    addEnzyme_(new DigestionEnzymeProtein(
      "Lys-N",
      "(?=[KX])",
      set<String>{"lys_n"},
      "Lys-N cuts before K.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "",
      "[X]|[K]",
      4,   // CometID
      4,   // MSGFID
      -1   // OMSSAID
    ));

    // Lys-C/P
    addEnzyme_(new DigestionEnzymeProtein(
      "Lys-C/P",
      "(?<=[KX])",
      set<String>(),
      "Lys-C/P cuts after K.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001310",
      "[K]|[X]",
      13,  // CometID
      3,   // MSGFID
      6    // OMSSAID
    ));

    // PepsinA
    addEnzyme_(new DigestionEnzymeProtein(
      "PepsinA",
      "(?<=[FLJX])",
      set<String>(),
      "PepsinA cuts after F or L(or J).",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001311",
      "[FL]|[X]",
      9,   // CometID
      -1,  // MSGFID
      7    // OMSSAID
    ));

    // TrypChymo
    addEnzyme_(new DigestionEnzymeProtein(
      "TrypChymo",
      "(?<=[FYWLJKRX])(?!P)",
      set<String>(),
      "TrypChymo cuts after F, Y, W, L(or J), K or R if not followed by P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001312",
      "[FYWLKR]|{P}",
      19,  // CometID
      -1,  // MSGFID
      9    // OMSSAID
    ));

    // Trypsin/P
    addEnzyme_(new DigestionEnzymeProtein(
      "Trypsin/P",
      "(?<=[KRX])",
      set<String>(),
      "Trypsin/P cuts after K or R.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001313",
      "[KR]|[X]",
      2,   // CometID
      1,   // MSGFID
      10   // OMSSAID
    ));

    // V8-DE
    addEnzyme_(new DigestionEnzymeProtein(
      "V8-DE",
      "(?<=[DBEZX])(?!P)",
      set<String>(),
      "V8-DE cuts after D(or B) or E(or Z) if not followed by P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001314",
      "[BDEZ]|{P}",
      20,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // V8-E
    addEnzyme_(new DigestionEnzymeProtein(
      "V8-E",
      "(?<=[EZX])(?!P)",
      set<String>(),
      "V8-E cuts after E(or Z) if not followed by P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001315",
      "[EZ]|{P}",
      21,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // leukocyte elastase
    addEnzyme_(new DigestionEnzymeProtein(
      "leukocyte elastase",
      "(?<=[ALIJVX])(?!P)",
      set<String>(),
      "leukocyte elastase cuts after A or L or I(or J) or V if not followed by P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001915",
      "[ALIV]|{P}",
      14,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // proline endopeptidase
    // Note: complex regex pattern that matches after H, K, or R followed by P, unless next residue is P
    // XTandemID is intentionally empty as this enzyme is not supported by X!Tandem
    addEnzyme_(new DigestionEnzymeProtein(
      "proline endopeptidase",
      "(?<=[HKRX][PX])(?!P)",
      set<String>(),
      "proline endopeptidase cuts after HP, KP or RP if not followed by P.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001916",
      "",  // No XTandemID available
      22,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // glutamyl endopeptidase
    addEnzyme_(new DigestionEnzymeProtein(
      "glutamyl endopeptidase",
      "(?<=[DBEZX])",
      set<String>{"Glu-C", "glu_c", "staphylococcal protease"},
      "glutamyl endopeptidase cuts after D(or B) or E(or Z).",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001917",
      "[DE]|[X]",
      8,   // CometID
      5,   // MSGFID
      20   // OMSSAID
    ));

    // Alpha-lytic protease
    addEnzyme_(new DigestionEnzymeProtein(
      "Alpha-lytic protease",
      "(?<=[TASVX])",
      set<String>(),
      "Alpha-lytic protease (aLP) cuts after T, A, S, or V.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "",
      "[ASTV]|[X]",
      23,  // CometID
      8,   // MSGFID
      -1   // OMSSAID
    ));

    // 2-iodobenzoate
    addEnzyme_(new DigestionEnzymeProtein(
      "2-iodobenzoate",
      "(?<=[WX])",
      set<String>(),
      "2-iodobenzoate cuts after W.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001918",
      "[W]|[X]",
      24,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // iodosobenzoate
    addEnzyme_(new DigestionEnzymeProtein(
      "iodosobenzoate",
      "(?<=W)",
      set<String>(),
      "iodosobenzoate cuts after W.",
      EmpiricalFormula(""),
      EmpiricalFormula(""),
      "",
      "",
      25,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // staphylococcal protease/D
    addEnzyme_(new DigestionEnzymeProtein(
      "staphylococcal protease/D",
      "(?<=[EZX])",
      set<String>{"Glu-C/D"},
      "staphylococcal protease/D cuts after E(or Z).",
      EmpiricalFormula(""),
      EmpiricalFormula(""),
      "",
      "",
      26,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // proline-endopeptidase/HKR
    addEnzyme_(new DigestionEnzymeProtein(
      "proline-endopeptidase/HKR",
      "(?<=[PX])",
      set<String>(),
      "proline-endopeptidase/HKR cuts after P.",
      EmpiricalFormula(""),
      EmpiricalFormula(""),
      "",
      "",
      27,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // Glu-C+P
    addEnzyme_(new DigestionEnzymeProtein(
      "Glu-C+P",
      "(?<=[DBEZX])(?!P)",
      set<String>{"staphylococcal protease+P"},
      "Glu-C+P cuts after D(or B) or E(or Z) unless followed by P.",
      EmpiricalFormula(""),
      EmpiricalFormula(""),
      "",
      "",
      28,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // PepsinA + P
    addEnzyme_(new DigestionEnzymeProtein(
      "PepsinA + P",
      "(?<=[FLJX])(?!P)",
      set<String>(),
      "PepsinA + P cuts after F or L(or J) unless followed by P.",
      EmpiricalFormula(""),
      EmpiricalFormula(""),
      "",
      "",
      29,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // cyanogen-bromide
    addEnzyme_(new DigestionEnzymeProtein(
      "cyanogen-bromide",
      "(?<=[MX])",
      set<String>(),
      "cyanogen-bromide cuts after M.",
      EmpiricalFormula(""),
      EmpiricalFormula(""),
      "",
      "",
      30,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // Clostripain/P
    addEnzyme_(new DigestionEnzymeProtein(
      "Clostripain/P",
      "(?<=[RX])",
      set<String>(),
      "Clostripain/P cuts after R.",
      EmpiricalFormula(""),
      EmpiricalFormula(""),
      "",
      "",
      31,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // elastase-trypsin-chymotrypsin
    addEnzyme_(new DigestionEnzymeProtein(
      "elastase-trypsin-chymotrypsin",
      "(?<=[ALIVKRWFYX])(?!P)",
      set<String>(),
      "elastase-trypsin-chymotrypsin cuts after A,L,I(or J),V,K,R,W,F,Y unless followed by P.",
      EmpiricalFormula(""),
      EmpiricalFormula(""),
      "",
      "",
      32,  // CometID
      -1,  // MSGFID
      -1   // OMSSAID
    ));

    // no cleavage
    addEnzyme_(new DigestionEnzymeProtein(
      "no cleavage",
      "()",
      set<String>{"no_cut"},
      "no cleavage.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001955",
      "[Z]|[Z]",
      11,  // CometID
      9,   // MSGFID
      11   // OMSSAID
    ));

    // unspecific cleavage
    addEnzyme_(new DigestionEnzymeProtein(
      "unspecific cleavage",
      "(?<=[A-Z])",
      set<String>{"no_enzyme", "nonspecific"},
      "unspecific cleavage cuts at every site.",
      EmpiricalFormula("H"),
      EmpiricalFormula("OH"),
      "MS:1001956",
      "[X]|[X]",
      0,   // CometID
      0,   // MSGFID
      17   // OMSSAID
    ));
  }

  void ProteaseDB::getAllXTandemNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if (!(*it)->getXTandemID().empty())
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllCometNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getCometID() != -1)
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllOMSSANames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getOMSSAID() != -1)
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllMSGFNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getMSGFID() != -1) // MS-GF+ starts enzyme numbering at 0
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::writeTSV(String const& filename)
  {
    std::ofstream ofs(filename, std::ofstream::out);
    ofs << "OpenMS_AllowedEnzymes" << "\n";
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      ofs << (*it)->getName() << "\n";
    }
  }
} // namespace OpenMS
