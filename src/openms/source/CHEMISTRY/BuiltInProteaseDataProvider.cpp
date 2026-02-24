// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sachsenberg $
// $Authors: Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/BuiltInProteaseDataProvider.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace std;

namespace OpenMS
{
  vector<unique_ptr<DigestionEnzymeProtein>> BuiltInProteaseDataProvider::loadEnzymes()
  {
    vector<unique_ptr<DigestionEnzymeProtein>> enzymes;

    // Trypsin
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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
    enzymes.push_back(make_unique<DigestionEnzymeProtein>(
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

    return enzymes;
  }

} // namespace OpenMS
