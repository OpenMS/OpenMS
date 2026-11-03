// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/MzPAF.h>
#include <OpenMS/CHEMISTRY/TheoreticalGlycanSpectrumGenerator.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/test_config.h>
#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>

using namespace OpenMS;
using Generator = TheoreticalGlycanSpectrumGenerator;
using Composition = Generator::Composition;
using Ion = Generator::IonType;

namespace
{
Composition comp(std::initializer_list<std::pair<std::string, int>> counts)
{
  Composition result;
  for (const auto& [symbol, count] : counts)
  {
    result.components.emplace_back(symbol, count);
  }
  return result;
}

const Generator::Fragment& findIon(const std::vector<Generator::Fragment>& fragments, std::string_view name, Int charge = 1)
{
  const auto it = std::find_if(fragments.begin(), fragments.end(), [&](const auto& f) { return f.name == name && f.charge == charge; });
  if (it == fragments.end()) { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing ion: " + std::string(name)); }
  return *it;
}

using GlypyKey = std::tuple<Ion, std::optional<Size>, std::vector<Size>>;
struct GlypyCase
{
  GlycanStructure tree;
  std::map<GlypyKey, double> masses;
};

struct GlypyNodeRow
{
  const char* fixture;
  Size index;
  Int parent;
  const char* symbol;
};

struct GlypyFragmentRow
{
  const char* fixture;
  char series;
  Int root;
  const char* branches;
  double mass;
};

// Reference glycan trees and fragment masses computed with GlyPy 1.0.17
// (https://github.com/mobiusklein/glypy, revision 8d129a8c950e8635165cda9b4b1af392d4e7289e).
// The three structures are GlyPy's own test glycans, read from tests/common.py at that
// revision as GlycoCT and converted to OpenMS composition symbols; G58143RL is a GlyTouCan
// accession. Masses are GlyPy's, verbatim to 11 decimal places. To re-derive them, load the
// same three GlycoCT strings with glypy==1.0.17 and enumerate glycosidic cleavages.
//
// Nodes: fixture, node index, parent (-1 = root), residue symbol or Formula:<formula>.
// clang-format off: one line group per fixture.
const std::vector<GlypyNodeRow> GLYPY_NODES {
  {"common_glycan", 0, -1, "Hex"}, {"common_glycan", 1, 0, "Hex"}, {"common_glycan", 2, 1, "HexNAc"}, {"common_glycan", 3, 2, "Fuc"},
  {"common_glycan", 4, 2, "Hex"}, {"common_glycan", 5, 4, "HexNAc"}, {"common_glycan", 6, 5, "Fuc"}, {"common_glycan", 7, 5, "Hex"},
  {"branchy_glycan", 0, -1, "HexNAc"}, {"branchy_glycan", 1, 0, "Hex"}, {"branchy_glycan", 2, 1, "Hex"}, {"branchy_glycan", 3, 2, "HexNAc"},
  {"branchy_glycan", 4, 3, "Hex"}, {"branchy_glycan", 5, 1, "Hex"}, {"branchy_glycan", 6, 5, "HexNAc"}, {"branchy_glycan", 7, 6, "Hex"},
  {"branchy_glycan", 8, 5, "HexNAc"}, {"branchy_glycan", 9, 8, "Hex"},
  {"G58143RL", 0, -1, "HexNAc"}, {"G58143RL", 1, 0, "HexNAc"}, {"G58143RL", 2, 1, "Hex"}, {"G58143RL", 3, 2, "Hex"}, {"G58143RL", 4, 3, "HexNAc"},
  {"G58143RL", 5, 4, "Hex"}, {"G58143RL", 6, 5, "Neu5Ac"}, {"G58143RL", 7, 2, "Hex"}, {"G58143RL", 8, 7, "Formula:C8H13N1O8S1"},
  {"G58143RL", 9, 8, "Hex"}, {"G58143RL", 10, 9, "Neu5Ac"}, {"G58143RL", 11, 0, "Fuc"},
};

// Fragments: fixture, series, root cleavage (-1 = none), branch cleavages ("-" = none), neutral mass (Da).
const std::vector<GlypyFragmentRow> GLYPY_FRAGMENTS {
  {"common_glycan", 'B', 1, "-", 1184.43303289240}, {"common_glycan", 'B', 1, "2", 162.05282341850}, {"common_glycan", 'B', 1, "3", 1038.37512409346},
  {"common_glycan", 'B', 1, "3,4", 365.13219593801}, {"common_glycan", 'B', 1, "3,5", 527.18501935651},
  {"common_glycan", 'B', 1, "3,6", 892.31721529452}, {"common_glycan", 'B', 1, "3,7", 876.32230067496},
  {"common_glycan", 'B', 1, "4", 511.19010473695}, {"common_glycan", 'B', 1, "5", 673.24292815545}, {"common_glycan", 'B', 1, "6", 1038.37512409346},
  {"common_glycan", 'B', 1, "6,7", 876.32230067496}, {"common_glycan", 'B', 1, "7", 1022.38020947390},
  {"common_glycan", 'B', 2, "-", 1022.38020947390}, {"common_glycan", 'B', 2, "3", 876.32230067496},
  {"common_glycan", 'B', 2, "3,4", 203.07937251951}, {"common_glycan", 'B', 2, "3,5", 365.13219593801},
  {"common_glycan", 'B', 2, "3,6", 730.26439187602}, {"common_glycan", 'B', 2, "3,7", 714.26947725646},
  {"common_glycan", 'B', 2, "4", 349.13728131845}, {"common_glycan", 'B', 2, "5", 511.19010473695}, {"common_glycan", 'B', 2, "6", 876.32230067496},
  {"common_glycan", 'B', 2, "6,7", 714.26947725646}, {"common_glycan", 'B', 2, "7", 860.32738605540}, {"common_glycan", 'B', 3, "-", 146.05790879894},
  {"common_glycan", 'B', 4, "-", 673.24292815545}, {"common_glycan", 'B', 4, "5", 162.05282341850}, {"common_glycan", 'B', 4, "6", 527.18501935651},
  {"common_glycan", 'B', 4, "6,7", 365.13219593801}, {"common_glycan", 'B', 4, "7", 511.19010473695}, {"common_glycan", 'B', 5, "-", 511.19010473695},
  {"common_glycan", 'B', 5, "6", 365.13219593801}, {"common_glycan", 'B', 5, "6,7", 203.07937251951}, {"common_glycan", 'B', 5, "7", 349.13728131845},
  {"common_glycan", 'B', 6, "-", 146.05790879894}, {"common_glycan", 'B', 7, "-", 162.05282341850}, {"common_glycan", 'C', 1, "-", 1202.44359757610},
  {"common_glycan", 'C', 1, "2", 180.06338810220}, {"common_glycan", 'C', 1, "3", 1056.38568877716},
  {"common_glycan", 'C', 1, "3,4", 383.14276062171}, {"common_glycan", 'C', 1, "3,5", 545.19558404021},
  {"common_glycan", 'C', 1, "3,6", 910.32777997822}, {"common_glycan", 'C', 1, "3,7", 894.33286535866},
  {"common_glycan", 'C', 1, "4", 529.20066942065}, {"common_glycan", 'C', 1, "5", 691.25349283915}, {"common_glycan", 'C', 1, "6", 1056.38568877716},
  {"common_glycan", 'C', 1, "6,7", 894.33286535866}, {"common_glycan", 'C', 1, "7", 1040.39077415760},
  {"common_glycan", 'C', 2, "-", 1040.39077415760}, {"common_glycan", 'C', 2, "3", 894.33286535866},
  {"common_glycan", 'C', 2, "3,4", 221.08993720321}, {"common_glycan", 'C', 2, "3,5", 383.14276062171},
  {"common_glycan", 'C', 2, "3,6", 748.27495655972}, {"common_glycan", 'C', 2, "3,7", 732.28004194016},
  {"common_glycan", 'C', 2, "4", 367.14784600215}, {"common_glycan", 'C', 2, "5", 529.20066942065}, {"common_glycan", 'C', 2, "6", 894.33286535866},
  {"common_glycan", 'C', 2, "6,7", 732.28004194016}, {"common_glycan", 'C', 2, "7", 878.33795073910}, {"common_glycan", 'C', 3, "-", 164.06847348264},
  {"common_glycan", 'C', 4, "-", 691.25349283915}, {"common_glycan", 'C', 4, "5", 180.06338810220}, {"common_glycan", 'C', 4, "6", 545.19558404021},
  {"common_glycan", 'C', 4, "6,7", 383.14276062171}, {"common_glycan", 'C', 4, "7", 529.20066942065}, {"common_glycan", 'C', 5, "-", 529.20066942065},
  {"common_glycan", 'C', 5, "6", 383.14276062171}, {"common_glycan", 'C', 5, "6,7", 221.08993720321}, {"common_glycan", 'C', 5, "7", 367.14784600215},
  {"common_glycan", 'C', 6, "-", 164.06847348264}, {"common_glycan", 'C', 7, "-", 180.06338810220}, {"common_glycan", 'Y', -1, "1", 180.06338810220},
  {"common_glycan", 'Y', -1, "2", 342.11621152070}, {"common_glycan", 'Y', -1, "3", 1218.43851219566},
  {"common_glycan", 'Y', -1, "3,4", 545.19558404021}, {"common_glycan", 'Y', -1, "3,5", 707.24840745871},
  {"common_glycan", 'Y', -1, "3,6", 1072.38060339672}, {"common_glycan", 'Y', -1, "3,6,7", 910.32777997822},
  {"common_glycan", 'Y', -1, "3,7", 1056.38568877716}, {"common_glycan", 'Y', -1, "4", 691.25349283915},
  {"common_glycan", 'Y', -1, "5", 853.30631625765}, {"common_glycan", 'Y', -1, "6", 1218.43851219566},
  {"common_glycan", 'Y', -1, "6,7", 1056.38568877716}, {"common_glycan", 'Y', -1, "7", 1202.44359757610},
  {"common_glycan", 'Z', -1, "1", 162.05282341850}, {"common_glycan", 'Z', -1, "2", 324.10564683700},
  {"common_glycan", 'Z', -1, "3", 1200.42794751196}, {"common_glycan", 'Z', -1, "3,4", 509.17445467281},
  {"common_glycan", 'Z', -1, "3,5", 671.22727809131}, {"common_glycan", 'Z', -1, "3,6", 1036.35947402932},
  {"common_glycan", 'Z', -1, "3,6,7", 856.29608592712}, {"common_glycan", 'Z', -1, "3,7", 1020.36455940976},
  {"common_glycan", 'Z', -1, "4", 673.24292815545}, {"common_glycan", 'Z', -1, "5", 835.29575157395},
  {"common_glycan", 'Z', -1, "6", 1200.42794751196}, {"common_glycan", 'Z', -1, "6,7", 1020.36455940976},
  {"common_glycan", 'Z', -1, "7", 1184.43303289240},
  {"branchy_glycan", 'B', 1, "-", 1581.55505806953}, {"branchy_glycan", 'B', 1, "2", 1054.37003871302},
  {"branchy_glycan", 'B', 1, "2,5", 162.05282341850}, {"branchy_glycan", 'B', 1, "2,6", 689.23784277501},
  {"branchy_glycan", 'B', 1, "2,7", 892.31721529452}, {"branchy_glycan", 'B', 1, "2,8", 689.23784277501},
  {"branchy_glycan", 'B', 1, "2,9", 892.31721529452}, {"branchy_glycan", 'B', 1, "3", 1216.42286213152},
  {"branchy_glycan", 'B', 1, "3,5", 324.10564683700}, {"branchy_glycan", 'B', 1, "3,6", 851.29066619351},
  {"branchy_glycan", 'B', 1, "3,7", 1054.37003871302}, {"branchy_glycan", 'B', 1, "3,8", 851.29066619351},
  {"branchy_glycan", 'B', 1, "3,9", 1054.37003871302}, {"branchy_glycan", 'B', 1, "4", 1419.50223465103},
  {"branchy_glycan", 'B', 1, "4,5", 527.18501935651}, {"branchy_glycan", 'B', 1, "4,6", 1054.37003871302},
  {"branchy_glycan", 'B', 1, "4,7", 1257.44941123253}, {"branchy_glycan", 'B', 1, "4,8", 1054.37003871302},
  {"branchy_glycan", 'B', 1, "4,9", 1257.44941123253}, {"branchy_glycan", 'B', 1, "5", 689.23784277501},
  {"branchy_glycan", 'B', 1, "6", 1216.42286213152}, {"branchy_glycan", 'B', 1, "6,8", 851.29066619351},
  {"branchy_glycan", 'B', 1, "6,9", 1054.37003871302}, {"branchy_glycan", 'B', 1, "7", 1419.50223465103},
  {"branchy_glycan", 'B', 1, "7,8", 1054.37003871302}, {"branchy_glycan", 'B', 1, "7,9", 1257.44941123253},
  {"branchy_glycan", 'B', 1, "8", 1216.42286213152}, {"branchy_glycan", 'B', 1, "9", 1419.50223465103},
  {"branchy_glycan", 'B', 2, "-", 527.18501935651}, {"branchy_glycan", 'B', 2, "3", 162.05282341850},
  {"branchy_glycan", 'B', 2, "4", 365.13219593801}, {"branchy_glycan", 'B', 3, "-", 365.13219593801},
  {"branchy_glycan", 'B', 3, "4", 203.07937251951}, {"branchy_glycan", 'B', 4, "-", 162.05282341850},
  {"branchy_glycan", 'B', 5, "-", 892.31721529452}, {"branchy_glycan", 'B', 5, "6", 527.18501935651},
  {"branchy_glycan", 'B', 5, "6,8", 162.05282341850}, {"branchy_glycan", 'B', 5, "6,9", 365.13219593801},
  {"branchy_glycan", 'B', 5, "7", 730.26439187602}, {"branchy_glycan", 'B', 5, "7,8", 365.13219593801},
  {"branchy_glycan", 'B', 5, "7,9", 568.21156845752}, {"branchy_glycan", 'B', 5, "8", 527.18501935651},
  {"branchy_glycan", 'B', 5, "9", 730.26439187602}, {"branchy_glycan", 'B', 6, "-", 365.13219593801},
  {"branchy_glycan", 'B', 6, "7", 203.07937251951}, {"branchy_glycan", 'B', 7, "-", 162.05282341850},
  {"branchy_glycan", 'B', 8, "-", 365.13219593801}, {"branchy_glycan", 'B', 8, "9", 203.07937251951},
  {"branchy_glycan", 'B', 9, "-", 162.05282341850}, {"branchy_glycan", 'C', 1, "-", 1599.56562275323},
  {"branchy_glycan", 'C', 1, "2", 1072.38060339672}, {"branchy_glycan", 'C', 1, "2,5", 180.06338810220},
  {"branchy_glycan", 'C', 1, "2,6", 707.24840745871}, {"branchy_glycan", 'C', 1, "2,7", 910.32777997822},
  {"branchy_glycan", 'C', 1, "2,8", 707.24840745871}, {"branchy_glycan", 'C', 1, "2,9", 910.32777997822},
  {"branchy_glycan", 'C', 1, "3", 1234.43342681522}, {"branchy_glycan", 'C', 1, "3,5", 342.11621152070},
  {"branchy_glycan", 'C', 1, "3,6", 869.30123087721}, {"branchy_glycan", 'C', 1, "3,7", 1072.38060339672},
  {"branchy_glycan", 'C', 1, "3,8", 869.30123087721}, {"branchy_glycan", 'C', 1, "3,9", 1072.38060339672},
  {"branchy_glycan", 'C', 1, "4", 1437.51279933473}, {"branchy_glycan", 'C', 1, "4,5", 545.19558404021},
  {"branchy_glycan", 'C', 1, "4,6", 1072.38060339672}, {"branchy_glycan", 'C', 1, "4,7", 1275.45997591623},
  {"branchy_glycan", 'C', 1, "4,8", 1072.38060339672}, {"branchy_glycan", 'C', 1, "4,9", 1275.45997591623},
  {"branchy_glycan", 'C', 1, "5", 707.24840745871}, {"branchy_glycan", 'C', 1, "6", 1234.43342681522},
  {"branchy_glycan", 'C', 1, "6,8", 869.30123087721}, {"branchy_glycan", 'C', 1, "6,9", 1072.38060339672},
  {"branchy_glycan", 'C', 1, "7", 1437.51279933473}, {"branchy_glycan", 'C', 1, "7,8", 1072.38060339672},
  {"branchy_glycan", 'C', 1, "7,9", 1275.45997591623}, {"branchy_glycan", 'C', 1, "8", 1234.43342681522},
  {"branchy_glycan", 'C', 1, "9", 1437.51279933473}, {"branchy_glycan", 'C', 2, "-", 545.19558404021},
  {"branchy_glycan", 'C', 2, "3", 180.06338810220}, {"branchy_glycan", 'C', 2, "4", 383.14276062171},
  {"branchy_glycan", 'C', 3, "-", 383.14276062171}, {"branchy_glycan", 'C', 3, "4", 221.08993720321},
  {"branchy_glycan", 'C', 4, "-", 180.06338810220}, {"branchy_glycan", 'C', 5, "-", 910.32777997822},
  {"branchy_glycan", 'C', 5, "6", 545.19558404021}, {"branchy_glycan", 'C', 5, "6,8", 180.06338810220},
  {"branchy_glycan", 'C', 5, "6,9", 383.14276062171}, {"branchy_glycan", 'C', 5, "7", 748.27495655972},
  {"branchy_glycan", 'C', 5, "7,8", 383.14276062171}, {"branchy_glycan", 'C', 5, "7,9", 586.22213314122},
  {"branchy_glycan", 'C', 5, "8", 545.19558404021}, {"branchy_glycan", 'C', 5, "9", 748.27495655972},
  {"branchy_glycan", 'C', 6, "-", 383.14276062171}, {"branchy_glycan", 'C', 6, "7", 221.08993720321},
  {"branchy_glycan", 'C', 7, "-", 180.06338810220}, {"branchy_glycan", 'C', 8, "-", 383.14276062171},
  {"branchy_glycan", 'C', 8, "9", 221.08993720321}, {"branchy_glycan", 'C', 9, "-", 180.06338810220},
  {"branchy_glycan", 'Y', -1, "1", 221.08993720321}, {"branchy_glycan", 'Y', -1, "2", 1275.45997591623},
  {"branchy_glycan", 'Y', -1, "2,5", 383.14276062171}, {"branchy_glycan", 'Y', -1, "2,6", 910.32777997822},
  {"branchy_glycan", 'Y', -1, "2,6,8", 545.19558404021}, {"branchy_glycan", 'Y', -1, "2,6,9", 748.27495655972},
  {"branchy_glycan", 'Y', -1, "2,7", 1113.40715249773}, {"branchy_glycan", 'Y', -1, "2,7,8", 748.27495655972},
  {"branchy_glycan", 'Y', -1, "2,7,9", 951.35432907923}, {"branchy_glycan", 'Y', -1, "2,8", 910.32777997822},
  {"branchy_glycan", 'Y', -1, "2,9", 1113.40715249773}, {"branchy_glycan", 'Y', -1, "3", 1437.51279933473},
  {"branchy_glycan", 'Y', -1, "3,5", 545.19558404021}, {"branchy_glycan", 'Y', -1, "3,6", 1072.38060339672},
  {"branchy_glycan", 'Y', -1, "3,6,8", 707.24840745871}, {"branchy_glycan", 'Y', -1, "3,6,9", 910.32777997822},
  {"branchy_glycan", 'Y', -1, "3,7", 1275.45997591623}, {"branchy_glycan", 'Y', -1, "3,7,8", 910.32777997822},
  {"branchy_glycan", 'Y', -1, "3,7,9", 1113.40715249773}, {"branchy_glycan", 'Y', -1, "3,8", 1072.38060339672},
  {"branchy_glycan", 'Y', -1, "3,9", 1275.45997591623}, {"branchy_glycan", 'Y', -1, "4", 1640.59217185424},
  {"branchy_glycan", 'Y', -1, "4,5", 748.27495655972}, {"branchy_glycan", 'Y', -1, "4,6", 1275.45997591623},
  {"branchy_glycan", 'Y', -1, "4,6,8", 910.32777997822}, {"branchy_glycan", 'Y', -1, "4,6,9", 1113.40715249773},
  {"branchy_glycan", 'Y', -1, "4,7", 1478.53934843574}, {"branchy_glycan", 'Y', -1, "4,7,8", 1113.40715249773},
  {"branchy_glycan", 'Y', -1, "4,7,9", 1316.48652501724}, {"branchy_glycan", 'Y', -1, "4,8", 1275.45997591623},
  {"branchy_glycan", 'Y', -1, "4,9", 1478.53934843574}, {"branchy_glycan", 'Y', -1, "5", 910.32777997822},
  {"branchy_glycan", 'Y', -1, "6", 1437.51279933473}, {"branchy_glycan", 'Y', -1, "6,8", 1072.38060339672},
  {"branchy_glycan", 'Y', -1, "6,9", 1275.45997591623}, {"branchy_glycan", 'Y', -1, "7", 1640.59217185424},
  {"branchy_glycan", 'Y', -1, "7,8", 1275.45997591623}, {"branchy_glycan", 'Y', -1, "7,9", 1478.53934843574},
  {"branchy_glycan", 'Y', -1, "8", 1437.51279933473}, {"branchy_glycan", 'Y', -1, "9", 1640.59217185424},
  {"branchy_glycan", 'Z', -1, "1", 203.07937251951}, {"branchy_glycan", 'Z', -1, "2", 1257.44941123253},
  {"branchy_glycan", 'Z', -1, "2,5", 347.12163125431}, {"branchy_glycan", 'Z', -1, "2,6", 874.30665061082},
  {"branchy_glycan", 'Z', -1, "2,6,8", 491.16388998911}, {"branchy_glycan", 'Z', -1, "2,6,9", 694.24326250862},
  {"branchy_glycan", 'Z', -1, "2,7", 1077.38602313033}, {"branchy_glycan", 'Z', -1, "2,7,8", 694.24326250862},
  {"branchy_glycan", 'Z', -1, "2,7,9", 897.32263502813}, {"branchy_glycan", 'Z', -1, "2,8", 874.30665061082},
  {"branchy_glycan", 'Z', -1, "2,9", 1077.38602313033}, {"branchy_glycan", 'Z', -1, "3", 1419.50223465103},
  {"branchy_glycan", 'Z', -1, "3,5", 509.17445467281}, {"branchy_glycan", 'Z', -1, "3,6", 1036.35947402932},
  {"branchy_glycan", 'Z', -1, "3,6,8", 653.21671340761}, {"branchy_glycan", 'Z', -1, "3,6,9", 856.29608592712},
  {"branchy_glycan", 'Z', -1, "3,7", 1239.43884654883}, {"branchy_glycan", 'Z', -1, "3,7,8", 856.29608592712},
  {"branchy_glycan", 'Z', -1, "3,7,9", 1059.37545844663}, {"branchy_glycan", 'Z', -1, "3,8", 1036.35947402932},
  {"branchy_glycan", 'Z', -1, "3,9", 1239.43884654883}, {"branchy_glycan", 'Z', -1, "4", 1622.58160717054},
  {"branchy_glycan", 'Z', -1, "4,5", 712.25382719232}, {"branchy_glycan", 'Z', -1, "4,6", 1239.43884654883},
  {"branchy_glycan", 'Z', -1, "4,6,8", 856.29608592712}, {"branchy_glycan", 'Z', -1, "4,6,9", 1059.37545844663},
  {"branchy_glycan", 'Z', -1, "4,7", 1442.51821906834}, {"branchy_glycan", 'Z', -1, "4,7,8", 1059.37545844663},
  {"branchy_glycan", 'Z', -1, "4,7,9", 1262.45483096614}, {"branchy_glycan", 'Z', -1, "4,8", 1239.43884654883},
  {"branchy_glycan", 'Z', -1, "4,9", 1442.51821906834}, {"branchy_glycan", 'Z', -1, "5", 892.31721529452},
  {"branchy_glycan", 'Z', -1, "6", 1419.50223465103}, {"branchy_glycan", 'Z', -1, "6,8", 1036.35947402932},
  {"branchy_glycan", 'Z', -1, "6,9", 1239.43884654883}, {"branchy_glycan", 'Z', -1, "7", 1622.58160717054},
  {"branchy_glycan", 'Z', -1, "7,8", 1239.43884654883}, {"branchy_glycan", 'Z', -1, "7,9", 1442.51821906834},
  {"branchy_glycan", 'Z', -1, "8", 1419.50223465103}, {"branchy_glycan", 'Z', -1, "9", 1622.58160717054},
  {"G58143RL", 'B', 1, "-", 2081.64988252265}, {"G58143RL", 'B', 1, "2", 203.07937251951}, {"G58143RL", 'B', 1, "3", 1263.36944665967},
  {"G58143RL", 'B', 1, "3,7", 365.13219593801}, {"G58143RL", 'B', 1, "3,8", 527.18501935651}, {"G58143RL", 'B', 1, "3,9", 810.22120673470},
  {"G58143RL", 'B', 1, "3,10", 972.27403015320}, {"G58143RL", 'B', 1, "4", 1425.42227007817}, {"G58143RL", 'B', 1, "4,7", 527.18501935651},
  {"G58143RL", 'B', 1, "4,8", 689.23784277501}, {"G58143RL", 'B', 1, "4,9", 972.27403015320}, {"G58143RL", 'B', 1, "4,10", 1134.32685357170},
  {"G58143RL", 'B', 1, "5", 1628.50164259768}, {"G58143RL", 'B', 1, "5,7", 730.26439187602}, {"G58143RL", 'B', 1, "5,8", 892.31721529452},
  {"G58143RL", 'B', 1, "5,9", 1175.35340267271}, {"G58143RL", 'B', 1, "5,10", 1337.40622609121}, {"G58143RL", 'B', 1, "6", 1790.55446601618},
  {"G58143RL", 'B', 1, "6,7", 892.31721529452}, {"G58143RL", 'B', 1, "6,8", 1054.37003871302}, {"G58143RL", 'B', 1, "6,9", 1337.40622609121},
  {"G58143RL", 'B', 1, "6,10", 1499.45904950971}, {"G58143RL", 'B', 1, "7", 1183.41263180099}, {"G58143RL", 'B', 1, "8", 1345.46545521949},
  {"G58143RL", 'B', 1, "9", 1628.50164259768}, {"G58143RL", 'B', 1, "10", 1790.55446601618}, {"G58143RL", 'B', 2, "-", 1878.57051000314},
  {"G58143RL", 'B', 2, "3", 1060.29007414016}, {"G58143RL", 'B', 2, "3,7", 162.05282341850}, {"G58143RL", 'B', 2, "3,8", 324.10564683700},
  {"G58143RL", 'B', 2, "3,9", 607.14183421519}, {"G58143RL", 'B', 2, "3,10", 769.19465763369}, {"G58143RL", 'B', 2, "4", 1222.34289755866},
  {"G58143RL", 'B', 2, "4,7", 324.10564683700}, {"G58143RL", 'B', 2, "4,8", 486.15847025550}, {"G58143RL", 'B', 2, "4,9", 769.19465763369},
  {"G58143RL", 'B', 2, "4,10", 931.24748105219}, {"G58143RL", 'B', 2, "5", 1425.42227007817}, {"G58143RL", 'B', 2, "5,7", 527.18501935651},
  {"G58143RL", 'B', 2, "5,8", 689.23784277501}, {"G58143RL", 'B', 2, "5,9", 972.27403015320}, {"G58143RL", 'B', 2, "5,10", 1134.32685357170},
  {"G58143RL", 'B', 2, "6", 1587.47509349667}, {"G58143RL", 'B', 2, "6,7", 689.23784277501}, {"G58143RL", 'B', 2, "6,8", 851.29066619351},
  {"G58143RL", 'B', 2, "6,9", 1134.32685357170}, {"G58143RL", 'B', 2, "6,10", 1296.37967699020}, {"G58143RL", 'B', 2, "7", 980.33325928148},
  {"G58143RL", 'B', 2, "8", 1142.38608269998}, {"G58143RL", 'B', 2, "9", 1425.42227007817}, {"G58143RL", 'B', 2, "10", 1587.47509349667},
  {"G58143RL", 'B', 3, "-", 818.28043586298}, {"G58143RL", 'B', 3, "4", 162.05282341850}, {"G58143RL", 'B', 3, "5", 365.13219593801},
  {"G58143RL", 'B', 3, "6", 527.18501935651}, {"G58143RL", 'B', 4, "-", 656.22761244448}, {"G58143RL", 'B', 4, "5", 203.07937251951},
  {"G58143RL", 'B', 4, "6", 365.13219593801}, {"G58143RL", 'B', 5, "-", 453.14823992497}, {"G58143RL", 'B', 5, "6", 162.05282341850},
  {"G58143RL", 'B', 6, "-", 291.09541650647}, {"G58143RL", 'B', 7, "-", 898.23725072166}, {"G58143RL", 'B', 7, "8", 162.05282341850},
  {"G58143RL", 'B', 7, "9", 445.08901079669}, {"G58143RL", 'B', 7, "10", 607.14183421519}, {"G58143RL", 'B', 8, "-", 736.18442730316},
  {"G58143RL", 'B', 8, "9", 283.03618737819}, {"G58143RL", 'B', 8, "10", 445.08901079669}, {"G58143RL", 'B', 9, "-", 453.14823992497},
  {"G58143RL", 'B', 9, "10", 162.05282341850}, {"G58143RL", 'B', 10, "-", 291.09541650647}, {"G58143RL", 'B', 11, "-", 146.05790879894},
  {"G58143RL", 'C', 1, "-", 2099.66044720635}, {"G58143RL", 'C', 1, "2", 221.08993720321}, {"G58143RL", 'C', 1, "3", 1281.38001134337},
  {"G58143RL", 'C', 1, "3,7", 383.14276062171}, {"G58143RL", 'C', 1, "3,8", 545.19558404021}, {"G58143RL", 'C', 1, "3,9", 828.23177141840},
  {"G58143RL", 'C', 1, "3,10", 990.28459483690}, {"G58143RL", 'C', 1, "4", 1443.43283476187}, {"G58143RL", 'C', 1, "4,7", 545.19558404021},
  {"G58143RL", 'C', 1, "4,8", 707.24840745871}, {"G58143RL", 'C', 1, "4,9", 990.28459483690}, {"G58143RL", 'C', 1, "4,10", 1152.33741825540},
  {"G58143RL", 'C', 1, "5", 1646.51220728138}, {"G58143RL", 'C', 1, "5,7", 748.27495655972}, {"G58143RL", 'C', 1, "5,8", 910.32777997822},
  {"G58143RL", 'C', 1, "5,9", 1193.36396735641}, {"G58143RL", 'C', 1, "5,10", 1355.41679077491}, {"G58143RL", 'C', 1, "6", 1808.56503069988},
  {"G58143RL", 'C', 1, "6,7", 910.32777997822}, {"G58143RL", 'C', 1, "6,8", 1072.38060339672}, {"G58143RL", 'C', 1, "6,9", 1355.41679077491},
  {"G58143RL", 'C', 1, "6,10", 1517.46961419341}, {"G58143RL", 'C', 1, "7", 1201.42319648469}, {"G58143RL", 'C', 1, "8", 1363.47601990319},
  {"G58143RL", 'C', 1, "9", 1646.51220728138}, {"G58143RL", 'C', 1, "10", 1808.56503069988}, {"G58143RL", 'C', 2, "-", 1896.58107468684},
  {"G58143RL", 'C', 2, "3", 1078.30063882386}, {"G58143RL", 'C', 2, "3,7", 180.06338810220}, {"G58143RL", 'C', 2, "3,8", 342.11621152070},
  {"G58143RL", 'C', 2, "3,9", 625.15239889889}, {"G58143RL", 'C', 2, "3,10", 787.20522231739}, {"G58143RL", 'C', 2, "4", 1240.35346224236},
  {"G58143RL", 'C', 2, "4,7", 342.11621152070}, {"G58143RL", 'C', 2, "4,8", 504.16903493920}, {"G58143RL", 'C', 2, "4,9", 787.20522231739},
  {"G58143RL", 'C', 2, "4,10", 949.25804573589}, {"G58143RL", 'C', 2, "5", 1443.43283476187}, {"G58143RL", 'C', 2, "5,7", 545.19558404021},
  {"G58143RL", 'C', 2, "5,8", 707.24840745871}, {"G58143RL", 'C', 2, "5,9", 990.28459483690}, {"G58143RL", 'C', 2, "5,10", 1152.33741825540},
  {"G58143RL", 'C', 2, "6", 1605.48565818037}, {"G58143RL", 'C', 2, "6,7", 707.24840745871}, {"G58143RL", 'C', 2, "6,8", 869.30123087721},
  {"G58143RL", 'C', 2, "6,9", 1152.33741825540}, {"G58143RL", 'C', 2, "6,10", 1314.39024167390}, {"G58143RL", 'C', 2, "7", 998.34382396518},
  {"G58143RL", 'C', 2, "8", 1160.39664738368}, {"G58143RL", 'C', 2, "9", 1443.43283476187}, {"G58143RL", 'C', 2, "10", 1605.48565818037},
  {"G58143RL", 'C', 3, "-", 836.29100054668}, {"G58143RL", 'C', 3, "4", 180.06338810220}, {"G58143RL", 'C', 3, "5", 383.14276062171},
  {"G58143RL", 'C', 3, "6", 545.19558404021}, {"G58143RL", 'C', 4, "-", 674.23817712818}, {"G58143RL", 'C', 4, "5", 221.08993720321},
  {"G58143RL", 'C', 4, "6", 383.14276062171}, {"G58143RL", 'C', 5, "-", 471.15880460867}, {"G58143RL", 'C', 5, "6", 180.06338810220},
  {"G58143RL", 'C', 6, "-", 309.10598119017}, {"G58143RL", 'C', 7, "-", 916.24781540536}, {"G58143RL", 'C', 7, "8", 180.06338810220},
  {"G58143RL", 'C', 7, "9", 463.09957548039}, {"G58143RL", 'C', 7, "10", 625.15239889889}, {"G58143RL", 'C', 8, "-", 754.19499198686},
  {"G58143RL", 'C', 8, "9", 301.04675206189}, {"G58143RL", 'C', 8, "10", 463.09957548039}, {"G58143RL", 'C', 9, "-", 471.15880460867},
  {"G58143RL", 'C', 9, "10", 180.06338810220}, {"G58143RL", 'C', 10, "-", 309.10598119017}, {"G58143RL", 'C', 11, "-", 164.06847348264},
  {"G58143RL", 'Y', -1, "1", 367.14784600215}, {"G58143RL", 'Y', -1, "1,11", 221.08993720321}, {"G58143RL", 'Y', -1, "2", 570.22721852166},
  {"G58143RL", 'Y', -1, "2,11", 424.16930972272}, {"G58143RL", 'Y', -1, "3", 1630.51729266182}, {"G58143RL", 'Y', -1, "3,7", 732.28004194016},
  {"G58143RL", 'Y', -1, "3,7,11", 586.22213314122}, {"G58143RL", 'Y', -1, "3,8", 894.33286535866}, {"G58143RL", 'Y', -1, "3,8,11", 748.27495655972},
  {"G58143RL", 'Y', -1, "3,9", 1177.36905273685}, {"G58143RL", 'Y', -1, "3,9,11", 1031.31114393791}, {"G58143RL", 'Y', -1, "3,10", 1339.42187615535},
  {"G58143RL", 'Y', -1, "3,10,11", 1193.36396735641}, {"G58143RL", 'Y', -1, "3,11", 1484.45938386288}, {"G58143RL", 'Y', -1, "4", 1792.57011608032},
  {"G58143RL", 'Y', -1, "4,7", 894.33286535866}, {"G58143RL", 'Y', -1, "4,7,11", 748.27495655972}, {"G58143RL", 'Y', -1, "4,8", 1056.38568877716},
  {"G58143RL", 'Y', -1, "4,8,11", 910.32777997822}, {"G58143RL", 'Y', -1, "4,9", 1339.42187615535}, {"G58143RL", 'Y', -1, "4,9,11", 1193.36396735641},
  {"G58143RL", 'Y', -1, "4,10", 1501.47469957385}, {"G58143RL", 'Y', -1, "4,10,11", 1355.41679077491},
  {"G58143RL", 'Y', -1, "4,11", 1646.51220728138}, {"G58143RL", 'Y', -1, "5", 1995.64948859983}, {"G58143RL", 'Y', -1, "5,7", 1097.41223787817},
  {"G58143RL", 'Y', -1, "5,7,11", 951.35432907923}, {"G58143RL", 'Y', -1, "5,8", 1259.46506129667}, {"G58143RL", 'Y', -1, "5,8,11", 1113.40715249773},
  {"G58143RL", 'Y', -1, "5,9", 1542.50124867486}, {"G58143RL", 'Y', -1, "5,9,11", 1396.44333987592}, {"G58143RL", 'Y', -1, "5,10", 1704.55407209336},
  {"G58143RL", 'Y', -1, "5,10,11", 1558.49616329442}, {"G58143RL", 'Y', -1, "5,11", 1849.59157980089}, {"G58143RL", 'Y', -1, "6", 2157.70231201833},
  {"G58143RL", 'Y', -1, "6,7", 1259.46506129667}, {"G58143RL", 'Y', -1, "6,7,11", 1113.40715249773}, {"G58143RL", 'Y', -1, "6,8", 1421.51788471517},
  {"G58143RL", 'Y', -1, "6,8,11", 1275.45997591623}, {"G58143RL", 'Y', -1, "6,9", 1704.55407209336},
  {"G58143RL", 'Y', -1, "6,9,11", 1558.49616329442}, {"G58143RL", 'Y', -1, "6,10", 1866.60689551186},
  {"G58143RL", 'Y', -1, "6,10,11", 1720.54898671292}, {"G58143RL", 'Y', -1, "6,11", 2011.64440321939}, {"G58143RL", 'Y', -1, "7", 1550.56047780314},
  {"G58143RL", 'Y', -1, "7,11", 1404.50256900420}, {"G58143RL", 'Y', -1, "8", 1712.61330122164}, {"G58143RL", 'Y', -1, "8,11", 1566.55539242270},
  {"G58143RL", 'Y', -1, "9", 1995.64948859983}, {"G58143RL", 'Y', -1, "9,11", 1849.59157980089}, {"G58143RL", 'Y', -1, "10", 2157.70231201833},
  {"G58143RL", 'Y', -1, "10,11", 2011.64440321939}, {"G58143RL", 'Y', -1, "11", 2302.73981972586}, {"G58143RL", 'Z', -1, "1", 349.13728131845},
  {"G58143RL", 'Z', -1, "1,11", 185.06880783581}, {"G58143RL", 'Z', -1, "2", 552.21665383796}, {"G58143RL", 'Z', -1, "2,11", 388.14818035532},
  {"G58143RL", 'Z', -1, "3", 1612.50672797812}, {"G58143RL", 'Z', -1, "3,7", 696.25891257276}, {"G58143RL", 'Z', -1, "3,7,11", 532.19043909012},
  {"G58143RL", 'Z', -1, "3,8", 858.31173599126}, {"G58143RL", 'Z', -1, "3,8,11", 694.24326250862}, {"G58143RL", 'Z', -1, "3,9", 1141.34792336945},
  {"G58143RL", 'Z', -1, "3,9,11", 977.27944988681}, {"G58143RL", 'Z', -1, "3,10", 1303.40074678795},
  {"G58143RL", 'Z', -1, "3,10,11", 1139.33227330531}, {"G58143RL", 'Z', -1, "3,11", 1448.43825449548}, {"G58143RL", 'Z', -1, "4", 1774.55955139662},
  {"G58143RL", 'Z', -1, "4,7", 858.31173599126}, {"G58143RL", 'Z', -1, "4,7,11", 694.24326250862}, {"G58143RL", 'Z', -1, "4,8", 1020.36455940976},
  {"G58143RL", 'Z', -1, "4,8,11", 856.29608592712}, {"G58143RL", 'Z', -1, "4,9", 1303.40074678795}, {"G58143RL", 'Z', -1, "4,9,11", 1139.33227330531},
  {"G58143RL", 'Z', -1, "4,10", 1465.45357020645}, {"G58143RL", 'Z', -1, "4,10,11", 1301.38509672381},
  {"G58143RL", 'Z', -1, "4,11", 1610.49107791398}, {"G58143RL", 'Z', -1, "5", 1977.63892391613}, {"G58143RL", 'Z', -1, "5,7", 1061.39110851077},
  {"G58143RL", 'Z', -1, "5,7,11", 897.32263502813}, {"G58143RL", 'Z', -1, "5,8", 1223.44393192927}, {"G58143RL", 'Z', -1, "5,8,11", 1059.37545844663},
  {"G58143RL", 'Z', -1, "5,9", 1506.48011930746}, {"G58143RL", 'Z', -1, "5,9,11", 1342.41164582482}, {"G58143RL", 'Z', -1, "5,10", 1668.53294272596},
  {"G58143RL", 'Z', -1, "5,10,11", 1504.46446924332}, {"G58143RL", 'Z', -1, "5,11", 1813.57045043349}, {"G58143RL", 'Z', -1, "6", 2139.69174733463},
  {"G58143RL", 'Z', -1, "6,7", 1223.44393192927}, {"G58143RL", 'Z', -1, "6,7,11", 1059.37545844663}, {"G58143RL", 'Z', -1, "6,8", 1385.49675534777},
  {"G58143RL", 'Z', -1, "6,8,11", 1221.42828186513}, {"G58143RL", 'Z', -1, "6,9", 1668.53294272596},
  {"G58143RL", 'Z', -1, "6,9,11", 1504.46446924332}, {"G58143RL", 'Z', -1, "6,10", 1830.58576614446},
  {"G58143RL", 'Z', -1, "6,10,11", 1666.51729266182}, {"G58143RL", 'Z', -1, "6,11", 1975.62327385199}, {"G58143RL", 'Z', -1, "7", 1532.54991311944},
  {"G58143RL", 'Z', -1, "7,11", 1368.48143963680}, {"G58143RL", 'Z', -1, "8", 1694.60273653794}, {"G58143RL", 'Z', -1, "8,11", 1530.53426305530},
  {"G58143RL", 'Z', -1, "9", 1977.63892391613}, {"G58143RL", 'Z', -1, "9,11", 1813.57045043349}, {"G58143RL", 'Z', -1, "10", 2139.69174733463},
  {"G58143RL", 'Z', -1, "10,11", 1975.62327385199}, {"G58143RL", 'Z', -1, "11", 2284.72925504216},
};
// clang-format on

std::map<std::string, GlypyCase> loadGlypyCases()
{
  std::map<std::string, GlypyCase> cases;
  for (const auto& row : GLYPY_NODES)
  {
    auto& reference = cases[row.fixture];
    const std::string symbol = row.symbol;
    const auto attachment = row.parent < 0 ? std::nullopt : std::optional<Size>(row.parent);
    const auto actual_index = symbol.starts_with("Formula:")
                                ? reference.tree.addMonosaccharide(ProForma::FormulaTag {symbol.substr(8), std::nullopt}, attachment)
                                : reference.tree.addMonosaccharide(symbol, attachment);
    if (actual_index != row.index) { throw std::runtime_error("Unexpected GlyPy node order: " + symbol); }
  }
  const std::map<char, Ion> ion_types {{'B', Ion::B}, {'C', Ion::C}, {'Y', Ion::Y}, {'Z', Ion::Z}};
  for (const auto& row : GLYPY_FRAGMENTS)
  {
    auto& reference = cases[row.fixture];
    std::vector<Size> cuts;
    if (std::string branches = row.branches; branches != "-")
    {
      std::replace(branches.begin(), branches.end(), ',', ' ');
      std::istringstream branch_stream(branches);
      for (Size cut; branch_stream >> cut;)
      {
        cuts.push_back(cut);
      }
    }
    const auto root_cut = row.root < 0 ? std::nullopt : std::optional<Size>(row.root);
    const auto inserted = reference.masses.emplace(GlypyKey {ion_types.at(row.series), root_cut, cuts}, row.mass);
    if (! inserted.second) { throw std::runtime_error("Duplicate GlyPy interpretation in " + std::string(row.fixture)); }
  }
  return cases;
}

// Reference fragment masses for `common_glycan`, as computed by GlyPy
// (https://github.com/mobiusklein/glypy, revision 8d129a8c950e8635165cda9b4b1af392d4e7289e,
// test_data/fragments-example.json). Only the kinds this generator supports are listed;
// repeated values are deliberate and record the multiplicity of that fragment.
// clang-format off: one line group per fragmentation kind.
const std::vector<std::pair<std::string, double>> GLYPY_STORED_MASSES {
  {"B", 146.05790880799998}, {"B", 146.05790880799998}, {"B", 162.05282343}, {"B", 511.190104769}, {"B", 673.2429281990001},
  {"B", 1022.3802095379999}, {"B", 1184.433032968},
  {"BY", 162.05282343}, {"BY", 162.05282343}, {"BY", 349.13728133899997}, {"BY", 349.13728133899997}, {"BY", 365.132195961}, {"BY", 511.190104769},
  {"BY", 511.190104769}, {"BY", 511.190104769}, {"BY", 527.185019391}, {"BY", 673.2429281990001}, {"BY", 860.327386108}, {"BY", 876.3223007299999},
  {"BY", 876.32230073}, {"BY", 1022.3802095379999}, {"BY", 1038.37512416}, {"BY", 1038.37512416},
  {"Z", 162.05282343}, {"Z", 324.10564686}, {"Z", 673.242928199}, {"Z", 835.2957516289999}, {"Z", 1184.433032968}, {"Z", 1200.4279475899998},
  {"Z", 1200.42794759},
  {"C", 164.068473494}, {"C", 164.068473494}, {"C", 180.063388116}, {"C", 529.200669455}, {"C", 691.253492885}, {"C", 1040.3907742239999},
  {"C", 1202.4435976539999},
  {"CY", 180.063388116}, {"CY", 180.063388116}, {"CY", 367.147846025}, {"CY", 367.147846025}, {"CY", 383.142760647}, {"CY", 529.200669455},
  {"CY", 529.200669455}, {"CY", 529.200669455}, {"CY", 545.1955840769999}, {"CY", 691.253492885}, {"CY", 878.337950794}, {"CY", 894.3328654159999},
  {"CY", 894.332865416}, {"CY", 1040.3907742239999}, {"CY", 1056.385688846}, {"CY", 1056.385688846},
  {"Y", 180.063388116}, {"Y", 342.116211546}, {"Y", 691.253492885}, {"Y", 853.306316315}, {"Y", 1202.4435976539999}, {"Y", 1218.438512276},
  {"Y", 1218.438512276},
  {"BYY", 203.079372531}, {"BYY", 203.079372531}, {"BYY", 365.132195961}, {"BYY", 365.132195961}, {"BYY", 365.132195961}, {"BYY", 527.185019391},
  {"BYY", 714.2694773000001}, {"BYY", 714.2694773000001}, {"BYY", 730.264391922}, {"BYY", 876.32230073}, {"BYY", 876.32230073},
  {"BYY", 892.317215352},
  {"CYY", 221.089937217}, {"CYY", 221.089937217}, {"CYY", 383.142760647}, {"CYY", 383.142760647}, {"CYY", 383.142760647}, {"CYY", 545.1955840769999},
  {"CYY", 732.280041986}, {"CYY", 732.280041986}, {"CYY", 748.2749566079999}, {"CYY", 894.332865416}, {"CYY", 894.332865416},
  {"CYY", 910.3277800379999},
  {"ZZ", 509.17445470499996}, {"ZZ", 671.227278135}, {"ZZ", 1020.3645594739999}, {"ZZ", 1020.3645594739999}, {"ZZ", 1036.359474096},
  {"YY", 545.1955840769999}, {"YY", 707.2484075069999}, {"YY", 1056.385688846}, {"YY", 1056.385688846}, {"YY", 1072.380603468},
  {"ZZZ", 856.2960859799999},
  {"YYY", 910.3277800379999},
};
// clang-format on
} // namespace

START_TEST(TheoreticalGlycanSpectrumGenerator, "$Id$")
TOLERANCE_ABSOLUTE(0.00001)
TOLERANCE_RELATIVE(1.00000001)

START_SECTION((diagnostic ions from a ProForma composition))
{
  Generator::Options options;
  options.add_b_ions = options.add_y_ions = false;
  Generator generator(options);
  const auto parsed = ProForma::parse("N[Glycan:HexNAc2Hex5NeuAc1]ST");
  const auto& residue = std::get<ProForma::SequenceElement>(parsed.sequence[0]);
  const auto glycan = std::get<Composition>(residue.modifications[0].alternatives[0].first);
  const auto ions = generator.getFragments(glycan);
  TEST_EQUAL(ions.size(), 12)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc").getMZ(), 204.08665)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-H2O").getMZ(), 186.07608)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-H4O2").getMZ(), 168.06552)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-C2H4O2").getMZ(), 144.06552)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-CH6O3").getMZ(), 138.05495)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-C2H6O3").getMZ(), 126.05495)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:Hex1HexNAc1").getMZ(), 366.13947)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:Neu5Ac").getMZ(), 292.10269)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:Neu5Ac-H2O").getMZ(), 274.09213)
  const auto aliases = generator.getFragments(comp({{"NeuAc", 1}, {"Neu5Ac", 2}, {"NeuGc", 1}, {"Fucose", 1}}));
  TEST_EQUAL(aliases.size(), 5)
  TEST_REAL_SIMILAR(findIon(aliases, "glycan:diagnostic:Neu5Gc").getMZ(), 308.09761)
  TEST_REAL_SIMILAR(findIon(aliases, "glycan:diagnostic:Fuc").getMZ(), 147.06519)
  for (const auto& ion : ions)
  {
    TEST_EQUAL(ion.charge, 1)
  }
  TEST_TRUE(generator.getFragments(Composition {}).empty())
}
END_SECTION

START_SECTION((deoxyhexose diagnostic aliases preserve canonical and stereochemical identity))
{
  Generator::Options options;
  options.add_b_ions = options.add_y_ions = false;
  Generator generator(options);
  for (const auto& symbol : {"dHex", "d-Hex"})
  {
    const auto ions = generator.getFragments(comp({{symbol, 1}}));
    TEST_EQUAL(ions.size(), 1)
    const auto& ion = findIon(ions, "glycan:diagnostic:d-Hex");
    TEST_REAL_SIMILAR(ion.getMZ(), 147.06519)
    TEST_EQUAL(ion.composition.components.size(), 1)
    TEST_EQUAL(std::get<std::string>(ion.composition.components[0].first), "d-Hex")
    TEST_EQUAL(*MzPAF::parse(ion.getAnnotation()).named_compound, ion.name)
  }
  // Aliases describe the same component; Fuc additionally specifies stereochemistry.
  TEST_EQUAL(generator.getFragments(comp({{"dHex", 1}, {"d-Hex", 2}})).size(), 1)
  const auto mixed = generator.getFragments(comp({{"Fuc", 1}, {"dHex", 1}}));
  TEST_EQUAL(mixed.size(), 2)
  TEST_REAL_SIMILAR(findIon(mixed, "glycan:diagnostic:Fuc").getMZ(), findIon(mixed, "glycan:diagnostic:d-Hex").getMZ())
}
END_SECTION

START_SECTION((bounded composition B / Y / C / Z ions, charge ranges and formulas))
{
  Generator::Options options;
  options.add_diagnostic_ions = false;
  options.max_composition_size = 1;
  options.add_c_ions = options.add_z_ions = true;
  Generator generator(options);
  const auto ions = generator.getFragments(comp({{"Hex", 2}, {"HexNAc", 1}}));
  // Two retained compositions * four series * two charges, plus Y0 (Z0 has zero neutral mass).
  TEST_EQUAL(ions.size(), 18)
  const auto& b = findIon(ions, "glycan:B:composition:Hex1");
  const auto& c = findIon(ions, "glycan:C:composition:Hex1");
  const auto& y = findIon(ions, "glycan:Y:composition:Hex1");
  const auto& z = findIon(ions, "glycan:Z:composition:Hex1");
  TEST_REAL_SIMILAR(b.getMZ(), 163.06010)
  TEST_REAL_SIMILAR(c.getMZ() - b.getMZ(), 18.010565)
  TEST_REAL_SIMILAR(y.getMZ() - z.getMZ(), 18.010565)
  TEST_REAL_SIMILAR(findIon(ions, b.name, 2).getMZ(), (b.getMZ() + Constants::PROTON_MASS_U) / 2)
  TEST_FALSE(b.root_cleavage.has_value())
  TEST_FALSE(b.attachment_position.has_value())
  Composition formula;
  formula.components.emplace_back(ProForma::FormulaTag {"C6H10O5", std::nullopt}, 2);
  const auto custom = generator.getFragments(formula);
  TEST_REAL_SIMILAR(findIon(custom, "glycan:B:composition:Formula(C6H10O5)1").getMZ(), b.getMZ())
  const auto duplicate = generator.getFragments(comp({{"Hex", 1}, {"Hex", 1}, {"HexNAc", 1}, {"Hex", 0}}));
  TEST_EQUAL(duplicate.size(), ions.size())
  for (Size i = 0; i < ions.size(); ++i)
  {
    TEST_EQUAL(duplicate[i].name, ions[i].name)
    TEST_REAL_SIMILAR(duplicate[i].getMZ(), ions[i].getMZ())
  }
}
END_SECTION

START_SECTION((single neutral losses are feasible, independent and residue - specific))
{
  Generator::Options options;
  options.add_diagnostic_ions = options.add_y_ions = false;
  options.max_composition_size = 1;
  options.max_charge = 1;
  options.neutral_losses = {EmpiricalFormula("H2O"), EmpiricalFormula("C100"), EmpiricalFormula("H2O")};
  options.specific_neutral_losses["HexNAc"] = {EmpiricalFormula("C2H4O2")};
  const auto ions = Generator(options).getFragments(comp({{"Hex", 1}, {"HexNAc", 1}}));
  TEST_EQUAL(ions.size(), 5)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:B:composition:HexNAc1-H2O1").getMZ(), 186.07608)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:B:composition:HexNAc1-C2H4O2").getMZ(), 144.06552)
  for (const auto& ion : ions)
  {
    TEST_TRUE(ion.name.find("C100") == std::string::npos)
  }
}
END_SECTION

START_SECTION((branched tree fragments preserve cleavage identity and connectivity))
{
  GlycanStructure tree;
  tree.addMonosaccharide(std::string("HexNAc")); // 0
  tree.addMonosaccharide(std::string("Hex"), 0); // 1
  tree.addMonosaccharide(std::string("Hex"), 0); // 2
  tree.addMonosaccharide(std::string("Fuc"), 1); // 3
  Generator::Options options;
  options.add_diagnostic_ions = false;
  options.max_charge = 1;
  options.add_c_ions = options.add_z_ions = true;
  const auto ions = Generator(options).getFragments(tree);
  const auto& internal = findIon(ions, "glycan:B:tree:root=1;cuts=3;retained=Hex1");
  TEST_EQUAL(*internal.root_cleavage, 1)
  TEST_EQUAL(internal.branch_cleavages.size(), 1)
  TEST_EQUAL(internal.branch_cleavages[0], 3)
  TEST_REAL_SIMILAR(internal.getMZ(), 163.06010)
  const auto& y = findIon(ions, "glycan:Y:tree:root=0;cuts=1,2;retained=HexNAc1");
  TEST_REAL_SIMILAR(y.neutral_mass, EmpiricalFormula("C8H15NO6").getMonoWeight())
  const auto& z = findIon(ions, "glycan:Z:tree:root=0;cuts=1,2;retained=HexNAc1");
  TEST_REAL_SIMILAR(y.neutral_mass - z.neutral_mass, 2 * EmpiricalFormula("H2O").getMonoWeight())
  // Node 2 is a distinct, isobaric terminal fragment; do not collapse its interpretation.
  TEST_REAL_SIMILAR(findIon(ions, "glycan:B:tree:root=2;cuts=;retained=Hex1").getMZ(), internal.getMZ())
  options.add_internal_fragments = false;
  for (const auto& ion : Generator(options).getFragments(tree))
  {
    if (ion.ion_type == Ion::B || ion.ion_type == Ion::C) { TEST_TRUE(ion.branch_cleavages.empty()) }
  }
  options.allow_structural = false;
  const auto fallback = Generator(options).getFragments(tree);
  const auto composition = Generator(options).getFragments(tree.getComposition());
  TEST_EQUAL(fallback.size(), composition.size())
  for (Size i = 0; i < fallback.size(); ++i)
  {
    TEST_EQUAL(fallback[i].name, composition[i].name)
  }
  options.allow_structural = true;
  options.max_cleavages = 1;
  for (const auto& ion : Generator(options).getFragments(tree))
  {
    TEST_TRUE(ion.branch_cleavages.size() + (ion.root_cleavage.has_value() ? 1 : 0) <= 1)
  }
}
END_SECTION

START_SECTION((GlyPy stored fragment masses match the supported glycosidic series))
{
  // Expected values are GlyPy's; see GLYPY_STORED_MASSES above for the pinned source.
  TOLERANCE_ABSOLUTE(0.0001)
  TOLERANCE_RELATIVE(1.0)
  Generator::Options options;
  options.add_diagnostic_ions = false;
  options.add_c_ions = options.add_z_ions = true;
  options.max_cleavages = 3;
  options.max_charge = 1;
  const auto cases = loadGlypyCases();
  std::map<std::string, std::vector<double>> observed, expected;
  const std::map<Ion, char> series {{Ion::B, 'B'}, {Ion::C, 'C'}, {Ion::Y, 'Y'}, {Ion::Z, 'Z'}};
  for (const auto& fragment : Generator(options).getFragments(cases.at("common_glycan").tree))
  {
    // GlyPy has no virtual bond below the reducing-end residue.
    if (fragment.root_cleavage == 0 || fragment.branch_cleavages == std::vector<Size> {0}) { continue; }
    const auto kind = fragment.root_cleavage ? std::string(1, series.at(fragment.ion_type)) + std::string(fragment.branch_cleavages.size(), 'Y')
                                             : std::string(fragment.branch_cleavages.size(), series.at(fragment.ion_type));
    observed[kind].push_back(fragment.neutral_mass);
  }
  for (const auto& [kind, mass] : GLYPY_STORED_MASSES)
  {
    expected[kind].push_back(mass);
  }
  TEST_EQUAL(GLYPY_STORED_MASSES.size(), 96)
  TEST_EQUAL(observed.size(), expected.size())
  for (auto& [kind, masses] : expected)
  {
    STATUS("GlyPy stored kind: " << kind)
    auto& actual = observed[kind];
    TEST_EQUAL(actual.size(), masses.size())
    std::sort(actual.begin(), actual.end());
    std::sort(masses.begin(), masses.end());
    for (Size i = 0; i < std::min(actual.size(), masses.size()); ++i)
    {
      TEST_REAL_SIMILAR(actual[i], masses[i])
    }
  }
  TOLERANCE_ABSOLUTE(0.00001)
  TOLERANCE_RELATIVE(1.00000001)
}
END_SECTION

START_SECTION((GlyPy tree fixtures preserve cleavage identities and masses at charges one to three))
{
  TOLERANCE_ABSOLUTE(0.0001)
  TOLERANCE_RELATIVE(1.0)
  const auto cases = loadGlypyCases();
  TEST_EQUAL(cases.size(), 3)
  Size comparisons = 0;
  for (const auto& [fixture, reference] : cases)
  {
    for (Size max_cleavages = 1; max_cleavages <= 3; ++max_cleavages)
    {
      STATUS("GlyPy fixture: " << fixture << "; maximum cleavages: " << max_cleavages)
      Generator::Options options;
      options.add_diagnostic_ions = false;
      options.add_c_ions = options.add_z_ions = true;
      options.max_charge = 3;
      options.max_cleavages = max_cleavages;
      using ChargedKey = std::pair<GlypyKey, Int>;
      std::map<ChargedKey, Generator::Fragment> observed;
      for (const auto& fragment : Generator(options).getFragments(reference.tree))
      {
        if (fragment.root_cleavage == 0 || fragment.branch_cleavages == std::vector<Size> {0}) { continue; }
        const GlypyKey key {fragment.ion_type, fragment.root_cleavage, fragment.branch_cleavages};
        const auto inserted = observed.emplace(ChargedKey {key, fragment.charge}, fragment);
        TEST_TRUE(inserted.second)
      }
      Size expected_count = 0;
      for (const auto& [key, mass] : reference.masses)
      {
        const auto& [series, root, branches] = key;
        if (branches.size() + root.has_value() > max_cleavages) { continue; }
        for (Int charge = 1; charge <= 3; ++charge)
        {
          ++expected_count;
          const auto found = observed.find(ChargedKey {key, charge});
          TEST_TRUE(found != observed.end())
          if (found == observed.end()) { continue; }
          TEST_REAL_SIMILAR(found->second.neutral_mass, mass)
          // Independent proton mass used by Pyteomics 5.0.1 (see fixture README).
          TEST_REAL_SIMILAR(found->second.getMZ(), mass / charge + 1.00727646677)
          ++comparisons;
        }
      }
      TEST_EQUAL(observed.size(), expected_count)
    }
  }
  TEST_EQUAL(comparisons, 2892)
  TOLERANCE_ABSOLUTE(0.00001)
  TOLERANCE_RELATIVE(1.00000001)
}
END_SECTION

START_SECTION((glycopeptide attachment, HCD / ETD / EThcD and configurable retention))
{
  const auto peptide = AASequence::fromString("ANST");
  const auto glycan = comp({{"HexNAc", 2}, {"Hex", 3}, {"Fuc", 1}});
  Generator generator;
  const auto hcd = generator.getGlycopeptideFragments(peptide, glycan, 1, Generator::FragmentationMethod::HCD);
  const double stub_mass = EmpiricalFormula("C8H13NO5").getMonoWeight();
  const double full_mass = 2 * stub_mass + 3 * EmpiricalFormula("C6H10O5").getMonoWeight() + EmpiricalFormula("C6H10O4").getMonoWeight();
  const auto& stub = findIon(hcd, "peptide:b2;glycan=HexNAc1;site=N2");
  TEST_REAL_SIMILAR(stub.neutral_mass, peptide.getPrefix(2).getMonoWeight(Residue::BIon) + stub_mass)
  TEST_EQUAL(*stub.attachment_position, 1)
  TEST_EQUAL(stub.attachment_residue, "N")
  TEST_REAL_SIMILAR(findIon(hcd, "peptide:b2;glycan=0;site=N2").neutral_mass, peptide.getPrefix(2).getMonoWeight(Residue::BIon))
  TEST_FALSE(findIon(hcd, "peptide:b1").attachment_position.has_value())
  TEST_FALSE(findIon(hcd, "peptide:y2").attachment_position.has_value())
  TEST_REAL_SIMILAR(findIon(hcd, "glycan:Y:composition:0;site=N2").neutral_mass, peptide.getMonoWeight())
  const auto etd = generator.getGlycopeptideFragments(peptide, glycan, 1, Generator::FragmentationMethod::ETD);
  TEST_EQUAL(etd.size(), 12)
  const auto& intact = findIon(etd, "peptide:c2;glycan=Fuc1Hex3HexNAc2;site=N2");
  TEST_REAL_SIMILAR(intact.neutral_mass, peptide.getPrefix(2).getMonoWeight(Residue::CIon) + full_mass)
  // The radical NST suffix has neutral formula C11H18N3O7.
  TEST_REAL_SIMILAR(findIon(etd, "peptide:z3;glycan=Fuc1Hex3HexNAc2;site=N2").neutral_mass, 304.11447493 + full_mass)
  for (const auto& ion : etd)
  {
    TEST_TRUE(ion.ion_type == Ion::PEPTIDE)
  }
  const auto ethcd = generator.getGlycopeptideFragments(peptide, glycan, 1, Generator::FragmentationMethod::ETHCD);
  const Size glycan_ions = std::count_if(hcd.begin(), hcd.end(), [](const auto& ion) { return ion.ion_type != Ion::PEPTIDE; });
  TEST_EQUAL(ethcd.size(), glycan_ions + 2 * etd.size())
  TEST_REAL_SIMILAR(findIon(ethcd, intact.name).getMZ(), intact.getMZ())
  TEST_REAL_SIMILAR(findIon(ethcd, "peptide:b2;glycan=Fuc1Hex3HexNAc2;site=N2").neutral_mass,
                    peptide.getPrefix(2).getMonoWeight(Residue::BIon) + full_mass)
  auto options = generator.getOptions();
  Generator::PeptideRetention rule;
  rule.stripped = false;
  rule.stubs = {comp({{"HexNAc", 1}, {"Fuc", 1}})};
  options.peptide_retention['b'] = rule;
  generator.setOptions(options);
  const auto custom = generator.getGlycopeptideFragments(peptide, glycan, 1, Generator::FragmentationMethod::HCD);
  TEST_REAL_SIMILAR(findIon(custom, "peptide:b2;glycan=Fuc1HexNAc1;site=N2").neutral_mass,
                    peptide.getPrefix(2).getMonoWeight(Residue::BIon) + stub_mass + EmpiricalFormula("C6H10O4").getMonoWeight())
  // No available fucose: explicit retention rules must fail instead of inventing mass.
  TEST_EXCEPTION(Exception::InvalidParameter,
                 generator.getGlycopeptideFragments(peptide, comp({{"HexNAc", 1}}), 1, Generator::FragmentationMethod::HCD))
  TEST_EXCEPTION(Exception::InvalidParameter, generator.getGlycopeptideFragments(peptide, glycan, 4, Generator::FragmentationMethod::HCD))
}
END_SECTION

START_SECTION((ETD and EThcD radical backbone masses match independent elemental expectations))
{
  Generator generator;
  const auto peptide = AASequence::fromString("ANST");
  const auto glycan = comp({{"HexNAc", 1}});
  // Fixed monoisotopic masses from elemental counts, without AASequence's ion conversions:
  // z1(T): C4H7O3; z3(NST)+HexNAc: C19H31N4O12; c2(AN)+HexNAc: C15H27N5O8.
  // Radical z = internal residue sum + O - N; charge adds one proton per unit.
  for (const auto method : {Generator::FragmentationMethod::ETD, Generator::FragmentationMethod::ETHCD})
  {
    const auto ions = generator.getGlycopeptideFragments(peptide, glycan, 1, method);
    const auto& z1 = findIon(ions, "peptide:z1");
    TEST_FALSE(z1.attachment_position.has_value())
    TEST_REAL_SIMILAR(z1.neutral_mass, 103.03951908)
    TEST_REAL_SIMILAR(z1.getMZ(), 104.04679555)
    TEST_REAL_SIMILAR(findIon(ions, "peptide:z1", 2).getMZ(), 52.52703601)
    const auto& z3 = findIon(ions, "peptide:z3;glycan=HexNAc1;site=N2");
    TEST_EQUAL(*z3.attachment_position, 1)
    TEST_REAL_SIMILAR(z3.neutral_mass, 507.19384745)
    TEST_REAL_SIMILAR(z3.getMZ(), 508.20112392)
    TEST_REAL_SIMILAR(findIon(ions, z3.name, 2).getMZ(), 254.60420019)
    TEST_REAL_SIMILAR(findIon(ions, "peptide:c2;glycan=HexNAc1;site=N2").getMZ(), 406.19323932)
    auto modified = peptide;
    modified.setCTerminalModification("Amidated"); // OH -> NH2: -0.98401558 Da.
    const auto amidated = generator.getGlycopeptideFragments(modified, glycan, 1, method);
    TEST_REAL_SIMILAR(findIon(amidated, "peptide:z1").getMZ(), 103.06277997)
    TEST_REAL_SIMILAR(findIon(amidated, z3.name).getMZ(), 507.21710834)
    TEST_REAL_SIMILAR(findIon(amidated, z3.name, 2).getMZ(), 254.11219240)
    TEST_REAL_SIMILAR(findIon(amidated, "peptide:c2;glycan=HexNAc1;site=N2").getMZ(), 406.19323932)
  }
}
END_SECTION

START_SECTION((fragment annotations reject names that cannot roundtrip through mzPAF))
{
  Generator::Fragment fragment;
  for (const auto& name :
       std::vector<std::string> {"", "glycan name", "glycan\tname", "glycan\nname", "glycan[name", "glycan]name", std::string("glycan\0name", 11)})
  {
    fragment.name = name;
    TEST_EXCEPTION(Exception::InvalidParameter, fragment.getAnnotation())
  }
  fragment.name = "glycan:diagnostic:d-Hex";
  fragment.charge = 2;
  fragment.neutral_mass = 146.0579088;
  const auto parsed = MzPAF::parse(fragment.getAnnotation());
  TEST_EQUAL(*parsed.named_compound, fragment.name)
  TEST_EQUAL(*parsed.charge, 2)
  // The public C++ fragment type can also reach serialization through spectrum export.
  fragment.name = "glycan name";
  TEST_EXCEPTION(Exception::InvalidParameter, Generator::toSpectrum({fragment}))
}
END_SECTION

START_SECTION((spectrum metadata and mzPAF roundtrip))
{
  auto ions = Generator().getFragments(comp({{"HexNAc", 2}, {"Hex", 1}}));
  std::reverse(ions.begin(), ions.end());
  const auto spectrum = Generator::toSpectrum(ions);
  TEST_EQUAL(spectrum.getMSLevel(), 2)
  TEST_EQUAL(spectrum.size(), ions.size())
  TEST_TRUE(spectrum.isSorted())
  TEST_EQUAL(spectrum.getStringDataArrays()[0].getName(), Constants::UserParam::IonNames)
  TEST_EQUAL(spectrum.getIntegerDataArrays()[0].getName(), "Charges")
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    const auto annotation = MzPAF::parse(spectrum.getStringDataArrays()[0][i]);
    TEST_TRUE(annotation.ion_series == MzPAFIonSeries::NAMED)
    TEST_EQUAL(annotation.charge.value_or(1), spectrum.getIntegerDataArrays()[0][i])
    const auto& ion = findIon(ions, *annotation.named_compound, annotation.charge.value_or(1));
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), ion.getMZ())
  }
}
END_SECTION

START_SECTION((validation and bounded work))
{
  Generator generator;
  TEST_EXCEPTION(Exception::InvalidParameter, generator.getFragments(comp({{"Hex", -1}})))
  TEST_EXCEPTION(Exception::ElementNotFound, generator.getFragments(comp({{"Unknown", 1}})))
  Generator::Options options;
  options.max_charge = 0;
  TEST_EXCEPTION(Exception::InvalidParameter, generator.setOptions(options))
  TEST_EQUAL(generator.getOptions().max_charge, 2)
  options = Generator::Options {};
  options.max_fragments = 1;
  TEST_EXCEPTION(Exception::InvalidParameter, Generator(options).getFragments(comp({{"HexNAc", 1}})))
  options = Generator::Options {};
  options.max_states = 1;
  TEST_EXCEPTION(Exception::InvalidParameter, Generator(options).getFragments(comp({{"Hex", 1000000}})))
  options = Generator::Options {};
  options.add_b_ions = options.add_y_ions = false;
  options.max_oxonium_charge = 2;
  const auto ions = Generator(options).getFragments(comp({{"HexNAc", 1000000}}));
  TEST_EQUAL(ions.size(), 12)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc", 2).getMZ(), (204.08665 + Constants::PROTON_MASS_U) / 2)
}
END_SECTION

END_TEST
