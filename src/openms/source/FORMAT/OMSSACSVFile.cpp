// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSSACSVFile.h>

#include <OpenMS/SYSTEM/File.h>

#include <fstream>

using namespace std;

namespace OpenMS
{

  OMSSACSVFile::OMSSACSVFile() = default;

  OMSSACSVFile::~OMSSACSVFile() = default;

  void OMSSACSVFile::load(const std::string & filename, ProteinIdentification & /* protein_identification */, PeptideIdentificationList & id_data) const
  {
    ifstream is(filename.c_str());
    if (!is)
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    std::string line;
    getline(is, line, '\n');
    if (!StringUtils::hasPrefix(line, "Spectrum"))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "first line does not contain a description", "");
    }

    // line number counter
    Size line_number = 0;

    // ignore first line
    Int actual_spectrum_number(-1);
    while (getline(is, line, '\n'))
    {
      ++line_number;

      // Spectrum number, Filename/id, Peptide, E-value, Mass, gi, Accession, Start, Stop, Defline, Mods, Charge, Theo Mass, P-value
      // 1,,MSHHWGYGK,0.00336754,1101.49,0,6599,1,9,CRHU2 carbonate dehydratase (EC 4.2.1.1) II [validated] - human,,1,1101.48,1.30819e-08
      StringUtils::trim(line);

      // replace ',' in protein name
      std::string::const_iterator it = find(line.begin(), line.end(), '"');
      UInt offset(0);
      if (it != line.end())
      {
        while (*(++it) != '"')
        {
          if (*it == ',')
          {
            offset++;
          }
        }
      }
      vector<std::string> split;
      StringUtils::split(line, ', ', split);
      if (split.size() != 14 && split.size() != 14 + offset)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "number of columns should be 14 in line " + StringUtils::toStr(line_number));
      }
      PeptideHit p;
      p.setSequence(AASequence::fromString(StringUtils::trim(split[2])));
      p.setScore(StringUtils::toDouble(StringUtils::trimmed(split[13 + offset])));
      p.setCharge(StringUtils::toInt32(StringUtils::trimmed(split[11 + offset])));

      if (actual_spectrum_number != StringUtils::toInt32(StringUtils::trimmed(split[0])))
      {
        // new id
        //id_data.push_back(IdentificationData());
        id_data.emplace_back();
        id_data.back().setScoreType("OMSSA");
        actual_spectrum_number = (UInt)StringUtils::toInt32(StringUtils::trimmed(split[0]));
      }

      id_data.back().insertHit(p);
    }

  }

} // namespace OpenMS
