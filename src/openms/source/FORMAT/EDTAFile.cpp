// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/EDTAFile.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace std;

namespace OpenMS
{

  EDTAFile::EDTAFile() = default;

  EDTAFile::~EDTAFile() = default;

  double EDTAFile::checkedToDouble_(const std::vector<String>& parts, Size index, double def)
  {
    if (index < parts.size() && parts[index] != "NA")
    {
      return parts[index].toDouble();
    }
    return def;
  }

  Int EDTAFile::checkedToInt_(const std::vector<String>& parts, Size index, Int def)
  {
    if (index < parts.size() && parts[index] != "NA")
    {
      return parts[index].toInt();
    }
    return def;
  }

  void EDTAFile::load(const String& filename, ConsensusMap& consensus_map)
  {
    // load input
    TextFile input(filename);
    TextFile::ConstIterator input_it = input.begin();

    // reset map
    consensus_map = ConsensusMap();
    consensus_map.setUniqueId();

    char separator = ' ';
    if (input_it->hasSubstring("\t"))
    {
      separator = '\t';
    }
    else if (input_it->hasSubstring(" "))
    {
      separator = ' ';
    }
    else if (input_it->hasSubstring(","))
    {
      separator = ',';
    }
    // parsing header line
    std::vector<String> headers;
    input_it->split(separator, headers);
    int offset = 0;
    for (Size i = 0; i < headers.size(); ++i)
    {
      headers[i].trim();
    }
    String header_trimmed = *input.begin();
    header_trimmed.trim();

    enum
    {
      TYPE_UNDEFINED,
      TYPE_OLD_NOCHARGE,
      TYPE_OLD_CHARGE,
      TYPE_CONSENSUS
    }
    input_type = TYPE_UNDEFINED;
    Size input_features = 1;

    double rt = 0.0;
    double mz = 0.0;
    double it = 0.0;
    Int ch = 0;

    if (headers.size() <= 2)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Failed parsing in line 1: not enough columns! Expected at least 3 columns!\nOffending line: '") + header_trimmed + "'  (line 1)\n");
    }
    else if (headers.size() == 3)
    {
      input_type = TYPE_OLD_NOCHARGE;
    }
    else if (headers.size() == 4)
    {
      input_type = TYPE_OLD_CHARGE;
    }
    // see if we have a header
    try
    {
      // try to convert... if not: that's a header
      rt = headers[0].toDouble();
      mz = headers[1].toDouble();
      it = headers[2].toDouble();
    }
    catch (Exception::BaseException&)
    {
      offset = 1;
      ++input_it;
      OPENMS_LOG_INFO << "Detected a header line.\n";
    }

    if (headers.size() >= 5)
    {
      if (String(headers[4].trim()).toUpper() == "RT1")
      {
        input_type = TYPE_CONSENSUS;
      }
      else
      {
        input_type = TYPE_OLD_CHARGE;
      }
    }
    if (input_type == TYPE_CONSENSUS)
    {
      // Every consensus style line includes features with four columns.
      // The remainder is meta data
      input_features = headers.size() / 4;
    }

    if (offset == 0 && (input_type == TYPE_OLD_CHARGE || input_type == TYPE_CONSENSUS))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Failed parsing in line 1: No HEADER provided. This is only allowed for three columns. You have more!\nOffending line: '") + header_trimmed + "'  (line 1)\n");
    }

    SignedSize input_size = input.end() - input.begin();

    ConsensusMap::ColumnHeader desc;
    desc.filename = filename;
    desc.size = (input_size) - offset;
    consensus_map.getColumnHeaders()[0] = desc;

    // parsing features
    consensus_map.reserve(input_size);

    for (; input_it != input.end(); ++input_it)
    {
      //do nothing for empty lines
      String line_trimmed = *input_it;
      line_trimmed.trim();
      if (line_trimmed.empty())
      {
        if ((input_it - input.begin()) < input_size - 1) OPENMS_LOG_WARN << "Notice: Empty line ignored (line " << ((input_it - input.begin()) + 1) << ").";
        continue;
      }

      //split line to tokens
      std::vector<String> parts;
      input_it->split(separator, parts);

      //abort if line does not contain enough fields
      if (parts.size() < 3)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                    String("Failed parsing in line ")
                                    + String((input_it - input.begin()) + 1)
                                    + ": At least three columns are needed! (got  "
                                    + String(parts.size())
                                    + ")\nOffending line: '"
                                    + line_trimmed
                                    + "'  (line "
                                    + String((input_it - input.begin()) + 1)
                                    + ")\n");
      }

      ConsensusFeature cf;
      cf.setUniqueId();

      try
      {
        // Convert values. Will return -1 if not available.
        rt = checkedToDouble_(parts, 0);
        mz = checkedToDouble_(parts, 1);
        it = checkedToDouble_(parts, 2);
        ch = checkedToInt_(parts, 3);

        cf.setRT(rt);
        cf.setMZ(mz);
        cf.setIntensity(it);
        if (input_type != TYPE_OLD_NOCHARGE)
          cf.setCharge(ch);
      }
      catch (Exception::BaseException&)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Failed parsing in line ") + String((input_it - input.begin()) + 1) + ": Could not convert the first three columns to a number!\nOffending line: '" + line_trimmed + "'  (line " + String((input_it - input.begin()) + 1) + ")\n");
      }

      // Check all features in one line
      for (Size j = 1; j < input_features; ++j)
      {
        try
        {
          Feature f;
          f.setUniqueId();

          // Convert values. Will return -1 if not available.
          rt = checkedToDouble_(parts, j * 4 + 0);
          mz = checkedToDouble_(parts, j * 4 + 1);
          it = checkedToDouble_(parts, j * 4 + 2);
          ch = checkedToInt_(parts, j * 4 + 3);

          // Only accept features with at least RT and MZ set
          if (rt != -1 && mz != -1)
          {
            f.setRT(rt);
            f.setMZ(mz);
            f.setIntensity(it);
            f.setCharge(ch);

            cf.insert(j - 1, f);
          }
        }
        catch (Exception::BaseException&)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", String("Failed parsing in line ") + String((input_it - input.begin()) + 1) + ": Could not convert one of the four sub-feature columns (starting at column " + (j * 4 + 1) + ") to a number! Is the correct separator specified?\nOffending line: '" + line_trimmed + "'  (line " + String((input_it - input.begin()) + 1) + ")\n");
        }
      }

      //parse meta data
      for (Size j = input_features * 4; j < parts.size(); ++j)
      {
        String part_trimmed = parts[j];
        part_trimmed.trim();
        if (!part_trimmed.empty())
        {
          //check if column name is ok
          if (headers.size() <= j || headers[j].empty())
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                        String("Error: Missing meta data header for column ") + (j + 1) + "!"
                                        + String("Offending header line: '") + header_trimmed + "'  (line 1)");
          }
          //add meta value
          cf.setMetaValue(headers[j], part_trimmed);
        }
      }

      //insert feature to map
      consensus_map.push_back(cf);
    }

    // register ColumnHeaders
    ConsensusMap::ColumnHeader fd;
    fd.filename = filename;
    fd.size = consensus_map.size();
    Size maps = std::max(input_features - 1, Size(1)); // its either a simple feature or a consensus map
    // (in this case the 'input_features' includes the centroid, which we do not count)
    for (Size i = 0; i < maps; ++i)
    {
      fd.label = String("EDTA_Map ") + String(i);
      consensus_map.getColumnHeaders()[i] = fd;
    }

  }

  void EDTAFile::store(const String& filename, const ConsensusMap& map) const
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::EDTA))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::EDTA) + "'");
    }

    TextFile tf;

    // search for maximum number of sub-features (since this determines the number of columns)
    Size max_sub(0);
    for (Size i = 0; i < map.size(); ++i)
    {
      max_sub = std::max(max_sub, map[i].getFeatures().size());
    }

    // write header
    String header("RT\tm/z\tintensity\tcharge");
    for (Size i = 1; i <= max_sub; ++i)
    {
      header += "\tRT" + String(i) + "\tm/z" + String(i) + "\tintensity" + String(i) + "\tcharge" + String(i);
    }
    tf.addLine(header);

    for (Size i = 0; i < map.size(); ++i)
    {
      ConsensusFeature f = map[i];
      // consensus
      String entry = String(f.getRT()) + "\t" + f.getMZ() + "\t" + f.getIntensity() + "\t" + f.getCharge();
      // sub-features
      const ConsensusFeature::HandleSetType& handle = f.getFeatures();
      for (ConsensusFeature::HandleSetType::const_iterator it = handle.begin(); it != handle.end(); ++it)
      {
        entry += String("\t") + it->getRT() + "\t" + it->getMZ() + "\t" + it->getIntensity() + "\t" + it->getCharge();
      }
      // missing sub-features
      for (Size j = handle.size(); j < max_sub; ++j)
      {
        entry += "\tNA\tNA\tNA\tNA";
      }
      tf.addLine(entry);
    }


    tf.store(filename);
  }

  void EDTAFile::store(const String& filename, const FeatureMap& map) const
  {
    TextFile tf;
    tf.addLine("RT\tm/z\tintensity\tcharge");

    for (Size i = 0; i < map.size(); ++i)
    {
      const Feature& f = map[i];
      tf.addLine(String(f.getRT()) + "\t" + f.getMZ() + "\t" + f.getIntensity() + "\t" + f.getCharge());
    }

    tf.store(filename);

  }

} // namespace OpenMS
