// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/KroenikFile.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/FORMAT/TextFile.h>

namespace OpenMS
{
  KroenikFile::KroenikFile() = default;

  KroenikFile::~KroenikFile() = default;

  void KroenikFile::load(const std::string& filename, FeatureMap& feature_map)
  {
    // load input
    TextFile input(filename);

    // reset map
    FeatureMap fmap;
    feature_map = fmap;

    TextFile::ConstIterator it = input.begin();
    if (it == input.end())
    {
      return; // no data to load
    }
    // skip header line
    ++it;
    // process content
    for (; it != input.end(); ++it)
    {
      const std::string& line = *it;

      //split lines: File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications
      std::vector<std::string> parts;
      StringUtils::split(line, '\t', parts);

      if (parts.size() != 14)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                    std::string("Failed parsing in line ")
                                    + StringUtils::toStr((it - input.begin()) + 1)
                                    + ": missing 14 tab-separated entries (got "
                                    + StringUtils::toStr(parts.size())
                                    + ")\nLine was: '"
                                    + line
                                    + "'");
      }
      //create feature
      Feature f;
      f.setCharge(StringUtils::toInt32(parts[4]));
      f.setMZ(StringUtils::toDouble(parts[5]) / f.getCharge() + Constants::PROTON_MASS_U);
      f.setRT(StringUtils::toDouble(parts[11]));
      f.setOverallQuality(StringUtils::toDouble(parts[12]));
      f.setIntensity(StringUtils::toDouble(parts[8]));
      ConvexHull2D hull;
      ConvexHull2D::PointType point;

      point.setX(StringUtils::toDouble(parts[9]));
      point.setY(f.getMZ());
      hull.addPoint(point);

      point.setX(StringUtils::toDouble(parts[9]));
      point.setY(f.getMZ() + 3.0 / static_cast<double>(f.getCharge()));
      hull.addPoint(point);

      point.setX(StringUtils::toDouble(parts[10]));
      point.setY(f.getMZ() + 3.0 / static_cast<double>(f.getCharge()));
      hull.addPoint(point);

      point.setX(StringUtils::toDouble(parts[10]));
      point.setY(f.getMZ());
      hull.addPoint(point);

      point.setX(StringUtils::toDouble(parts[9]));
      point.setY(f.getMZ());
      hull.addPoint(point);

      std::vector<ConvexHull2D> hulls;
      hulls.push_back(hull);
      f.setConvexHulls(hulls);
      f.setMetaValue("Mass", StringUtils::toDouble(parts[5]));
      f.setMetaValue("FirstScan", StringUtils::toDouble(parts[1]));
      f.setMetaValue("LastScan", StringUtils::toInt32(parts[2]));
      f.setMetaValue("NumOfScans", StringUtils::toDouble(parts[3]));
      f.setMetaValue("AveragineModifications", parts[13]);
      feature_map.push_back(f);
    }

    OPENMS_LOG_INFO << "Hint: The convex hulls are approximated in m/z dimension (Kroenik lacks this information)!\n";
  }

} // namespace OpenMS
