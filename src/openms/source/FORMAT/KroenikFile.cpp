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

  void KroenikFile::load(const String& filename, FeatureMap& feature_map)
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
      String line = *it;

      //split lines: File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications
      std::vector<String> parts;
      line.split('\t', parts);

      if (parts.size() != 14)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "",
                                    String("Failed parsing in line ")
                                    + String((it - input.begin()) + 1)
                                    + ": missing 14 tab-separated entries (got "
                                    + String(parts.size())
                                    + ")\nLine was: '"
                                    + line
                                    + "'");
      }
      //create feature
      Feature f;
      f.setCharge(parts[4].toInt());
      f.setMZ(parts[5].toDouble() / f.getCharge() + Constants::PROTON_MASS_U);
      f.setRT(parts[11].toDouble());
      f.setOverallQuality(parts[12].toDouble());
      f.setIntensity(parts[8].toDouble());
      ConvexHull2D hull;
      ConvexHull2D::PointType point;

      point.setX(parts[9].toDouble());
      point.setY(f.getMZ());
      hull.addPoint(point);

      point.setX(parts[9].toDouble());
      point.setY(f.getMZ() + 3.0 / static_cast<double>(f.getCharge()));
      hull.addPoint(point);

      point.setX(parts[10].toDouble());
      point.setY(f.getMZ() + 3.0 / static_cast<double>(f.getCharge()));
      hull.addPoint(point);

      point.setX(parts[10].toDouble());
      point.setY(f.getMZ());
      hull.addPoint(point);

      point.setX(parts[9].toDouble());
      point.setY(f.getMZ());
      hull.addPoint(point);

      std::vector<ConvexHull2D> hulls;
      hulls.push_back(hull);
      f.setConvexHulls(hulls);
      f.setMetaValue("Mass", parts[5].toDouble());
      f.setMetaValue("FirstScan", parts[1].toDouble());
      f.setMetaValue("LastScan", parts[2].toInt());
      f.setMetaValue("NumOfScans", parts[3].toDouble());
      f.setMetaValue("AveragineModifications", parts[13]);
      feature_map.push_back(f);
    }

    OPENMS_LOG_INFO << "Hint: The convex hulls are approximated in m/z dimension (Kroenik lacks this information)!\n";
  }

} // namespace OpenMS
