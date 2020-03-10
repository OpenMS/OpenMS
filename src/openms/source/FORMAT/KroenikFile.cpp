// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/KroenikFile.h>

#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/FORMAT/TextFile.h>

namespace OpenMS
{
  KroenikFile::KroenikFile()
  {
  }

  KroenikFile::~KroenikFile()
  {
  }

  void KroenikFile::load(const String& filename, FeatureMap& feature_map)
  {
    // load input
    TextFile input(filename);

    // reset map
    FeatureMap fmap;
    feature_map = fmap;

    TextFile::ConstIterator it = input.begin();
    if (it == input.end()) return; // no data to load

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
