// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/MATH/MISC/LinearResampling.h>
#include <initializer_list>
#include <sstream>
#include <utility>

using namespace OpenMS;

namespace
{
MSChromatogram makeChromatogram(std::initializer_list<std::pair<double, float>> points)
{
  MSChromatogram chrom;
  for (const auto& [rt, intensity] : points)
  {
    ChromatogramPeak peak;
    peak.setRT(rt);
    peak.setIntensity(intensity);
    chrom.push_back(peak);
  }
  return chrom;
}
} // namespace

START_TEST(LinearResampling, "$Id$")

START_SECTION([EXTRA] empty input and a single data point preserve metadata and intensity)
{
  Internal::LinearResampling resampler(0.5);
  MSChromatogram chrom;
  chrom.setNativeID("TIC");
  resampler.raster(chrom);
  TEST_TRUE(chrom.empty())
  TEST_EQUAL(chrom.getNativeID(), "TIC")

  chrom = makeChromatogram({{10.0, 7.0f}});
  chrom.setNativeID("TIC");
  resampler.raster(chrom);
  TEST_EQUAL(chrom.size(), 1)
  ABORT_IF(chrom.size() != 1)
  TEST_EQUAL(chrom[0].getRT(), 10.0)
  TEST_EQUAL(chrom[0].getIntensity(), 7.0f)
  TEST_EQUAL(chrom.getNativeID(), "TIC")
}
END_SECTION

START_SECTION([EXTRA] irregular samples use the original grid origin and conserve intensity)
{
  auto chrom = makeChromatogram({{5.0, 3.0f}, {5.5, 6.0f}, {6.0, 8.0f}, {6.6, 2.0f}, {6.8, 1.0f}});
  Internal::LinearResampling(0.75).raster(chrom);
  TEST_EQUAL(chrom.size(), 4)
  ABORT_IF(chrom.size() != 4)
  TEST_EQUAL(chrom[0].getRT(), 5.0)
  TEST_EQUAL(chrom[1].getRT(), 5.75)
  TEST_EQUAL(chrom[2].getRT(), 6.5)
  TEST_EQUAL(chrom[3].getRT(), 7.25)
  TEST_REAL_SIMILAR(chrom[0].getIntensity(), 5.0)
  TEST_REAL_SIMILAR(chrom[1].getIntensity(), 4.0 + 16.0 / 3.0)
  TEST_REAL_SIMILAR(chrom[2].getIntensity(), 5.0)
  TEST_REAL_SIMILAR(chrom[3].getIntensity(), 2.0 / 3.0)
}
END_SECTION

START_SECTION([EXTRA] iterator resampling retains existing intensities and accumulates outside the grid)
{
  const auto raw = makeChromatogram({{-1.0, 2.0f}, {0.5, 4.0f}, {2.0, 6.0f}});
  auto output = makeChromatogram({{0.0, 10.0f}, {1.0, 20.0f}});
  Internal::LinearResampling resampler(1.0);
  resampler.raster(raw.begin(), raw.end(), output.begin(), output.end());
  TEST_EQUAL(output[0].getIntensity(), 14.0f)
  TEST_EQUAL(output[1].getIntensity(), 28.0f)

  output = makeChromatogram({{0.0, 1.0f}});
  resampler.raster(raw.begin(), raw.end(), output.begin(), output.end());
  TEST_EQUAL(output[0].getIntensity(), 13.0f)
}
END_SECTION

START_SECTION([EXTRA] ppm grid retains its relative spacing and exclusive upper endpoint)
{
  std::vector<ChromatogramPeak> grid;
  Internal::LinearResampling(10000.0, true).populateRaster(grid, 100.0, 104.0, 99);
  TEST_EQUAL(grid.size(), 4)
  ABORT_IF(grid.size() != 4)
  TEST_REAL_SIMILAR(grid[0].getRT(), 100.0)
  TEST_REAL_SIMILAR(grid[1].getRT(), 101.0)
  TEST_REAL_SIMILAR(grid[2].getRT(), 102.01)
  TEST_REAL_SIMILAR(grid[3].getRT(), 103.0301)
  for (const auto& peak : grid)
  {
    TEST_EQUAL(peak.getIntensity(), 0.0f)
  }
}
END_SECTION

START_SECTION([EXTRA] spacing warnings respect the shared suppression flag and ppm units)
{
  const auto raw = makeChromatogram({{100.0, 1.0f}, {102.0, 2.0f}});
  auto access = [](auto it) { return it->getPos(); };
  auto& logger = getThreadLocalLogWarn();
  logger.rdbuf()->clearCache();
  std::ostringstream messages;
  logger.insert(messages);
  const bool previous = suppress_resampling_spacing_warning.exchange(false);

  Internal::LinearResampling(1.0).verifySpacing(raw.begin(), raw.end(), access);
  logger.flushIncomplete();
  logger.rdbuf()->clearCache();
  TEST_TRUE(messages.str().find("Resampling spacing (1)") != std::string::npos)

  messages.str("");
  suppress_resampling_spacing_warning.store(true);
  Internal::LinearResampling(1.0).verifySpacing(raw.begin(), raw.end(), access);
  logger.flushIncomplete();
  TEST_TRUE(messages.str().empty())

  suppress_resampling_spacing_warning.store(false);
  Internal::LinearResampling(1.0, true).verifySpacing(raw.begin(), raw.end(), access);
  logger.flushIncomplete();
  TEST_TRUE(messages.str().empty())

  suppress_resampling_spacing_warning.store(previous);
  logger.remove(messages);
}
END_SECTION

END_TEST
