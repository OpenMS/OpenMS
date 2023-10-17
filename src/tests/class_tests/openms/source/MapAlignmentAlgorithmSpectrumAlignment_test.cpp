// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>

using namespace std;
using namespace OpenMS;

/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmSpectrumAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmSpectrumAlignment* ptr = nullptr;
MapAlignmentAlgorithmSpectrumAlignment* nullPointer = nullptr;
START_SECTION((MapAlignmentAlgorithmSpectrumAlignment()))
	ptr = new MapAlignmentAlgorithmSpectrumAlignment();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentAlgorithmSpectrumAlignment()))
	delete ptr;
END_SECTION

START_SECTION((virtual void align(std::vector<PeakMap >&, std::vector<TransformationDescription>&)))
{
  MapAlignmentAlgorithmSpectrumAlignment ma;
  std::vector<PeakMap > maps;
	PeakMap map1;
	PeakMap map2;	
	for (UInt i = 0; i < 15; ++i)
	{
		for (UInt j = 1 ; j < 5; ++j)
		{
			PeakSpectrum spectrum;
			spectrum.setRT(i);
			spectrum.setMSLevel(j);
		
			for (float mz = 500.0; mz <= 900; mz += 100.0)
      {
				Peak1D peak;
				peak.setMZ(mz + i);
				peak.setIntensity(mz + i);
				spectrum.push_back(peak);  
      }
      map1.addSpectrum(spectrum);
		}
	}
	for (UInt i = 0; i < 15; ++i)
  {
    for (UInt j = 1 ; j < 5; ++j)
    {
      PeakSpectrum spectrum;
      spectrum.setRT(i * 1.2 + 200);
      spectrum.setMSLevel(j);

      for (float mz = 500.0; mz <= 900; mz += 100.0)
      {
        Peak1D peak;
        peak.setMZ(mz + i);
        peak.setIntensity(mz + i);
        spectrum.push_back(peak);  
      }
      map2.addSpectrum(spectrum);
    }
  }
	
	maps.push_back(map1);
	maps.push_back(map2);
	std::vector<TransformationDescription> transformations;
  ma.align(maps, transformations);
	String model_type = "interpolated";
	Param params;
	params.setValue("interpolation_type", "cspline");
	transformations[0].fitModel(model_type, params);
	transformations[1].fitModel(model_type, params);
  MapAlignmentTransformer::transformRetentionTimes(maps[0], transformations[0]);
  MapAlignmentTransformer::transformRetentionTimes(maps[1], transformations[1]);
	maps[0].updateRanges(-1);
	maps[1].updateRanges(-1);
  for (Size i = 0; i < maps[0].size(); ++i)
  {
		if (maps[0][i].getMSLevel() < 2)
		{
			TEST_REAL_SIMILAR(maps[0][i].getRT(), maps[1][i].getRT());
		}
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

