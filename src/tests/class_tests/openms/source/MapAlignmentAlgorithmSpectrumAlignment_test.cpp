// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

