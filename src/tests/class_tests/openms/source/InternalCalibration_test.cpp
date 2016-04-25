// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Authors: Alexandra Zerck, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(InternalCalibration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

InternalCalibration* ptr = 0;
InternalCalibration* nullPointer = 0;
START_SECTION(InternalCalibration())
{
	ptr = new InternalCalibration();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~InternalCalibration())
{
	delete ptr;
}
END_SECTION

ptr = new InternalCalibration();

MSExperiment<> exp;
MzDataFile file;
file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_test.mzData"),exp);
std::vector<double> ref_masses;
ref_masses.push_back(1296.68476942);
ref_masses.push_back(2465.19833942);
Param param;
param.setValue("mz_tolerance", 100.);
param.setValue("mz_tolerance_unit", "ppm");

START_SECTION((template < typename InputPeakType > void calibrateMapSpectrumwise(MSExperiment< InputPeakType >& exp, const std::vector< double >& ref_masses)))
{
  TOLERANCE_ABSOLUTE(0.000001)

  ptr->setParameters(param);
  ptr->calibrateMapSpectrumwise(exp, ref_masses);
	
  TEST_REAL_SIMILAR(exp[0][14].getMZ(), 1296.68476942)
  TEST_REAL_SIMILAR(exp[0][77].getMZ(), 2465.19833942)

}
END_SECTION

START_SECTION((template < typename InputPeakType > void calibrateMapGlobally(MSExperiment< InputPeakType >& exp, const std::vector< double >& ref_masses, const String& trafo_file_name="")))
{
  TOLERANCE_ABSOLUTE(0.000001)
  ptr->setParameters(param);
  ptr->calibrateMapGlobally(exp, ref_masses);
  
  TEST_REAL_SIMILAR(exp[0][14].getMZ(), 1296.68476942)
  TEST_REAL_SIMILAR(exp[1][40].getMZ(), 1296.68476942)
  TEST_REAL_SIMILAR(exp[0][77].getMZ(), 2465.19833942)
  TEST_REAL_SIMILAR(exp[1][90].getMZ(), 2465.19833942)
}
END_SECTION

IdXMLFile id_file;
std::vector<ProteinIdentification> prot_ids;
std::vector<PeptideIdentification> pep_ids;
id_file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_1.idXML"), prot_ids, pep_ids);

START_SECTION((template < typename InputPeakType > void calibrateMapGlobally(MSExperiment< InputPeakType >& exp, const std::vector< PeptideIdentification >& ref_ids, const String& trafo_file_name="")))
{
  TOLERANCE_ABSOLUTE(0.000001)
  ptr->setParameters(param);
  ptr->calibrateMapGlobally(exp, pep_ids);
  
	TEST_REAL_SIMILAR(exp[0][14].getMZ(), 1296.68476942)
  TEST_REAL_SIMILAR(exp[1][40].getMZ(), 1296.68476942)
  TEST_REAL_SIMILAR(exp[0][77].getMZ(), 2465.19833942)
  TEST_REAL_SIMILAR(exp[1][90].getMZ(), 2465.19833942)
}
END_SECTION

FeatureMap f_map;
FeatureXMLFile f_file;
f_file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_annotated.featureXML"), f_map);
START_SECTION((void calibrateMapGlobally(FeatureMap& feature_map, const String& trafo_file_name="")))
{
  FeatureMap fmap2 = f_map;
  ptr->calibrateMapGlobally(fmap2);
  TEST_REAL_SIMILAR(fmap2[0].getMZ(), 687.841430243171)
  TEST_REAL_SIMILAR(fmap2[1].getMZ(), 720.005082366204)
  TEST_REAL_SIMILAR(fmap2[2].getMZ(), 927.493444113771)
  TEST_REAL_SIMILAR(fmap2[3].getMZ(), 1052.06529617992)
  TEST_REAL_SIMILAR(fmap2[4].getMZ(), 1224.59976809287)
  TEST_REAL_SIMILAR(fmap2[5].getMZ(), 998.486309862771)
}
END_SECTION
id_file.load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_2.idXML") ,prot_ids, pep_ids);
START_SECTION((void calibrateMapGlobally(FeatureMap& feature_map, const std::vector< PeptideIdentification >& ref_ids, const String& trafo_file_name="")))
{
  FeatureMap f_map2 = f_map;
  ptr->calibrateMapGlobally(f_map2, pep_ids);
  TEST_REAL_SIMILAR(f_map2[0].getMZ(),687.841430243171)
  TEST_REAL_SIMILAR(f_map2[1].getMZ(),720.005082366204)
  TEST_REAL_SIMILAR(f_map2[2].getMZ(),927.493444113771)
  TEST_REAL_SIMILAR(f_map2[3].getMZ(),1052.06529617992)
  TEST_REAL_SIMILAR(f_map2[4].getMZ(),1224.59976809287)
  TEST_REAL_SIMILAR(f_map2[5].getMZ(),998.486309862771)
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


