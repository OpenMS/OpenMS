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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

using namespace OpenMS;
using namespace std;

START_TEST(ChromatogramExtractorAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromatogramExtractorAlgorithm* ptr = nullptr;
ChromatogramExtractorAlgorithm* nullPointer = nullptr;

START_SECTION(ChromatogramExtractorAlgorithm())
{
	ptr = new ChromatogramExtractorAlgorithm();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ChromatogramExtractorAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION(void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, std::vector< OpenSwath::ChromatogramPtr > &output, std::vector< ExtractionCoordinates > extraction_coordinates, double mz_extraction_window, bool ppm, String filter))
{
  double extract_window = 0.05;
  boost::shared_ptr<PeakMap > exp(new PeakMap);
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.mzML"), *exp);
  OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

  ///////////////////////////////////////////////////////////////////////////
  ChromatogramExtractorAlgorithm extractor;

  std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > coordinates;
  std::vector< OpenSwath::ChromatogramPtr > out_exp;
  for (int i = 0; i < 3; i++)
  {
    OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
    out_exp.push_back(s);
  }

  {
    ChromatogramExtractorAlgorithm::ExtractionCoordinates coord;
    coord.mz = 618.31; coord.rt_start = 0; coord.rt_end = -1; coord.id = "tr1";
    coordinates.push_back(coord);
    coord.mz = 628.45; coord.rt_start = 0; coord.rt_end = -1; coord.id = "tr2";
    coordinates.push_back(coord);
    coord.mz = 654.38; coord.rt_start = 0; coord.rt_end = -1; coord.id = "tr3";
    coordinates.push_back(coord);
  }
  extractor.extractChromatograms(expptr, out_exp, coordinates, extract_window, false, "tophat");

  OpenSwath::ChromatogramPtr chrom = out_exp[0];

  TEST_EQUAL(chrom->getTimeArray()->data.size(), 59);
  TEST_EQUAL(chrom->getIntensityArray()->data.size(), 59);
  // we sort/reorder 
  int firstchromat  = 1;
  int secondchromat = 2;
  int thirdchromat  = 0;

  double max_value = -1; double foundat = -1;
  chrom = out_exp[firstchromat];
  for(Size i = 0; i < chrom->getTimeArray()->data.size(); i++)
  {
    double rt = chrom->getTimeArray()->data[i];
    double in = chrom->getIntensityArray()->data[i];
    if(in > max_value)
    {
      max_value = in;
      foundat = rt;
    }
  }
  TEST_REAL_SIMILAR(max_value, 169.792);
  TEST_REAL_SIMILAR(foundat, 3120.26);

  max_value = -1; foundat = -1;
  chrom = out_exp[secondchromat];
  for(Size i = 0; i < chrom->getTimeArray()->data.size(); i++)
  {
    double rt = chrom->getTimeArray()->data[i];
    double in = chrom->getIntensityArray()->data[i];
    if(in > max_value)
    {
      max_value = in;
      foundat = rt;
    }
  }

  TEST_REAL_SIMILAR(max_value, 577.33);
  TEST_REAL_SIMILAR(foundat, 3120.26);

  max_value = -1; foundat = -1;
  chrom = out_exp[thirdchromat];
  for(Size i = 0; i < chrom->getTimeArray()->data.size(); i++)
  {
    double rt = chrom->getTimeArray()->data[i];
    double in = chrom->getIntensityArray()->data[i];
    if(in > max_value)
    {
      max_value = in;
      foundat = rt;
    }
  }

  TEST_REAL_SIMILAR(max_value, 35.593);
  TEST_REAL_SIMILAR(foundat, 3055.16);

}
END_SECTION

///////////////////////////////////////////////////////////////////////////
/// Private functions
///////////////////////////////////////////////////////////////////////////

//  mz_a = [400+0.01*i for i in range(20)]
//  int_a = [0 + i*100.0 for i in range(10)] + [900 - i*100.0 for i in range(10)]
static const double mz_arr[] = {
  400.0 ,
  400.01,
  400.02,
  400.03,
  400.04,
  400.05,
  400.06,
  400.07,
  400.08,
  400.09,
  400.1 ,
  400.11,
  400.12,
  400.13,
  400.14,
  400.15,
  400.16,
  400.17,
  400.18,
  400.19,
  450.0,
  500.0,
};
static const double int_arr[] = {
  8.0  , 
  100.0,
  200.0,
  300.0,
  400.0,
  500.0,
  600.0,
  700.0,
  800.0,
  900.0,
  900.0,
  800.0,
  700.0,
  600.0,
  500.0,
  400.0,
  300.0,
  200.0,
  100.0,
  0.0, 
  10.0, 
  10.0, 
};

START_SECTION(void extract_value_tophat(const std::vector< double >::const_iterator &mz_start, std::vector< double >::const_iterator &mz_it, const std::vector< double >::const_iterator &mz_end, std::vector< double >::const_iterator &int_it, const double &mz, double &integrated_intensity, const double &mz_extraction_window, bool ppm))
{ 
  std::vector<double> mz (mz_arr, mz_arr + sizeof(mz_arr) / sizeof(mz_arr[0]) );
  std::vector<double> intensities (int_arr, int_arr + sizeof(int_arr) / sizeof(int_arr[0]) );

  // convert the data into a spectrum
  MSSpectrum spectrum;
  for(Size i=0; i<mz.size(); ++i)
  {
    Peak1D peak;
    peak.setMZ(mz[i]);
    peak.setIntensity(intensities[i]);
    spectrum.push_back(peak);
  }

  std::vector<double>::const_iterator mz_start = mz.begin();
  std::vector<double>::const_iterator mz_it_end = mz.end();
  std::vector<double>::const_iterator mz_it = mz.begin();
  std::vector<double>::const_iterator int_it = intensities.begin();

  double integrated_intensity = 0;
  double extract_window = 0.2; // +/- 0.1

  // If we use monotonically increasing m/z values then everything should work fine
  ChromatogramExtractorAlgorithm extractor;

  // test the zero first value
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 399.805, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 0.0); // test very first data point

  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 399.91, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 108.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.0, integrated_intensity, extract_window, false);
  // print(sum([0 + i*100.0 for i in range(10)] + 8) )
  TEST_REAL_SIMILAR( integrated_intensity, 4508.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.05,  integrated_intensity, extract_window, false);
  //print(sum([0 + i*100.0 for i in range(10)]) + sum([900 - i*100.0 for i in range(6)])  )
  TEST_REAL_SIMILAR( integrated_intensity, 8400.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.1, integrated_intensity, extract_window, false);
  //print(sum([0 + i*100.0 for i in range(10)]) + sum([900 - i*100.0 for i in range(10)])  )
  TEST_REAL_SIMILAR( integrated_intensity, 9000.0);
  TEST_EQUAL((int)integrated_intensity, 9000);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.28, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity,100.0);

  // test the very last value
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 500.0, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity, 10.0);

  // this is to document the situation of using m/z values that are not monotonically increasing:
  //  --> it might not give the correct result (9000) if we try to extract 400.1 AFTER 500.0 
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.1, integrated_intensity, extract_window, false);
  TEST_NOT_EQUAL((int)integrated_intensity,9000);

  /// use ppm extraction windows
  //

  mz_it = mz.begin();
  int_it = intensities.begin();
  integrated_intensity = 0;
  extract_window = 500; // 500 ppm == 0.2 Da @ 400 m/z

  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 399.89, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity, 0.0);  // below 400, 500ppm is below 0.2 Da...
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 399.91, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity, 8.0);  // very first value
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 399.92, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,108.0); 
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.0, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,4508.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.05, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,8400.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.1, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,9000.0);

}
END_SECTION

START_SECTION([EXTRA]void extract_value_tophat(const std::vector< double >::const_iterator &mz_start, std::vector< double >::const_iterator &mz_it, const std::vector< double >::const_iterator &mz_end, std::vector< double >::const_iterator &int_it, const double &mz, double &integrated_intensity, const double &mz_extraction_window, bool ppm))
{ 
  std::vector<double> mz (mz_arr, mz_arr + sizeof(mz_arr) / sizeof(mz_arr[0]) );
  std::vector<double> intensities (int_arr, int_arr + sizeof(int_arr) / sizeof(int_arr[0]) );

  // convert the data into a spectrum
  MSSpectrum spectrum;
  for(Size i=0; i<mz.size(); ++i)
  {
    Peak1D peak;
    peak.setMZ(mz[i]);
    peak.setIntensity(intensities[i]);
    spectrum.push_back(peak);
  }

  std::vector<double>::const_iterator mz_start = mz.begin();
  std::vector<double>::const_iterator mz_it_end = mz.end();
  std::vector<double>::const_iterator mz_it = mz.begin();
  std::vector<double>::const_iterator int_it = intensities.begin();

  double integrated_intensity = 0;
  double extract_window = 0.2; // +/- 0.1

  // If we use monotonically increasing m/z values then everything should work fine
  ChromatogramExtractorAlgorithm extractor;

  // test the zero first value
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 399.805, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 0.0); // test very first data point

  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.0001, integrated_intensity, 0.001, false);
  TEST_REAL_SIMILAR(integrated_intensity, 8.0);
}
END_SECTION

START_SECTION( [ChromatogramExtractorAlgorithm::ExtractionCoordinates] static bool SortExtractionCoordinatesByMZ(const ChromatogramExtractorAlgorithm::ExtractionCoordinates &left, const ChromatogramExtractorAlgorithm::ExtractionCoordinates &right))    
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION([ChromatogramExtractorAlgorithm::ExtractionCoordinates] static bool SortExtractionCoordinatesReverseByMZ(const ChromatogramExtractorAlgorithm::ExtractionCoordinates &left, const ChromatogramExtractorAlgorithm::ExtractionCoordinates &right))    
{
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

