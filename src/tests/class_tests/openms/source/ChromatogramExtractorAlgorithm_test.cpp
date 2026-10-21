// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>

#include <algorithm>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;


void find_max_helper(const OpenSwath::ChromatogramPtr& chrom, double &max_value, double &foundat)
{
  max_value = -1;
  foundat = -1;
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
}

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

START_SECTION(void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, std::vector< OpenSwath::ChromatogramPtr > &output, std::vector< ExtractionCoordinates >& extraction_coordinates, double mz_extraction_window, bool ppm, std::string filter))
{
  double extract_window = 0.05;
  std::shared_ptr<PeakMap > exp(new PeakMap);
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
  extractor.extractChromatograms(expptr, out_exp, coordinates, extract_window, false, -1, "tophat");

  OpenSwath::ChromatogramPtr chrom = out_exp[0];

  TEST_EQUAL(chrom->getTimeArray()->data.size(), 59);
  TEST_EQUAL(chrom->getIntensityArray()->data.size(), 59);
  // we sort/reorder 
  int firstchromat  = 1;
  int secondchromat = 2;
  int thirdchromat  = 0;

  double max_value = -1; double foundat = -1;
  find_max_helper(out_exp[firstchromat], max_value, foundat);
  TEST_REAL_SIMILAR(max_value, 169.792);
  TEST_REAL_SIMILAR(foundat, 3120.26);

  find_max_helper(out_exp[secondchromat], max_value, foundat);
  TEST_REAL_SIMILAR(max_value, 577.33);
  TEST_REAL_SIMILAR(foundat, 3120.26);

  find_max_helper(out_exp[thirdchromat], max_value, foundat);
  TEST_REAL_SIMILAR(max_value, 35.593);
  TEST_REAL_SIMILAR(foundat, 3055.16);

  // there is no ion mobility, so this should not work
  TEST_EXCEPTION(Exception::IllegalArgument, extractor.extractChromatograms(expptr, out_exp, coordinates, extract_window, false, 1, "tophat"))
}
END_SECTION

START_SECTION([EXTRA] void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, std::vector< OpenSwath::ChromatogramPtr > &output, std::vector< ExtractionCoordinates >& extraction_coordinates, double mz_extraction_window, bool ppm, std::string filter))
{
  typedef OpenMS::DataArrays::FloatDataArray FloatDataArray;

  double extract_window = 0.10;
  std::shared_ptr<PeakMap > exp(new PeakMap);
  // MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.mzML"), *exp);
  for (int i = 0; i < 4; i++)
  {
    MSSpectrum s;
    s.setRT(i);
    FloatDataArray fda;
    for (int k = 0; k < 10; k++)
    {
      Peak1D p;
      p.setMZ(618.3 + k * 0.01);
      p.setIntensity(100 * i + k *2);
      s.push_back(p);
      fda.push_back(100 + k *10);
    }
    for (int k = 0; k < 10; k++)
    {
      Peak1D p;
      p.setMZ(628.4 + k * 0.01);
      p.setIntensity(100 * i + k*2);
      s.push_back(p);
      fda.push_back(100 + k *10);
      std::cout << " ion mobility  " << 100 + k*10 << " : " << p << std::endl;
    }
    fda.setName(Constants::UserParam::ION_MOBILITY);
    s.getFloatDataArrays().push_back(fda);
    exp->addSpectrum(s);
  }
  OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

  ///////////////////////////////////////////////////////////////////////////
  ChromatogramExtractorAlgorithm extractor;

  std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > coordinates;
  {
    ChromatogramExtractorAlgorithm::ExtractionCoordinates coord;
    coord.mz = 618.31; coord.rt_start = 0; coord.rt_end = -1; coord.id = "tr1"; coord.ion_mobility = 120;
    coordinates.push_back(coord);
    coord.mz = 628.45; coord.rt_start = 0; coord.rt_end = -1; coord.id = "tr2"; coord.ion_mobility = 170;
    coordinates.push_back(coord);
  }

  // no IM window
  {
    std::vector< OpenSwath::ChromatogramPtr > out_exp;
    for (int i = 0; i < 2; i++)
    {
      OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
      out_exp.push_back(s);
    }

    extractor.extractChromatograms(expptr, out_exp, coordinates, extract_window, false, -1, "tophat");
    OpenSwath::ChromatogramPtr chrom = out_exp[0];

    TEST_EQUAL(chrom->getTimeArray()->data.size(), 4);
    TEST_EQUAL(chrom->getIntensityArray()->data.size(), 4);

    double max_value = -1; double foundat = -1;
    find_max_helper(out_exp[0], max_value, foundat);
    TEST_REAL_SIMILAR(max_value, 1830)
    TEST_REAL_SIMILAR(foundat, 3)

    find_max_helper(out_exp[1], max_value, foundat);
    TEST_REAL_SIMILAR(max_value, 2790)
    TEST_REAL_SIMILAR(foundat, 3)
  }

  // small IM window
  {
    std::vector< OpenSwath::ChromatogramPtr > out_exp;
    for (int i = 0; i < 2; i++)
    {
      OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
      out_exp.push_back(s);
    }

    extractor.extractChromatograms(expptr, out_exp, coordinates, extract_window, false, 15, "tophat");
    OpenSwath::ChromatogramPtr chrom = out_exp[0];

    TEST_EQUAL(chrom->getTimeArray()->data.size(), 4);
    TEST_EQUAL(chrom->getIntensityArray()->data.size(), 4);

    double max_value = -1; double foundat = -1;
    find_max_helper(out_exp[0], max_value, foundat);
    TEST_REAL_SIMILAR(max_value, 304)
    TEST_REAL_SIMILAR(foundat, 3)

    find_max_helper(out_exp[1], max_value, foundat);
    TEST_REAL_SIMILAR(max_value, 314)
    TEST_REAL_SIMILAR(foundat, 3)
  }

  // larger IM window
  {
    std::vector< OpenSwath::ChromatogramPtr > out_exp;
    for (int i = 0; i < 2; i++)
    {
      OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
      out_exp.push_back(s);
    }

    extractor.extractChromatograms(expptr, out_exp, coordinates, extract_window, false, 30, "tophat");
    OpenSwath::ChromatogramPtr chrom = out_exp[0];

    TEST_EQUAL(chrom->getTimeArray()->data.size(), 4);
    TEST_EQUAL(chrom->getIntensityArray()->data.size(), 4);

    double max_value = -1; double foundat = -1;
    find_max_helper(out_exp[0], max_value, foundat);
    TEST_REAL_SIMILAR(max_value, 303 + 304 + 305)
    TEST_REAL_SIMILAR(foundat, 3)

    find_max_helper(out_exp[1], max_value, foundat);
    TEST_REAL_SIMILAR(max_value, 313 + 314 + 315)
    TEST_REAL_SIMILAR(foundat, 3)
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////
/// Private functions
///////////////////////////////////////////////////////////////////////////

//  mz_a = [400+0.01*i for i in range(20)]
//  int_a = [0 + i*100.0 for i in range(10)] + [900 - i*100.0 for i in range(10)]
//  im_a = [100+0.01*i +100*(i%2) for i in range(20)]
//  zip_a = [ (a,b,c) for a,b,c in zip(mz_a, int_a, im_a) ]
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
  // additional values
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
  // additional values
  10.0, 
  10.0, 
};
static const double im_arr[] = {
  100.0,
  200.01,
  100.02,
  200.03,
  100.04,
  200.05,
  100.06,
  200.07,
  100.08,
  200.09,
  100.1,
  200.11,
  100.12,
  200.13,
  100.14,
  200.15,
  100.16,
  200.17,
  100.18,
  200.19,
  // additional values
  300.1,
  300.2
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
  //print(sum([0 + i*100.0 for i in range(10)]) + sum([900 - i*100.0 for i in range(6)]) + 8 )
  TEST_REAL_SIMILAR( integrated_intensity, 8408.0);
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
  TEST_REAL_SIMILAR( integrated_intensity,8408.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, 400.1, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,9008.0); // 500 ppm is slightly more than 0.1 Da: includes the very first value

}
END_SECTION

START_SECTION([EXTRA IM]void extract_value_tophat(const std::vector< double >::const_iterator &mz_start, std::vector< double >::const_iterator &mz_it, const std::vector< double >::const_iterator &mz_end, std::vector< double >::const_iterator &int_it, const double &mz, double &integrated_intensity, const double &mz_extraction_window, bool ppm))
{ 
  std::vector<double> mz (mz_arr, mz_arr + sizeof(mz_arr) / sizeof(mz_arr[0]) );
  std::vector<double> intensities (int_arr, int_arr + sizeof(int_arr) / sizeof(int_arr[0]) );
  std::vector<double> ion_mobility (im_arr, im_arr + sizeof(im_arr) / sizeof(im_arr[0]) );

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
  std::vector<double>::const_iterator im_it = ion_mobility.begin();

  double integrated_intensity = 0;
  double extract_window = 0.2; // +/- 0.1
  double im_extract_window = 0.3; // +/- 0.15

  // If we use monotonically increasing m/z values then everything should work fine
  ChromatogramExtractorAlgorithm extractor;

  // extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 399.805, 100, integrated_intensity, extract_window, im_extract_window, false);
  // test the zero first value
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 399.805, 100, integrated_intensity, extract_window, im_extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 0.0); // test very first data point

  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 399.91, 100, integrated_intensity, extract_window, im_extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 8.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 400.0, 100, integrated_intensity, extract_window, im_extract_window, false);
  // sum([i for m,i,im in zip_a if im < 100.15 and m < 400.1]) + 8
  TEST_REAL_SIMILAR( integrated_intensity, 2008.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 400.05,  100, integrated_intensity, extract_window, im_extract_window, false);
  // sum([i for m,i,im in zip_a if im < 100.15 and m < 400.15]) + 8
  TEST_REAL_SIMILAR( integrated_intensity, 4108.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 400.1, 100, integrated_intensity, extract_window, im_extract_window, false);
  // sum([i for m,i,im in zip_a if im < 100.15 and m < 400.2])
  TEST_REAL_SIMILAR( integrated_intensity, 4100.0);
  TEST_EQUAL((int)integrated_intensity, 4100);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 400.28, 100, integrated_intensity, extract_window, im_extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity, 0.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 400.28, 200, integrated_intensity, extract_window, im_extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity, 0.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 400.28, 200.1, integrated_intensity, extract_window, im_extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity, 0.0);

  // test the very last value
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 500.0, 300.0, integrated_intensity, extract_window, im_extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity, 0.0);
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 500.0, 300.1, integrated_intensity, extract_window, im_extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity, 10.0);

  // this is to document the situation of using m/z values that are not monotonically increasing:
  //  --> it might not give the correct result (9000) if we try to extract 400.1 AFTER 500.0 
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 400.1, 100, integrated_intensity, extract_window, im_extract_window, false);
  TEST_NOT_EQUAL((int)integrated_intensity, 9000);
}
END_SECTION

START_SECTION([EXTRA]void extract_value_tophat(const std::vector< double >::const_iterator &mz_start, std::vector< double >::const_iterator &mz_it, const std::vector< double >::const_iterator &mz_end, std::vector< double >::const_iterator &int_it, const double &mz, double &integrated_intensity, const double &mz_extraction_window, bool ppm))
{ 
  std::vector<double> mz (mz_arr, mz_arr + sizeof(mz_arr) / sizeof(mz_arr[0]) );
  std::vector<double> intensities (int_arr, int_arr + sizeof(int_arr) / sizeof(int_arr[0]) );

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

START_SECTION([EXTRA]void extract_value_tophat(const std::vector< double >::const_iterator &mz_start, std::vector< double >::const_iterator &mz_it, const std::vector< double >::const_iterator &mz_end, std::vector< double >::const_iterator &int_it, const double &mz, double &integrated_intensity, const double &mz_extraction_window, bool ppm))
{ 
  std::vector<double> mz (mz_arr, mz_arr + sizeof(mz_arr) / sizeof(mz_arr[0]) );
  std::vector<double> intensities (int_arr, int_arr + sizeof(int_arr) / sizeof(int_arr[0]) );
  std::vector<double> ion_mobility (im_arr, im_arr + sizeof(im_arr) / sizeof(im_arr[0]) );

  std::vector<double>::const_iterator mz_start = mz.begin();
  std::vector<double>::const_iterator mz_it_end = mz.end();
  std::vector<double>::const_iterator mz_it = mz.begin();
  std::vector<double>::const_iterator int_it = intensities.begin();
  std::vector<double>::const_iterator im_it = ion_mobility.begin();

  double integrated_intensity = 0;
  double extract_window = 0.2; // +/- 0.1
  double im_extract_window = 0.3; // +/- 0.15

  // If we use monotonically increasing m/z values then everything should work fine
  ChromatogramExtractorAlgorithm extractor;

  // test the zero first value
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 399.805, 100, integrated_intensity, extract_window, im_extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 0.0); // test very first data point
  extractor.extract_value_tophat(mz_start, mz_it, mz_it_end, int_it, im_it, 400.0001, 100, integrated_intensity, 0.001, im_extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 8.0);
}
END_SECTION

START_SECTION([EXTRA]void extract_value_tophat(const std::vector< double >::const_iterator &mz_start, std::vector< double >::const_iterator &mz_it, const std::vector< double >::const_iterator &mz_end, std::vector< double >::const_iterator &int_it, const double &mz, double &integrated_intensity, const double &mz_extraction_window, bool ppm))
{
  // Targets far ahead of the current position (e.g. the precursors of one isolation window in a
  // full MS1 spectrum): the iterators have to end up where a peak-by-peak walk would leave them.
  std::vector<double> mz, intensities, ion_mobility;
  for (Size i = 0; i < 10000; ++i)
  {
    mz.push_back(300.0 + 0.1 * static_cast<double>(i));
    intensities.push_back(static_cast<double>(i % 7 + 1));
    ion_mobility.push_back(100.0);
  }
  // some duplicated m/z values around a target
  mz[5001] = mz[5000];
  mz[5002] = mz[5000];

  // (the last target lies beyond the spectrum, with the last peak outside of its window)
  const double targets[] = {300.0, 300.35, 420.04, 800.0, 800.0, 800.05, 1299.9, 1500.0};
  ChromatogramExtractorAlgorithm extractor;

  std::vector<double>::const_iterator mz_it = mz.begin();
  std::vector<double>::const_iterator int_it = intensities.begin();
  std::vector<double>::const_iterator mz_it_im = mz.begin();
  std::vector<double>::const_iterator int_it_im = intensities.begin();
  std::vector<double>::const_iterator im_it = ion_mobility.begin();
  for (const double target : targets)
  {
    // reference: first peak that is not below the target, and the summed peaks inside the window
    const std::ptrdiff_t expected_pos = std::lower_bound(mz.begin(), mz.end(), target) - mz.begin();
    double expected_intensity = 0;
    for (Size i = 0; i < mz.size(); ++i)
    {
      if (mz[i] > target - 0.1 && mz[i] < target + 0.1) expected_intensity += intensities[i];
    }

    double integrated_intensity = 0;
    extractor.extract_value_tophat(mz.begin(), mz_it, mz.end(), int_it, target, integrated_intensity, 0.2, false);
    TEST_EQUAL(mz_it - mz.begin(), expected_pos)
    TEST_EQUAL(int_it - intensities.begin(), expected_pos)
    TEST_REAL_SIMILAR(integrated_intensity, expected_intensity)

    extractor.extract_value_tophat(mz.begin(), mz_it_im, mz.end(), int_it_im, im_it, target, 100, integrated_intensity, 0.2, 0.3, false);
    TEST_EQUAL(mz_it_im - mz.begin(), expected_pos)
    TEST_EQUAL(int_it_im - intensities.begin(), expected_pos)
    TEST_EQUAL(im_it - ion_mobility.begin(), expected_pos)
    TEST_REAL_SIMILAR(integrated_intensity, expected_intensity)
  }
}
END_SECTION

START_SECTION([EXTRA]void extract_value_tophat(const std::vector< double >::const_iterator &mz_start, std::vector< double >::const_iterator &mz_it, const std::vector< double >::const_iterator &mz_end, std::vector< double >::const_iterator &int_it, const double &mz, double &integrated_intensity, const double &mz_extraction_window, bool ppm))
{
  // The first and the last peak of a spectrum count once if they are inside the window: the first
  // one also when the walk to the left reaches it after several steps, the last one also when the
  // target lies beyond the end of the spectrum. Both overloads, window of 1 Th (+/- 0.5).
  ChromatogramExtractorAlgorithm extractor;
  auto extract = [&extractor](const std::vector<double>& mz, const std::vector<double>& intensities, const double target)
  {
    const std::vector<double> ion_mobility(mz.size(), 100.0);
    std::vector<double>::const_iterator mz_it = mz.begin();
    std::vector<double>::const_iterator int_it = intensities.begin();
    double intensity = 0;
    extractor.extract_value_tophat(mz.begin(), mz_it, mz.end(), int_it, target, intensity, 1.0, false);

    std::vector<double>::const_iterator mz_it_im = mz.begin();
    std::vector<double>::const_iterator int_it_im = intensities.begin();
    std::vector<double>::const_iterator im_it = ion_mobility.begin();
    double intensity_im = 0;
    extractor.extract_value_tophat(mz.begin(), mz_it_im, mz.end(), int_it_im, im_it, target, 100.0, intensity_im, 1.0, 0.3, false);
    return std::make_pair(intensity, intensity_im);
  };

  // the walk to the left passes three peaks, the last of them is the first peak of the spectrum
  std::pair<double, double> result = extract({100.0, 100.1, 100.2, 100.3}, {1.0, 10.0, 100.0, 1000.0}, 100.25);
  TEST_REAL_SIMILAR(result.first, 1111.0)
  TEST_REAL_SIMILAR(result.second, 1111.0)

  // target beyond the last peak
  result = extract({100.0, 100.1, 100.2}, {1.0, 10.0, 100.0}, 100.3);
  TEST_REAL_SIMILAR(result.first, 111.0)
  TEST_REAL_SIMILAR(result.second, 111.0)

  // a single peak, target beyond and before it
  result = extract({100.0}, {5.0}, 100.1);
  TEST_REAL_SIMILAR(result.first, 5.0)
  TEST_REAL_SIMILAR(result.second, 5.0)
  result = extract({100.0}, {5.0}, 99.9);
  TEST_REAL_SIMILAR(result.first, 5.0)
  TEST_REAL_SIMILAR(result.second, 5.0)
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

START_SECTION([EXTRA] extractChromatograms of many coordinates with and without ion mobility (in parallel ranges))
{
  // peak i has m/z 100 + i, intensity i + 1 and ion mobility i; every second coordinate has no ion mobility
  const int nr_peaks = 600;
  std::shared_ptr<PeakMap> exp(new PeakMap);
  for (int scan = 0; scan < 3; ++scan)
  {
    MSSpectrum s;
    s.setRT(scan);
    DataArrays::FloatDataArray fda;
    for (int i = 0; i < nr_peaks; ++i)
    {
      s.push_back(Peak1D(100.0 + i, i + 1.0));
      fda.push_back(i);
    }
    fda.setName(Constants::UserParam::ION_MOBILITY);
    s.getFloatDataArrays().push_back(fda);
    exp->addSpectrum(s);
  }
  // spectra in memory: extracted in ranges of the coordinates by several threads
  OpenSwath::SpectrumAccessPtr in_memory(new SpectrumAccessOpenMSInMemory(*SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp)));

  std::vector<ChromatogramExtractorAlgorithm::ExtractionCoordinates> coordinates(nr_peaks);
  for (int i = 0; i < nr_peaks; ++i)
  {
    coordinates[i].mz = 100.0 + i;
    coordinates[i].rt_start = 0;
    coordinates[i].rt_end = -1;
    coordinates[i].ion_mobility = (i % 2 == 0) ? i : -1;
    coordinates[i].id = std::to_string(i);
  }

  auto extract = [&]([[maybe_unused]] int threads)
  {
    std::vector<OpenSwath::ChromatogramPtr> output;
    for (int i = 0; i < nr_peaks; ++i)
    {
      output.emplace_back(new OpenSwath::Chromatogram);
    }
#ifdef _OPENMP
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads(threads);
#endif
    ChromatogramExtractorAlgorithm().extractChromatograms(in_memory, output, coordinates, 0.5, false, 0.5, "tophat");
#ifdef _OPENMP
    omp_set_num_threads(max_threads);
#endif
    return output;
  };

  for (int threads : {1, 4})
  {
    std::vector<OpenSwath::ChromatogramPtr> output = extract(threads);
    Size wrong_chromatograms = 0;
    for (int i = 0; i < nr_peaks; ++i)
    {
      const std::vector<double>& intensities = output[i]->getIntensityArray()->data;
      if (intensities != std::vector<double>(3, i + 1.0))
      {
        ++wrong_chromatograms;
      }
    }
    TEST_EQUAL(wrong_chromatograms, 0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

