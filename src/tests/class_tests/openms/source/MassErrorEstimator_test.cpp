#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/PROCESSING/CALIBRATION/MassErrorEstimator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>

using namespace OpenMS;

START_TEST(MassErrorEstimator, "$Id$")

/////////////////////////////////////////////////////////////

MassErrorEstimator* ptr = nullptr;
MassErrorEstimator* nullPointer = nullptr;

START_SECTION((MassErrorEstimator())) 
{
  ptr = new MassErrorEstimator();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~MassErrorEstimator())) 
{
  delete ptr;
}
END_SECTION

START_SECTION((double calculateMean(const std::vector<MarkerPeak>&, double)))  // checking that mean = 0 when all peaks are at reference m/z
{
  MassErrorEstimator est;

  std::vector<MarkerPeak> peaks;
  peaks.push_back({445.12003, 10.0, 100.0});
  peaks.push_back({445.12003, 20.0, 200.0});

  double mean = est.calculateMean(peaks, 445.12003);

  TEST_REAL_SIMILAR(mean, 0.0)
}
END_SECTION

START_SECTION((double calculateSD(const std::vector<MarkerPeak>&, double))) // checking that sd = 0 when all mass errors are the same
{
  MassErrorEstimator est;

  std::vector<MarkerPeak> peaks;
  peaks.push_back({445.12003, 10.0, 100.0});
  peaks.push_back({445.12003, 20.0, 200.0});

  double sd = est.calculateSD(peaks, 445.12003);

  TEST_REAL_SIMILAR(sd, 0.0)
}
END_SECTION

START_SECTION((std::vector<MarkerPeak> RTdistribution(const std::vector<MarkerPeak>&))) // checking if points in distribution are returned if input is not empty
{
  MassErrorEstimator est;

  std::vector<MarkerPeak> candidates;
  candidates.push_back({445.12003, 10.0, 100.0});
  candidates.push_back({445.12004, 20.0, 200.0});
  candidates.push_back({445.12005, 30.0, 300.0});

  std::vector<MarkerPeak> selected = est.RTdistribution(candidates);

  TEST_EQUAL(selected.empty(), false)
}
END_SECTION

START_SECTION(std::vector<MarkerPeak> findPolysiloxaneCandidates(const MSExperiment&)) // checking that a strong polysiloxane peak is found with low-intensity background
{
  MassErrorEstimator est;

  MSExperiment exp;
  MSSpectrum spec;
  spec.setMSLevel(1);
  spec.setRT(100.0);

  for (Size i = 0; i < 20; ++i) // background peaks added for median-based intensity filter
  {
    Peak1D noise;
    noise.setMZ(444.2 + i * 0.05);
    noise.setIntensity(10.0);
    spec.push_back(noise);
  }

  Peak1D signal; // one strong polysiloxane peak added at reference m/z
  signal.setMZ(445.12003);
  signal.setIntensity(1000.0);
  spec.push_back(signal);

  exp.addSpectrum(spec);

  std::vector<MarkerPeak> result = est.findPolysiloxaneCandidates(exp);

  TEST_EQUAL(result.size(), 1)
  TEST_REAL_SIMILAR(result[0].mz, 445.12003)
  TEST_REAL_SIMILAR(result[0].rt, 100.0)
}
END_SECTION

// Edge cases

START_SECTION((calculateMean() with empty input)) // checking if algorithm doesn't crash when calculateMean has an empty input --> /0 -> mathematically undefined behaviour
{
  MassErrorEstimator est;
  std::vector<MarkerPeak> peaks;

  double mean = est.calculateMean(peaks, 445.12003);

  TEST_REAL_SIMILAR(mean, 0.0)
}
END_SECTION

START_SECTION((calculateSD() with one point)) // checking if algorithm doesn't crash when calculateSD has only 1 input point --> n-1 turns into 0 so /0 -> undefined behaviour
{
  MassErrorEstimator est;

  std::vector<MarkerPeak> peaks;
  peaks.push_back({445.12003, 10.0, 100.0});

  double sd = est.calculateSD(peaks, 445.12003);

  TEST_REAL_SIMILAR(sd, 0.0)
}
END_SECTION

START_SECTION((RTdistribution() with empty input)) // checking that RTdistribution doesn't crash when input is empty 
{
  MassErrorEstimator est;
  std::vector<MarkerPeak> candidates;

  std::vector<MarkerPeak> selected = est.RTdistribution(candidates);

  TEST_EQUAL(selected.empty(), true)
}
END_SECTION

START_SECTION(findPolysiloxaneCandidates ignores MS2 spectra) // checking that algorithm only works with MS1 spectra
{
  MassErrorEstimator est;

  MSExperiment exp;
  MSSpectrum spec;
  spec.setMSLevel(2);
  spec.setRT(100.0);

  Peak1D signal;
  signal.setMZ(445.12003);
  signal.setIntensity(1000.0);
  spec.push_back(signal);

  exp.addSpectrum(spec);

  std::vector<MarkerPeak> result = est.findPolysiloxaneCandidates(exp);

  TEST_EQUAL(result.empty(), true)
}
END_SECTION
/////////////////////////////////////////////////////////////

END_TEST
