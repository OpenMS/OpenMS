// Timing harness for ThermoRawFile::transform with an MS-level selection (counting consumer).
// Usage: thermo_time <raw> <level|all> <threads> [profile]
#include <OpenMS/FORMAT/ThermoRawFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace OpenMS;
struct Count : Interfaces::IMSDataConsumer
{
  std::size_t spectra = 0, peaks = 0;
  double sum = 0.0;
  void consumeSpectrum(MSSpectrum& s) override
  {
    ++spectra;
    peaks += s.size();
    for (const auto& p : s) sum += p.getMZ() * 1e-3 + std::log1p(p.getIntensity());
  }
  void consumeChromatogram(MSChromatogram&) override {}
  void setExpectedSize(Size, Size) override {}
  void setExperimentalSettings(const ExperimentalSettings&) override {}
};

int main(int argc, char** argv)
{
  if (argc < 4) return 2;
  omp_set_num_threads(std::atoi(argv[3]));
  ThermoRawFile::Options o;
  o.centroid = !(argc > 4 && std::strcmp(argv[4], "profile") == 0);
  o.checksum = false;
  o.skip_chromatograms = true;
  if (std::strcmp(argv[2], "all") != 0) o.ms_levels = {std::atoi(argv[2])};
  ThermoRawFile f;
  f.setOptions(o);
  Count c;
  const auto t = std::chrono::steady_clock::now();
  f.transform(argv[1], &c);
  const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
  std::printf("%s level=%s threads=%s %s: %.2f s, spectra %zu, peaks %zu, checksum %.6f\n", argv[1], argv[2], argv[3],
              o.centroid ? "centroid" : "profile", s, c.spectra, c.peaks, c.sum);
  return 0;
}
