// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MzPeakFile.h>
///////////////////////////

#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <cmath>
#include <limits>

using namespace OpenMS;
using namespace std;

///////////////////////////

/// The spectra of @p got are those of @p expected (as MzMLFile reads MzMLFile_1.mzML): what the mzPeak writer stores
void compareSpectra(const PeakMap& expected, const PeakMap& got)
{
  TEST_EQUAL(got.getNrSpectra(), expected.getNrSpectra())
  if (got.getNrSpectra() != expected.getNrSpectra()) { return; }
  for (Size i = 0; i < expected.getNrSpectra(); ++i)
  {
    const MSSpectrum& e = expected[i];
    const MSSpectrum& g = got[i];
    TEST_STRING_EQUAL(g.getNativeID(), e.getNativeID())
    TEST_EQUAL(g.getMSLevel(), e.getMSLevel())
    TEST_EQUAL(g.getRT(), e.getRT()) // exact: the writer picks the stored minutes so that the seconds come back exactly
    TEST_EQUAL(static_cast<int>(g.getInstrumentSettings().getPolarity()), static_cast<int>(e.getInstrumentSettings().getPolarity()))
    TEST_EQUAL(g.size(), e.size())
    if (g.size() != e.size()) { continue; }
    for (Size k = 0; k < e.size(); ++k)
    {
      TEST_EQUAL(g[k].getMZ(), e[k].getMZ())
      TEST_EQUAL(g[k].getIntensity(), e[k].getIntensity())
    }
    TEST_EQUAL(g.getPrecursors().size(), e.getPrecursors().size())
    for (Size k = 0; k < std::min(g.getPrecursors().size(), e.getPrecursors().size()); ++k)
    {
      TEST_EQUAL(g.getPrecursors()[k].getMZ(), e.getPrecursors()[k].getMZ())
      TEST_EQUAL(g.getPrecursors()[k].getCharge(), e.getPrecursors()[k].getCharge())
      TEST_REAL_SIMILAR(g.getPrecursors()[k].getIntensity(), e.getPrecursors()[k].getIntensity())
      TEST_REAL_SIMILAR(g.getPrecursors()[k].getIsolationWindowLowerOffset(), e.getPrecursors()[k].getIsolationWindowLowerOffset())
      TEST_REAL_SIMILAR(g.getPrecursors()[k].getIsolationWindowUpperOffset(), e.getPrecursors()[k].getIsolationWindowUpperOffset())
    }
    TEST_EQUAL(g.getDriftTime(), e.getDriftTime())
    if (e.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET) { TEST_EQUAL(static_cast<int>(g.getDriftTimeUnit()), static_cast<int>(e.getDriftTimeUnit())) }
  }
}

START_TEST(MzPeakFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzPeakFile* ptr = nullptr;
MzPeakFile* nullPointer = nullptr;
START_SECTION((MzPeakFile()))
  ptr = new MzPeakFile;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzPeakFile()))
  delete ptr;
END_SECTION

START_SECTION((static bool isMzPeak(const std::string& path)))
{
  TEST_EQUAL(MzPeakFile::isMzPeak(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML")), false)
  TEST_EQUAL(MzPeakFile::isMzPeak(OPENMS_GET_TEST_DATA_PATH("does_not_exist.mzpeak")), false)
  TEST_EQUAL(MzPeakFile::isMzPeak(OPENMS_GET_TEST_DATA_PATH("")), false) // a directory without the index
  TEST_EQUAL(MzPeakFile::isMzPeak(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak")), true)
  TEST_EQUAL(MzPeakFile::isMzPeak(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_small.mzpeak")), true)
  TEST_EQUAL(MzPeakFile::isMzPeak(OPENMS_GET_TEST_DATA_PATH("MultiEntry.mzML.zip")), false) // a ZIP archive, but no mzPeak run
}
END_SECTION

START_SECTION((static bool isSupported()))
{
#ifdef WITH_MZPEAK
  TEST_EQUAL(MzPeakFile::isSupported(), true)
#else
  TEST_EQUAL(MzPeakFile::isSupported(), false)
#endif
}
END_SECTION

#ifdef WITH_MZPEAK

// MzMLFile_1.mzML as MzMLFile reads it: the reference for the files written from it
PeakMap reference;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), reference);

START_SECTION((void load(const std::string& filename, PeakMap& exp)))
{
  // MzPeakFile_1.mzpeak: MzMLFile_1.mzML written with store() (FileConverter)
  PeakMap exp;
  MzPeakFile file;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), exp);
  compareSpectra(reference, exp);
  TEST_EQUAL(exp.getSourceFiles().size(), 1) // (FileConverter rewrote the source file when it wrote the fixture)
  TEST_EQUAL(exp.getLoadedFileType(), FileTypes::MZPEAK)
  // the run's start time, stored as "2007-06-27 15:23:45" (the form store() wrote before it wrote ISO 8601)
  TEST_EQUAL(exp.getDateTime().isValid(), true)
  TEST_STRING_EQUAL(exp.getDateTime().get(), "2007-06-27 15:23:45")
  // the options that select spectra: MS level ...
  file.getOptions().addMSLevel(2);
  file.load(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), exp);
  Size ms2 = 0;
  for (const MSSpectrum& s : reference) { ms2 += s.getMSLevel() == 2; }
  TEST_EQUAL(exp.getNrSpectra(), ms2)
  for (const MSSpectrum& s : exp) { TEST_EQUAL(s.getMSLevel(), 2) }
  // ... RT range (half-open, as for mzML) ...
  file.getOptions() = PeakFileOptions();
  file.getOptions().setRTRange(DRange<1>(reference[1].getRT(), reference[3].getRT()));
  file.load(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), exp);
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_STRING_EQUAL(exp[0].getNativeID(), reference[1].getNativeID())
  // ... and the peaks: an m/z range, metadata only
  file.getOptions() = PeakFileOptions();
  file.getOptions().setMZRange(DRange<1>(reference[0][2].getMZ(), reference[0][6].getMZ()));
  file.load(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), exp);
  TEST_EQUAL(exp.getNrSpectra(), reference.getNrSpectra())
  TEST_EQUAL(exp[0].size(), 4)
  TEST_EQUAL(exp[0][0].getMZ(), reference[0][2].getMZ())
  file.getOptions() = PeakFileOptions();
  file.getOptions().setMetadataOnly(true);
  file.load(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), exp);
  TEST_EQUAL(exp.getNrSpectra(), reference.getNrSpectra())
  TEST_EQUAL(exp[0].size(), 0)
  TEST_EQUAL(exp[0].getRT(), reference[0].getRT())

  // errors
  TEST_EXCEPTION(Exception::FileNotFound, file.load(OPENMS_GET_TEST_DATA_PATH("does_not_exist.mzpeak"), exp))
  TEST_EXCEPTION(Exception::ParseError, file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp))
}
END_SECTION

START_SECTION((void store(const std::string& filename, const PeakMap& exp) const))
{
  std::string tmp;
  NEW_TMP_FILE(tmp);
  tmp += ".mzpeak";
  MzPeakFile file;
  // a spectrum with a drift time, as a timsTOF frame of an mzML file carries it
  PeakMap with_im = reference;
  with_im[0].setDriftTime(0.9876543);
  with_im[0].setDriftTimeUnit(DriftTimeUnit::VSSC);
  with_im[1].setDriftTime(-35.0);
  with_im[1].setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  file.store(tmp, with_im);
  TEST_EQUAL(File::exists(tmp), true)
  TEST_EQUAL(MzPeakFile::isMzPeak(tmp), true)
  PeakMap exp;
  file.load(tmp, exp);
  compareSpectra(with_im, exp);
  TEST_EQUAL(exp[0].getDriftTime(), 0.9876543)
  TEST_EQUAL(static_cast<int>(exp[0].getDriftTimeUnit()), static_cast<int>(DriftTimeUnit::VSSC))
  TEST_EQUAL(exp[1].getDriftTime(), -35.0)
  TEST_EQUAL(static_cast<int>(exp[1].getDriftTimeUnit()), static_cast<int>(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE))
  TEST_EQUAL(exp[2].getDriftTime(), IMTypes::DRIFTTIME_NOT_SET)
  // retention times come back exactly wherever a stored minutes value exists whose * 60 is the seconds (for about 4 % of
  // values none does, the format storing minutes: those differ by one ulp, far below the 15 digits OpenMS writes)
  PeakMap rts;
  for (Size i = 0; i < 1000; ++i)
  {
    MSSpectrum s;
    s.setRT(0.1 + 1.7 * static_cast<double>(i));
    s.setNativeID("scan=" + std::to_string(i));
    s.emplace_back(100.0 + static_cast<double>(i), 1.0f); // (a run whose spectra all have no peaks has no data table and reads as empty)
    rts.addSpectrum(s);
  }
  file.store(tmp, rts);
  file.load(tmp, exp);
  TEST_EQUAL(exp.getNrSpectra(), 1000)
  Size exact = 0, within_ulp = 0;
  for (Size i = 0; i < exp.getNrSpectra(); ++i)
  {
    exact += exp[i].getRT() == rts[i].getRT();
    within_ulp += std::abs(exp[i].getRT() - rts[i].getRT()) <= 2.0 * std::numeric_limits<double>::epsilon() * rts[i].getRT();
  }
  TEST_EQUAL(exact >= 950, true) // 962 with this series; 750 or so with a plain division
  TEST_EQUAL(within_ulp, 1000)   // (the rest within 2 ulps)
  // a run whose spectra have no peaks keeps its spectra (the writer emits an empty peaks table for it)
  PeakMap empty_spectra;
  {
    MSSpectrum s;
    s.setRT(10.0);
    s.setNativeID("scan=1");
    s.setMSLevel(2);
    Precursor prec;
    prec.setMZ(500.25);
    prec.setCharge(2);
    s.getPrecursors().push_back(prec);
    empty_spectra.addSpectrum(s);
    empty_spectra.addSpectrum(s);
    empty_spectra.getSpectra().back().setNativeID("scan=2");
  }
  file.store(tmp, empty_spectra);
  file.load(tmp, exp);
  TEST_EQUAL(exp.getNrSpectra(), 2)
  ABORT_IF(exp.getNrSpectra() != 2)
  TEST_EQUAL(exp[1].size(), 0)
  TEST_STRING_EQUAL(exp[1].getNativeID(), "scan=2")
  TEST_EQUAL(exp[1].getPrecursors().size(), 1)
  // an isolation target that differs from the selected ion m/z survives (as the mzML writer keeps them apart)
  PeakMap distinct = reference;
  Size with_precursor = 0;
  for (MSSpectrum& s : distinct)
  {
    if (s.getPrecursors().empty()) { continue; }
    s.getPrecursors()[0].setMZ(500.123456789);
    s.getPrecursors()[0].setMetaValue("isolation window target m/z", 500.0);
    ++with_precursor;
  }
  TEST_EQUAL(with_precursor > 0, true)
  file.store(tmp, distinct);
  file.load(tmp, exp);
  for (Size i = 0; i < exp.getNrSpectra(); ++i)
  {
    if (exp[i].getPrecursors().empty()) { continue; }
    TEST_EQUAL(exp[i].getPrecursors()[0].getMZ(), 500.123456789)
    TEST_EQUAL(exp[i].getPrecursors()[0].metaValueExists("isolation window target m/z"), true)
    TEST_REAL_SIMILAR(double(exp[i].getPrecursors()[0].getMetaValue("isolation window target m/z")), 500.0)
  }
  // the run's start time and its source file: checksum, file format and native ID format
  PeakMap with_source = reference;
  DateTime start;
  start.set("2007-06-27T15:23:45");
  with_source.setDateTime(start);
  SourceFile source;
  source.setNameOfFile("run.RAW");
  source.setPathToFile("file:///data");
  source.setFileType("Thermo RAW format");
  source.setNativeIDType("Thermo nativeID format");
  source.setNativeIDTypeAccession("MS:1000768");
  source.setChecksum("0123456789abcdef0123456789abcdef01234567", SourceFile::ChecksumType::SHA1);
  with_source.getSourceFiles() = {source};
  file.store(tmp, with_source);
  file.load(tmp, exp);
  TEST_STRING_EQUAL(exp.getDateTime().get(), "2007-06-27 15:23:45")
  TEST_EQUAL(exp.getSourceFiles().size(), 1)
  ABORT_IF(exp.getSourceFiles().size() != 1)
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNameOfFile(), "run.RAW")
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), "0123456789abcdef0123456789abcdef01234567")
  TEST_EQUAL(static_cast<int>(exp.getSourceFiles()[0].getChecksumType()), static_cast<int>(SourceFile::ChecksumType::SHA1))
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getFileType(), "Thermo RAW format")
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNativeIDType(), "Thermo nativeID format")
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNativeIDTypeAccession(), "MS:1000768")
  // a failed write leaves the archive that was there (the new one is written beside it and moved when finished)
  file.store(tmp, reference);
  PeakMap bad = reference;
  bad[0][0].setMZ(std::numeric_limits<double>::quiet_NaN());
  TEST_EXCEPTION(Exception::UnableToCreateFile, file.store(tmp, bad))
  file.load(tmp, exp);
  compareSpectra(reference, exp);
  // via FileHandler
  FileHandler().storeExperiment(tmp, reference, {FileTypes::MZPEAK});
  FileHandler().loadExperiment(tmp, exp, {FileTypes::MZPEAK});
  compareSpectra(reference, exp);
}
END_SECTION

START_SECTION((void transform(const std::string& filename, Interfaces::IMSDataConsumer* consumer)))
{
  PeakMap exp;
  MSDataTransformingConsumer consumer;
  consumer.setSpectraProcessingFunc([&exp](MSSpectrum& s) { exp.addSpectrum(s); });
  MzPeakFile file;
  file.transform(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), &consumer);
  compareSpectra(reference, exp);
  // with options, and through FileHandler (which forwards to the consumer)
  exp.clear(true);
  FileHandler fh;
  fh.getOptions().addMSLevel(1);
  fh.transformExperiment(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), &consumer, {FileTypes::MZPEAK});
  Size ms1 = 0;
  for (const MSSpectrum& s : reference) { ms1 += s.getMSLevel() == 1; }
  TEST_EQUAL(exp.getNrSpectra(), ms1)
  TEST_EXCEPTION(Exception::NullPointer, file.transform(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), nullptr))
}
END_SECTION

START_SECTION([EXTRA] FileHandler type detection and loading)
{
  TEST_EQUAL(FileHandler::getTypeByFileName("x.mzpeak"), FileTypes::MZPEAK)
  TEST_EQUAL(FileHandler::getType(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak")), FileTypes::MZPEAK)
  TEST_EQUAL(FileHandler::getType(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_small.mzpeak")), FileTypes::MZPEAK)
  TEST_STRING_EQUAL(FileTypes::typeToName(FileTypes::MZPEAK), "mzpeak")
  PeakMap exp;
  FileHandler fh;
  fh.getOptions().addMSLevel(2);
  fh.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_1.mzpeak"), exp, {FileTypes::MZPEAK});
  Size ms2 = 0;
  for (const MSSpectrum& s : reference) { ms2 += s.getMSLevel() == 2; }
  TEST_EQUAL(exp.getNrSpectra(), ms2)
  TEST_EQUAL(exp.getSourceFiles().size(), 1) // rewritten by FileHandler
}
END_SECTION

START_SECTION([EXTRA] a file of the reference implementation)
{
  // MzPeakFile_small.mzpeak: test/files/small.mzpeak of https://github.com/OpenMS/mzpeak (written by the Rust reference
  // implementation from the small.mzML of ProteoWizard: 48 profile spectra of a Thermo LTQ FT, 5 s to 9 min). The
  // values are those the library's own tests (spectra_test.cpp) check.
  PeakMap exp;
  MzPeakFile file;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_small.mzpeak"), exp);
  TEST_EQUAL(exp.getNrSpectra(), 48)
  ABORT_IF(exp.getNrSpectra() != 48)
  // (the reference writer drops the zero-intensity points of the profile spectra: 13589 of the 19913 of spectrum 0 remain)
  TEST_EQUAL(exp[0].size(), 13589)
  TEST_EQUAL(exp[1].size(), 18177)
  TEST_EQUAL(exp[7].size(), 10329)
  TEST_REAL_SIMILAR(exp[0][8].getMZ(), 204.7593349)
  TEST_EQUAL(exp[0].getMSLevel(), 1)
  TEST_REAL_SIMILAR(exp[0].getRT(), 0.2961) // seconds: the file stores 0.004935 minutes
  TEST_EQUAL(static_cast<int>(exp[0].getType()), static_cast<int>(SpectrumSettings::SpectrumType::PROFILE))
  TEST_EQUAL(static_cast<int>(exp[0].getInstrumentSettings().getPolarity()), static_cast<int>(IonSource::Polarity::POSITIVE))
  TEST_STRING_EQUAL(exp[0].getNativeID(), "controllerType=0 controllerNumber=1 scan=1")
  Size ms2 = 0;
  for (const MSSpectrum& s : exp) { ms2 += s.getMSLevel() == 2; }
  TEST_EQUAL(ms2, 34)
  // the first MS2 spectrum (scan=3): CID at 35 eV of 810.79 m/z, isolation window +-1
  TEST_EQUAL(exp[2].getMSLevel(), 2)
  TEST_STRING_EQUAL(exp[2].getNativeID(), "controllerType=0 controllerNumber=1 scan=3")
  TEST_REAL_SIMILAR(exp[2].getRT(), 0.011218333333 * 60.0)
  TEST_EQUAL(exp[2].getPrecursors().size(), 1)
  ABORT_IF(exp[2].getPrecursors().empty())
  TEST_REAL_SIMILAR(exp[2].getPrecursors()[0].getMZ(), 810.79)
  TEST_REAL_SIMILAR(exp[2].getPrecursors()[0].getIntensity(), 1.994039182404e06)
  TEST_REAL_SIMILAR(exp[2].getPrecursors()[0].getIsolationWindowLowerOffset(), 1.0)
  TEST_REAL_SIMILAR(exp[2].getPrecursors()[0].getIsolationWindowUpperOffset(), 1.0)
  TEST_EQUAL(exp[2].getPrecursors()[0].getActivationMethods().count(Precursor::ActivationMethod::CID), 1)
  // MS:1000045 collision energy is the meta value "collision energy", as MzMLHandler keeps it (the member is MS:1000509)
  TEST_EQUAL(exp[2].getPrecursors()[0].metaValueExists("collision energy"), true)
  TEST_REAL_SIMILAR(double(exp[2].getPrecursors()[0].getMetaValue("collision energy")), 35.0)
  TEST_REAL_SIMILAR(exp[2].getPrecursors()[0].getActivationEnergy(), 0.0)
  TEST_EQUAL(exp.getChromatograms().size(), 1)
  ABORT_IF(exp.getChromatograms().empty())
  TEST_STRING_EQUAL(exp.getChromatograms()[0].getNativeID(), "TIC")
  TEST_EQUAL(exp.getChromatograms()[0].size(), 48)
  TEST_EQUAL(static_cast<int>(exp.getChromatograms()[0].getChromatogramType()), static_cast<int>(ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM))
  TEST_REAL_SIMILAR(exp.getChromatograms()[0][0].getRT(), 0.2961)
  // the run's start time (ISO 8601 with zone: "2005-07-20T19:44:22Z") and source file (by its CV terms)
  TEST_STRING_EQUAL(exp.getDateTime().get(), "2005-07-20 19:44:22")
  TEST_EQUAL(exp.getSourceFiles().size(), 1)
  ABORT_IF(exp.getSourceFiles().size() != 1)
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNameOfFile(), "small.RAW")
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getFileType(), "Thermo RAW format")
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNativeIDType(), "Thermo nativeID format")
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNativeIDTypeAccession(), "MS:1000768")
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), "b43e9286b40e8b5dbc0dfa2e428495769ca96a96")
  // FileHandler with rewrite_source_file (and no hash): the source file is the archive, without the checksum recorded
  // for small.RAW
  FileHandler().loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzPeakFile_small.mzpeak"), exp, {FileTypes::MZPEAK}, ProgressLogger::NONE, true, false);
  TEST_EQUAL(exp.getSourceFiles().size(), 1)
  ABORT_IF(exp.getSourceFiles().size() != 1)
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNameOfFile(), "MzPeakFile_small.mzpeak")
  TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), "")
}
END_SECTION

#endif // WITH_MZPEAK

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
