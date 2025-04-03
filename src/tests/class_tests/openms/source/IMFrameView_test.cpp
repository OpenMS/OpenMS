// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/IMFrameView.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/METADATA/DataArrays.h>

using namespace OpenMS;
using namespace std;

MSSpectrum createTestSpectrum()
{
  MSSpectrum spectrum;
  
  // Add peaks
  Peak1D peak1;
  peak1.setMZ(100.0);
  peak1.setIntensity(1000.0f);
  spectrum.push_back(peak1);
  
  Peak1D peak2;
  peak2.setMZ(200.0);
  peak2.setIntensity(2000.0f);
  spectrum.push_back(peak2);
  
  Peak1D peak3;
  peak3.setMZ(300.0);
  peak3.setIntensity(3000.0f);
  spectrum.push_back(peak3);
  
  // Set spectrum as centroided
  spectrum.setType(SpectrumSettings::CENTROID);
  
  // Add drift time array using controlled vocabulary
  DataArrays::FloatDataArray drift_times;
  IMDataConverter::setIMUnit(drift_times, DriftTimeUnit::MILLISECOND);
  drift_times.push_back(10.0f);
  drift_times.push_back(20.0f);
  drift_times.push_back(30.0f);
  
  spectrum.getFloatDataArrays().push_back(drift_times);
  
  // Set other metadata
  spectrum.setRT(42.0);
  spectrum.setMSLevel(1);
  spectrum.setNativeID("test_spectrum");
  
  // Verify that the spectrum has the correct IMFormat
  TEST_EQUAL(IMTypes::determineIMFormat(spectrum), IMFormat::CONCATENATED);
  
  return spectrum;
}

// Helper function to create a test spectrum without drift time array
MSSpectrum createInvalidSpectrum()
{
  MSSpectrum spectrum;
  
  // Add peaks
  Peak1D peak1;
  peak1.setMZ(100.0);
  peak1.setIntensity(1000.0f);
  spectrum.push_back(peak1);
  
  // Verify that the spectrum has the correct IMFormat
  TEST_EQUAL(IMTypes::determineIMFormat(spectrum), IMFormat::NONE);
  
  // Set spectrum as centroided
  spectrum.setType(SpectrumSettings::CENTROID);
  
  return spectrum;
}

// Helper function to create a test spectrum with mismatched drift time array size
MSSpectrum createMismatchedSpectrum()
{
  MSSpectrum spectrum;
  
  // Add peaks
  Peak1D peak1;
  peak1.setMZ(100.0);
  peak1.setIntensity(1000.0f);
  spectrum.push_back(peak1);
  
  Peak1D peak2;
  peak2.setMZ(200.0);
  peak2.setIntensity(2000.0f);
  spectrum.push_back(peak2);
  
  // Set spectrum as centroided
  spectrum.setType(SpectrumSettings::CENTROID);
  
  // Add drift time array with wrong size
  DataArrays::FloatDataArray drift_times;
  IMDataConverter::setIMUnit(drift_times, DriftTimeUnit::MILLISECOND);
  drift_times.push_back(10.0f);
  
  spectrum.getFloatDataArrays().push_back(drift_times);
  
  // Verify that the spectrum has the correct IMFormat
  TEST_EQUAL(IMTypes::determineIMFormat(spectrum), IMFormat::CONCATENATED);
  
  return spectrum;
}

// Helper function to create a test spectrum that is not centroided
MSSpectrum createNonCentroidedSpectrum()
{
  MSSpectrum spectrum;
  
  // Add peaks
  Peak1D peak1;
  peak1.setMZ(100.0);
  peak1.setIntensity(1000.0f);
  spectrum.push_back(peak1);
  
  // Set spectrum as profile
  spectrum.setType(SpectrumSettings::PROFILE);
  
  // Add drift time array
  DataArrays::FloatDataArray drift_times;
  IMDataConverter::setIMUnit(drift_times, DriftTimeUnit::MILLISECOND);
  drift_times.push_back(10.0f);
  
  spectrum.getFloatDataArrays().push_back(drift_times);
  
  // Verify that the spectrum has the correct IMFormat
  TEST_EQUAL(IMTypes::determineIMFormat(spectrum), IMFormat::CONCATENATED);
  
  return spectrum;
}

// Helper function to create a test spectrum with MULTIPLE_SPECTRA format
MSSpectrum createMultipleSpectraFormatSpectrum()
{
  MSSpectrum spectrum;
  
  // Add peaks
  Peak1D peak1;
  peak1.setMZ(100.0);
  peak1.setIntensity(1000.0f);
  spectrum.push_back(peak1);
  
  // Set spectrum as centroided
  spectrum.setType(SpectrumSettings::CENTROID);
  
  // Set drift time directly on the spectrum (MULTIPLE_SPECTRA format)
  spectrum.setDriftTime(15.0);
  spectrum.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  
  // Verify that the spectrum has the correct IMFormat
  TEST_EQUAL(IMTypes::determineIMFormat(spectrum), IMFormat::MULTIPLE_SPECTRA);
  
  return spectrum;
}

START_TEST(IMFrameView, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMFrameView* imframe_ptr = nullptr;
IMFrameView* nullPointer = nullptr;

START_SECTION((explicit IMFrameView(MSSpectrum* spectrum_ptr)))
{
  MSSpectrum spectrum = createTestSpectrum();
  imframe_ptr = new IMFrameView(&spectrum);
  TEST_NOT_EQUAL(imframe_ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~IMFrameView()))
{
  delete imframe_ptr;
}
END_SECTION

START_SECTION((const MSSpectrum* getSpectrum() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  TEST_EQUAL(imframe.getSpectrum(), &spectrum)
}
END_SECTION

START_SECTION((MSSpectrum* getSpectrum()))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  TEST_EQUAL(imframe.getSpectrum(), &spectrum)
}
END_SECTION

START_SECTION((std::vector<IMPeak> getPeaks() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  std::vector<IMFrameView::IMPeak> peaks = imframe.getPeaks();
  
  TEST_EQUAL(peaks.size(), 3)
  TEST_REAL_SIMILAR(peaks[0].mz, 100.0)
  TEST_REAL_SIMILAR(peaks[0].intensity, 1000.0)
  TEST_REAL_SIMILAR(peaks[0].drift_time, 10.0)
  
  TEST_REAL_SIMILAR(peaks[1].mz, 200.0)
  TEST_REAL_SIMILAR(peaks[1].intensity, 2000.0)
  TEST_REAL_SIMILAR(peaks[1].drift_time, 20.0)
  
  TEST_REAL_SIMILAR(peaks[2].mz, 300.0)
  TEST_REAL_SIMILAR(peaks[2].intensity, 3000.0)
  TEST_REAL_SIMILAR(peaks[2].drift_time, 30.0)
}
END_SECTION

START_SECTION((std::vector<double> getMZArray() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  std::vector<double> mzs = imframe.getMZArray();
  
  TEST_EQUAL(mzs.size(), 3)
  TEST_REAL_SIMILAR(mzs[0], 100.0)
  TEST_REAL_SIMILAR(mzs[1], 200.0)
  TEST_REAL_SIMILAR(mzs[2], 300.0)
}
END_SECTION

START_SECTION((std::vector<float> getIntensityArray() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  std::vector<float> intensities = imframe.getIntensityArray();
  
  TEST_EQUAL(intensities.size(), 3)
  TEST_REAL_SIMILAR(intensities[0], 1000.0)
  TEST_REAL_SIMILAR(intensities[1], 2000.0)
  TEST_REAL_SIMILAR(intensities[2], 3000.0)
}
END_SECTION

START_SECTION((std::vector<double> getDriftTimeArray() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  std::vector<double> drift_times = imframe.getDriftTimeArray();
  
  TEST_EQUAL(drift_times.size(), 3)
  TEST_REAL_SIMILAR(drift_times[0], 10.0)
  TEST_REAL_SIMILAR(drift_times[1], 20.0)
  TEST_REAL_SIMILAR(drift_times[2], 30.0)
}
END_SECTION

START_SECTION((Size size() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  TEST_EQUAL(imframe.size(), 3)
}
END_SECTION

START_SECTION((bool empty() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  TEST_EQUAL(imframe.empty(), false)
  
  MSSpectrum empty_spectrum;
  empty_spectrum.setType(SpectrumSettings::CENTROID);
  DataArrays::FloatDataArray empty_drift_times;
  IMDataConverter::setIMUnit(empty_drift_times, DriftTimeUnit::MILLISECOND);
  empty_spectrum.getFloatDataArrays().push_back(empty_drift_times);
  
  IMFrameView empty_imframe(&empty_spectrum);
  TEST_EQUAL(empty_imframe.empty(), true)
}
END_SECTION

START_SECTION((bool isValid() const))
{
  MSSpectrum valid_spectrum = createTestSpectrum();
  IMFrameView valid_imframe(&valid_spectrum); // Constructor succeeds
  TEST_EQUAL(valid_imframe.isValid(), true)
  
  MSSpectrum invalid_spectrum = createInvalidSpectrum();
  // For spectra with IMFormat::NONE, the constructor should throw an exception
  TEST_EXCEPTION(Exception::InvalidValue, IMFrameView invalid_imframe(&invalid_spectrum))
  
  // For spectra with IMFormat::MULTIPLE_SPECTRA, the constructor should throw an exception
  MSSpectrum multiple_spectra_spectrum = createMultipleSpectraFormatSpectrum();
  TEST_EXCEPTION(Exception::InvalidValue, IMFrameView multiple_spectra_imframe(&multiple_spectra_spectrum))
  
  // For non-centroided spectra, the constructor should throw an exception
  MSSpectrum non_centroided_spectrum = createNonCentroidedSpectrum();
  TEST_EXCEPTION(Exception::InvalidValue, IMFrameView non_centroided_imframe(&non_centroided_spectrum))
  
  // For spectra with mismatched drift time array size, the constructor should throw an exception
  MSSpectrum mismatched_spectrum = createMismatchedSpectrum();
  TEST_EXCEPTION(Exception::InvalidValue, IMFrameView mismatched_imframe(&mismatched_spectrum))
  
  // Note: isValid() always returns true if the constructor succeeds
  // If the spectrum is invalid, the constructor will throw an exception
}
END_SECTION

START_SECTION((const MSSpectrum::FloatDataArray* getDriftTimeDataArrayPtr() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  const MSSpectrum::FloatDataArray* dt_array = imframe.getDriftTimeDataArrayPtr();
  
  TEST_NOT_EQUAL(dt_array, nullPointer)
  TEST_EQUAL(IMTypes::determineIMFormat(spectrum), IMFormat::CONCATENATED)
  TEST_EQUAL(dt_array->size(), 3)
  TEST_REAL_SIMILAR((*dt_array)[0], 10.0)
  TEST_REAL_SIMILAR((*dt_array)[1], 20.0)
  TEST_REAL_SIMILAR((*dt_array)[2], 30.0)
}
END_SECTION

START_SECTION((double getRT() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  TEST_REAL_SIMILAR(imframe.getRT(), 42.0)
}
END_SECTION

START_SECTION((void setRT(double rt)))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  imframe.setRT(99.9);
  TEST_REAL_SIMILAR(imframe.getRT(), 99.9)
  TEST_REAL_SIMILAR(spectrum.getRT(), 99.9)
}
END_SECTION

START_SECTION((UInt getMSLevel() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  TEST_EQUAL(imframe.getMSLevel(), 1)
}
END_SECTION

START_SECTION((void setMSLevel(UInt level)))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  imframe.setMSLevel(2);
  TEST_EQUAL(imframe.getMSLevel(), 2)
  TEST_EQUAL(spectrum.getMSLevel(), 2)
}
END_SECTION

START_SECTION((const std::string& getNativeID() const))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  TEST_EQUAL(imframe.getNativeID(), "test_spectrum")
}
END_SECTION

START_SECTION((void setNativeID(const std::string& id)))
{
  MSSpectrum spectrum = createTestSpectrum();
  IMFrameView imframe(&spectrum);
  
  imframe.setNativeID("new_id");
  TEST_EQUAL(imframe.getNativeID(), "new_id")
  TEST_EQUAL(spectrum.getNativeID(), "new_id")
}
END_SECTION

START_SECTION((void sortByMZ()))
{
  MSSpectrum spectrum;
  spectrum.setType(SpectrumSettings::CENTROID);
  
  // Add peaks in unsorted order
  Peak1D peak1;
  peak1.setMZ(300.0);
  peak1.setIntensity(3000.0f);
  spectrum.push_back(peak1);
  
  Peak1D peak2;
  peak2.setMZ(100.0);
  peak2.setIntensity(1000.0f);
  spectrum.push_back(peak2);
  
  Peak1D peak3;
  peak3.setMZ(200.0);
  peak3.setIntensity(2000.0f);
  spectrum.push_back(peak3);
  
  // Add drift time array in corresponding order
  DataArrays::FloatDataArray drift_times;
  IMDataConverter::setIMUnit(drift_times, DriftTimeUnit::MILLISECOND);
  drift_times.push_back(30.0f);
  drift_times.push_back(10.0f);
  drift_times.push_back(20.0f);
  
  spectrum.getFloatDataArrays().push_back(drift_times);
  
  // Create IMFrameView and sort
  IMFrameView imframe(&spectrum);
  imframe.sortByMZ();
  
  // Check if peaks are sorted by m/z
  TEST_EQUAL(spectrum.size(), 3)
  TEST_REAL_SIMILAR(spectrum[0].getMZ(), 100.0)
  TEST_REAL_SIMILAR(spectrum[1].getMZ(), 200.0)
  TEST_REAL_SIMILAR(spectrum[2].getMZ(), 300.0)
  
  // Check if drift times are sorted accordingly
  const MSSpectrum::FloatDataArray* dt_array = imframe.getDriftTimeDataArrayPtr();
  TEST_EQUAL(dt_array->size(), 3)
  TEST_REAL_SIMILAR((*dt_array)[0], 10.0)
  TEST_REAL_SIMILAR((*dt_array)[1], 20.0)
  TEST_REAL_SIMILAR((*dt_array)[2], 30.0)
}
END_SECTION

START_SECTION((void sortByIntensity(bool reverse = true)))
{
  MSSpectrum spectrum;
  spectrum.setType(SpectrumSettings::CENTROID);
  
  // Add peaks
  Peak1D peak1;
  peak1.setMZ(100.0);
  peak1.setIntensity(1000.0f);
  spectrum.push_back(peak1);
  
  Peak1D peak2;
  peak2.setMZ(200.0);
  peak2.setIntensity(3000.0f);
  spectrum.push_back(peak2);
  
  Peak1D peak3;
  peak3.setMZ(300.0);
  peak3.setIntensity(2000.0f);
  spectrum.push_back(peak3);
  
  // Add drift time array
  DataArrays::FloatDataArray drift_times;
  IMDataConverter::setIMUnit(drift_times, DriftTimeUnit::MILLISECOND);
  drift_times.push_back(10.0f);
  drift_times.push_back(30.0f);
  drift_times.push_back(20.0f);
  
  spectrum.getFloatDataArrays().push_back(drift_times);
  
  // Create IMFrameView and sort by intensity (descending by default)
  IMFrameView imframe(&spectrum);
  imframe.sortByIntensity();
  
  // Check if peaks are sorted by intensity (descending)
  TEST_EQUAL(spectrum.size(), 3)
  TEST_REAL_SIMILAR(spectrum[0].getIntensity(), 3000.0)
  TEST_REAL_SIMILAR(spectrum[1].getIntensity(), 2000.0)
  TEST_REAL_SIMILAR(spectrum[2].getIntensity(), 1000.0)
  
  // Check if drift times are sorted accordingly
  const MSSpectrum::FloatDataArray* dt_array = imframe.getDriftTimeDataArrayPtr();
  TEST_EQUAL(dt_array->size(), 3)
  TEST_REAL_SIMILAR((*dt_array)[0], 30.0)
  TEST_REAL_SIMILAR((*dt_array)[1], 20.0)
  TEST_REAL_SIMILAR((*dt_array)[2], 10.0)
  
  // Sort by intensity ascending
  imframe.sortByIntensity(false);
  
  // Check if peaks are sorted by intensity (ascending)
  TEST_EQUAL(spectrum.size(), 3)
  TEST_REAL_SIMILAR(spectrum[0].getIntensity(), 1000.0)
  TEST_REAL_SIMILAR(spectrum[1].getIntensity(), 2000.0)
  TEST_REAL_SIMILAR(spectrum[2].getIntensity(), 3000.0)
  
  // Check if drift times are sorted accordingly
  dt_array = imframe.getDriftTimeDataArrayPtr();
  TEST_EQUAL(dt_array->size(), 3)
  TEST_REAL_SIMILAR((*dt_array)[0], 10.0)
  TEST_REAL_SIMILAR((*dt_array)[1], 20.0)
  TEST_REAL_SIMILAR((*dt_array)[2], 30.0)
}
END_SECTION

START_SECTION((void sortByDriftTime()))
{
  MSSpectrum spectrum;
  spectrum.setType(SpectrumSettings::CENTROID);
  
  // Add peaks
  Peak1D peak1;
  peak1.setMZ(100.0);
  peak1.setIntensity(1000.0f);
  spectrum.push_back(peak1);
  
  Peak1D peak2;
  peak2.setMZ(200.0);
  peak2.setIntensity(2000.0f);
  spectrum.push_back(peak2);
  
  Peak1D peak3;
  peak3.setMZ(300.0);
  peak3.setIntensity(3000.0f);
  spectrum.push_back(peak3);
  
  // Add drift time array in unsorted order
  DataArrays::FloatDataArray drift_times;
  IMDataConverter::setIMUnit(drift_times, DriftTimeUnit::MILLISECOND);
  drift_times.push_back(30.0f);
  drift_times.push_back(10.0f);
  drift_times.push_back(20.0f);
  
  spectrum.getFloatDataArrays().push_back(drift_times);
  
  // Create IMFrameView and sort by drift time
  IMFrameView imframe(&spectrum);
  imframe.sortByDriftTime();
  
  // Check if drift times are sorted
  const MSSpectrum::FloatDataArray* dt_array = imframe.getDriftTimeDataArrayPtr();
  TEST_EQUAL(dt_array->size(), 3)
  TEST_REAL_SIMILAR((*dt_array)[0], 10.0)
  TEST_REAL_SIMILAR((*dt_array)[1], 20.0)
  TEST_REAL_SIMILAR((*dt_array)[2], 30.0)
  
  // Check if peaks are sorted accordingly
  TEST_EQUAL(spectrum.size(), 3)
  TEST_REAL_SIMILAR(spectrum[0].getMZ(), 200.0)
  TEST_REAL_SIMILAR(spectrum[1].getMZ(), 300.0)
  TEST_REAL_SIMILAR(spectrum[2].getMZ(), 100.0)
}
END_SECTION

START_SECTION((bool operator==(const IMFrameView& rhs) const))
{
  MSSpectrum spectrum1 = createTestSpectrum();
  IMFrameView imframe1(&spectrum1);
  
  MSSpectrum spectrum2 = createTestSpectrum();
  IMFrameView imframe2(&spectrum2);
  
  // Different objects but same content
  TEST_EQUAL(imframe1 == imframe2, true)
  
  // Same object
  TEST_EQUAL(imframe1 == imframe1, true)
  
  // Different content
  spectrum2[0].setIntensity(999.0f);
  TEST_EQUAL(imframe1 == imframe2, false)
}
END_SECTION

START_SECTION((bool operator!=(const IMFrameView& rhs) const))
{
  MSSpectrum spectrum1 = createTestSpectrum();
  IMFrameView imframe1(&spectrum1);
  
  MSSpectrum spectrum2 = createTestSpectrum();
  IMFrameView imframe2(&spectrum2);
  
  // Different objects but same content
  TEST_EQUAL(imframe1 != imframe2, false)
  
  // Same object
  TEST_EQUAL(imframe1 != imframe1, false)
  
  // Different content
  spectrum2[0].setIntensity(999.0f);
  TEST_EQUAL(imframe1 != imframe2, true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST