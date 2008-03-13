// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerCWT, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerCWT* ptr = 0;
CHECK((PeakPickerCWT()))
  ptr = new PeakPickerCWT();
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~PeakPickerCWT()))
  delete ptr;
RESULT


MzDataFile mz_data_file;
MSExperiment<RawDataPoint1D > exp_raw;
mz_data_file.load("data/PeakPicker_test.mzData",exp_raw);
CHECK((template<typename InputPeakContainer, typename OutputPeakContainer > void pick(const InputPeakContainer& input_peak_container, OutputPeakContainer& picked_peaks_container, int ms_level = 1)))
  MSSpectrum<Peak1D > peaks;
  PeakPickerCWT pp;
    
  pp.pick(exp_raw[0],peaks);
  MSSpectrum<Peak1D >::const_iterator it = peaks.begin();
  TEST_REAL_EQUAL(peaks.size() == pp.getPeakShapes().size(), true)  
  TEST_REAL_EQUAL(it->getMZ(),pp.getPeakShapes()[0].mz_position)
  TEST_REAL_EQUAL(it->getIntensity(),pp.getPeakShapes()[0].height)
RESULT

CHECK((template<typename InputPeakIterator, typename OutputPeakContainer  > void pick(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& picked_peak_container, int ms_level = 1)))
  MSSpectrum<PickedPeak1D > peaks;
  PeakPickerCWT pp;
  
  pp.pick(exp_raw[0].begin(),exp_raw[0].end(),peaks,1);
  MSSpectrum<PickedPeak1D >::const_iterator it = peaks.begin();
  TEST_REAL_EQUAL(peaks.size() == pp.getPeakShapes().size(), true)   
  TEST_REAL_EQUAL(it->getMZ(),pp.getPeakShapes()[0].mz_position)
  TEST_REAL_EQUAL(it->getIntensity(),pp.getPeakShapes()[0].height)
  TEST_REAL_EQUAL(it->getRValue(),pp.getPeakShapes()[0].r_value)
  TEST_REAL_EQUAL(it->getArea(),pp.getPeakShapes()[0].area)
  TEST_REAL_EQUAL(it->getFWHM(),pp.getPeakShapes()[0].getFWHM())
  TEST_REAL_EQUAL(it->getLeftWidthParameter(),pp.getPeakShapes()[0].left_width)
  TEST_REAL_EQUAL(it->getRightWidthParameter(),pp.getPeakShapes()[0].right_width)
  TEST_REAL_EQUAL(it->getPeakShape(),pp.getPeakShapes()[0].type)
  TEST_REAL_EQUAL(it->getSN(),pp.getPeakShapes()[0].signal_to_noise)
RESULT

CHECK((template<typename InputPeakType, typename OutputPeakType > void pickExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_peaks)))
  MSExperiment<PickedPeak1D > peaks;
  PeakPickerCWT pp;
  pp.setPeakBound(1500);
   
  pp.pickExperiment(exp_raw,peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
  ExperimentalSettings e = peaks;
  TEST_EQUAL(e == exp_raw, true)
RESULT

CHECK((template <typename InputSpectrumIterator, typename OutputPeakType> void pickExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperimentExtern< OutputPeakType > &ms_exp_peaks)))
  MSExperimentExtern<PickedPeak1D > peaks;
  PeakPickerCWT pp;
  pp.setPeakBound(1500);
   
  pp.pickExperiment(exp_raw.begin(),exp_raw.end(),peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)   
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
RESULT
		
CHECK((template<typename InputSpectrumIterator, typename OutputPeakType > void pickExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_peaks)))
  MSExperiment<PickedPeak1D > peaks;
  PeakPickerCWT pp;
  pp.setPeakBound(1500);
   
  pp.pickExperiment(exp_raw.begin(),exp_raw.end(),peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)   
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
RESULT

CHECK((template <typename InputPeakType, typename OutputPeakType> void pickExperiment(const MSExperimentExtern< InputPeakType > &ms_exp_raw, MSExperimentExtern< OutputPeakType > &ms_exp_peaks)))
	MSExperimentExtern<RawDataPoint1D > exp_raw_ext;
  mz_data_file.load("data/PeakPicker_test.mzData",exp_raw_ext);
	MSExperimentExtern<PickedPeak1D > peaks;
  PeakPickerCWT pp;
  pp.setPeakBound(1500);
   
  pp.pickExperiment(exp_raw_ext.begin(),exp_raw_ext.end(),peaks);
  TEST_EQUAL(peaks.size() == exp_raw_ext.size(), true)   
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
RESULT

CHECK((template <typename OutputPeakType> void fillPeak(const PeakShape &, OutputPeakType &)))

RESULT

CHECK((template <typename OutputPeakType> void fillPeak(const PeakShape &peak_shape, PickedPeak1D &picked_peak)))
  double height = 100.0;
  double mz_position = 0.0;
  double left_width = 4.0;
  double right_width = 4.0;
  double area = 100;
  PeakShapeType::Enum type = PeakShapeType::LORENTZ_PEAK;
    
  PeakShape p(height,
							mz_position,
							left_width,
							right_width,
							area,
							std::vector<RawDataPoint1D>::iterator(),
							std::vector<RawDataPoint1D>::iterator(),
							type);
  PeakPickerCWT pp;
  PickedPeak1D peak;
  pp.fillPeak(p,peak);
  TEST_REAL_EQUAL(p.r_value,peak.getRValue())
  TEST_REAL_EQUAL(p.getFWHM(),peak.getFWHM())
	TEST_REAL_EQUAL(p.left_width,peak.getLeftWidthParameter())
	TEST_REAL_EQUAL(p.right_width,peak.getRightWidthParameter())
	TEST_REAL_EQUAL(p.area,peak.getArea())
	TEST_EQUAL(p.type,peak.getPeakShape())
 	TEST_EQUAL(p.signal_to_noise,peak.getSN())
RESULT
	
CHECK((const ContinuousWaveletTransformNumIntegration& getWaveletTransform() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getWaveletTransform().getSpacing(), 0.0)
RESULT

CHECK((const bool& getOptimizationFlag() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getOptimizationFlag(),false)
RESULT

CHECK((const bool& getDeconvolutionFlag() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getDeconvolutionFlag(),false)
RESULT

CHECK((const bool& get2DOptimizationFlag() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.get2DOptimizationFlag(),false)
RESULT

	
CHECK((Real getNoiseLevel() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getNoiseLevel(),0.1)
RESULT

CHECK((Real getPeakBoundCWT() const))
  PeakPickerCWT pp;
  TEST_REAL_EQUAL(pp.getPeakBoundCWT(),0.0)
RESULT
	
CHECK((Real getPeakBound() const))
  PeakPickerCWT pp;
  pp.setPeakBound(320);
  TEST_REAL_EQUAL(pp.getPeakBound(),320)
RESULT
	
CHECK((Real getPeakBoundMs2LevelCWT() const))
  PeakPickerCWT pp;
  TEST_REAL_EQUAL(pp.getPeakBoundMs2LevelCWT(),0.0)
RESULT

CHECK((Real getPeakCorrBound() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getPeakCorrBound(),0.5)
RESULT

CHECK((void setPeakCorrBound(Real peak_corr_bound)))
  PeakPickerCWT pp;
  pp.setPeakCorrBound(0.6);
  TEST_REAL_EQUAL(pp.getPeakCorrBound(),0.6)
RESULT
	
CHECK((Real getWaveletScale() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.15)
RESULT

CHECK((UInt getSearchRadius() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getSearchRadius(),3)
RESULT

CHECK((const std::vector<PeakShape>& getPeakShapes() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getPeakShapes().size(),0)
RESULT

CHECK((void setNoiseLevel(Real noise_level)))
  PeakPickerCWT pp;
  
  pp.setNoiseLevel(12);
  TEST_REAL_EQUAL(pp.getNoiseLevel(),12)
RESULT

CHECK((void setOptimizationFlag(const bool& optimization)))
  PeakPickerCWT pp;
  
  pp.setOptimizationFlag(true);
  TEST_REAL_EQUAL(pp.getOptimizationFlag(),true)
RESULT

CHECK((void setDeconvolutionFlag(const bool& deconvolution)))
  PeakPickerCWT pp;
  
  pp.setDeconvolutionFlag(true);
  TEST_REAL_EQUAL(pp.getDeconvolutionFlag(),true)
RESULT

CHECK((void set2DOptimizationFlag(const bool &two_d_optimization)))
  PeakPickerCWT pp;
  
  pp.set2DOptimizationFlag(true);
  TEST_REAL_EQUAL(pp.get2DOptimizationFlag(),true)
RESULT
	
	
CHECK((void setWaveletScale(Real scale)))
  PeakPickerCWT pp;
  
  pp.setWaveletScale(0.1);
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.1)
RESULT

CHECK((void setSearchRadius(UInt radius)))
  PeakPickerCWT pp;
  pp.setSearchRadius(5);
  TEST_REAL_EQUAL(pp.getSearchRadius(),5)
RESULT		
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



