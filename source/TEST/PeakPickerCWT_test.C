// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

CHECK((~PeakPickerCWT()))
  delete ptr;
RESULT

CHECK(setParam)
  Param param;
  param.setValue("thresholds:correlation",0.8);
  param.setValue("wavelet_transform:scale",0.3);
  param.setValue("thresholds:noise_level",9);
  param.setValue("thresholds:search_radius",2);
  param.setValue("deconvolution:skip_deconvolution","no");
  param.setValue("Optimization:optimization","two_dimensional");
  
  PeakPickerCWT pp;
  pp.setParam(param);
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.3)
  TEST_REAL_EQUAL(pp.getPeakCorrBound(),0.8)
  TEST_REAL_EQUAL(pp.getNoiseLevel(),9)
  TEST_EQUAL(pp.getOptimizationFlag() == false, true)
  TEST_REAL_EQUAL(pp.getSearchRadius(),2)
	TEST_EQUAL(pp.getDeconvolutionFlag() == true, true)
	TEST_EQUAL(pp.get2DOptimizationFlag() == true, true)
RESULT

CHECK((PeakPickerCWT& operator=(const PeakPickerCWT& pp)))
  Param param;
  param.setValue("thresholds:correlation",0.8);
  param.setValue("wavelet_transform:scale",0.3);
  param.setValue("thresholds:noise_level",9);
  param.setValue("thresholds:search_radius",2);
  param.setValue("deconvolution:skip_deconvolution","no");
  param.setValue("Optimization:optimization","one_dimensional");
  
  PeakPickerCWT pp;
  pp.setParam(param);
  PeakPickerCWT pp_copy;
  pp_copy = pp;
  TEST_REAL_EQUAL(pp_copy.getWaveletScale(),0.3)
  TEST_REAL_EQUAL(pp_copy.getPeakCorrBound(),0.8)
  TEST_REAL_EQUAL(pp_copy.getNoiseLevel(),9)
  TEST_EQUAL(pp_copy.getOptimizationFlag() == true, true)
  TEST_REAL_EQUAL(pp_copy.getSearchRadius(),2)
	TEST_EQUAL(pp.getDeconvolutionFlag() == true, true)
	TEST_EQUAL(pp.get2DOptimizationFlag() == false, true)
RESULT

CHECK((PeakPickerCWT(const PeakPickerCWT& pp)))
  Param param;
  param.setValue("thresholds:correlation",0.8);
  param.setValue("wavelet_transform:scale",0.3);
  param.setValue("thresholds:noise_level",9);
  param.setValue("thresholds:search_radius",2);
  param.setValue("deconvolution:skip_deconvolution","no");
  param.setValue("Optimization:optimization","no");
  
  PeakPickerCWT pp;
  pp.setParam(param);
  
  PeakPickerCWT pp_copy(pp);
  TEST_REAL_EQUAL(pp_copy.getWaveletScale(),0.3)
  TEST_REAL_EQUAL(pp_copy.getPeakCorrBound(),0.8)
  TEST_REAL_EQUAL(pp_copy.getNoiseLevel(),9)
  TEST_EQUAL(pp_copy.getOptimizationFlag() == false, true)
  TEST_REAL_EQUAL(pp_copy.getSearchRadius(),2)
	TEST_EQUAL(pp.getDeconvolutionFlag() == true, true)
	TEST_EQUAL(pp.get2DOptimizationFlag() == false, true)	
RESULT

MzDataFile mz_data_file;
MSExperiment<DRawDataPoint<1> > exp_raw;
mz_data_file.load("data/PeakPicker_test.mzData",exp_raw);
CHECK((template<typename InputPeakContainer, typename OutputPeakContainer > void pick(const InputPeakContainer& input_peak_container, OutputPeakContainer& picked_peaks_container, int ms_level = 1)))
  MSSpectrum<DPeak<1> > peaks;
  PeakPickerCWT pp;
    
  pp.pick(exp_raw[0],peaks);
  MSSpectrum<DPeak<1> >::const_iterator it = peaks.begin();
  TEST_REAL_EQUAL(peaks.size() == pp.getPeakShapes().size(), true)  
  TEST_REAL_EQUAL(it->getPos(),pp.getPeakShapes()[0].mz_position)
  TEST_REAL_EQUAL(it->getIntensity(),pp.getPeakShapes()[0].height)
RESULT

CHECK((template<typename InputPeakIterator, typename OutputPeakContainer  > void pick(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& picked_peak_container, int ms_level = 1)))
  MSSpectrum<DPickedPeak<1> > peaks;
  PeakPickerCWT pp;
  
  pp.pick(exp_raw[0].begin(),exp_raw[0].end(),peaks,1);
  MSSpectrum<DPickedPeak<1> >::const_iterator it = peaks.begin();
  TEST_REAL_EQUAL(peaks.size() == pp.getPeakShapes().size(), true)   
  TEST_REAL_EQUAL(it->getPos(),pp.getPeakShapes()[0].mz_position)
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
  MSExperiment<DPickedPeak<1> > peaks;
  PeakPickerCWT pp;
  pp.setPeakBound(1500);
   
  pp.pickExperiment(exp_raw,peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
  ExperimentalSettings e = peaks;
  TEST_EQUAL(e == exp_raw, true)
RESULT

CHECK((template<typename InputSpectrumIterator, typename OutputPeakType > void pickExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_peaks)))
  MSExperiment<DPickedPeak<1> > peaks;
  PeakPickerCWT pp;
  pp.setPeakBound(1500);
   
  pp.pickExperiment(exp_raw.begin(),exp_raw.end(),peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)   
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
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

	
CHECK((const float& getNoiseLevel() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getNoiseLevel(),0.1)
RESULT

CHECK((const float& getPeakBoundCWT() const))
  PeakPickerCWT pp;
  TEST_REAL_EQUAL(pp.getPeakBoundCWT(),0.0)
RESULT

CHECK((const float& getPeakBoundMs2LevelCWT() const))
  PeakPickerCWT pp;
  TEST_REAL_EQUAL(pp.getPeakBoundMs2LevelCWT(),0.0)
RESULT

CHECK((const float& getPeakCorrBound() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getPeakCorrBound(),0.5)
RESULT

CHECK((const float& getWaveletScale() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.15)
RESULT

CHECK((const unsigned int& getSearchRadius() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getSearchRadius(),3)
RESULT

CHECK((const std::vector<PeakShape>& getPeakShapes() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getPeakShapes().size(),0)
RESULT

CHECK((void setNoiseLevel(const float& noise_level)))
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

CHECK((void set2DOptimizationFlag(const bool& twod_optimization)))
  PeakPickerCWT pp;
  
  pp.set2DOptimizationFlag(true);
  TEST_REAL_EQUAL(pp.get2DOptimizationFlag(),true)
RESULT
	
	
CHECK((void setWaveletScale(const float& scale)))
  PeakPickerCWT pp;
  
  pp.setWaveletScale(0.1);
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



