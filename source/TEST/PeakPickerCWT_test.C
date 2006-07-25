// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

CHECK((PeakPickerCWT(const String& filename)))
  Param param;
  param.setValue("Thresholds:Correlation",0.8);
  param.setValue("Optimization:SkipOptimization","yes");
  param.setValue("WaveletTransform:Scale",0.3);
  param.setValue("Thresholds:NoiseLevel",9);
  param.setValue("Thresholds:SearchRadius",2);
  
  String file("bla.xml");
  NEW_TMP_FILE(file)
  param.store(file);
    
  PeakPickerCWT pp(file);
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.3)
  TEST_REAL_EQUAL(pp.getPeakCorrBound(),0.8)
  TEST_REAL_EQUAL(pp.getNoiseLevel(),9)
  TEST_EQUAL(pp.getOptimizationValue() == false, true)
  TEST_REAL_EQUAL(pp.getSearchRadius(),2)
  std::cout << "(PeakPickerCWT(const String& filename)" << std::endl;
RESULT

CHECK((PeakPickerCWT(const Param& parameters)))
  Param param;
  param.setValue("Thresholds:Correlation",0.8);
  param.setValue("Optimization:SkipOptimization","yes");
  param.setValue("WaveletTransform:Scale",0.3);
  param.setValue("Thresholds:NoiseLevel",9);
  param.setValue("Thresholds:SearchRadius",2);
  
  PeakPickerCWT pp(param);
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.3)
  TEST_REAL_EQUAL(pp.getPeakCorrBound(),0.8)
  TEST_REAL_EQUAL(pp.getNoiseLevel(),9)
  TEST_EQUAL(pp.getOptimizationValue() == false, true)
  TEST_REAL_EQUAL(pp.getSearchRadius(),2)
  TEST_EQUAL(pp.getParam() == param, true)
  std::cout << "(PeakPickerCWT(const param)" << std::endl;
RESULT

CHECK((ContinuousWaveletTransform& getWaveletTransform()))
  ContinuousWaveletTransformNumIntegration cwt;
  cwt.setSpacing(0.1);
  
  PeakPickerCWT pp;
  pp.getWaveletTransform() = cwt;
  
  TEST_REAL_EQUAL(cwt.getSpacing(),pp.getWaveletTransform().getSpacing())
  std::cout << "getWaveletScale" << std::endl;
RESULT

CHECK((PeakPickerCWT& operator=(const PeakPickerCWT& pp)))
  Param param;
  param.setValue("Thresholds:Correlation",0.8);
  param.setValue("Optimization:SkipOptimization","yes");
  param.setValue("WaveletTransform:Scale",0.3);
  param.setValue("Thresholds:NoiseLevel",9);
  param.setValue("Thresholds:SearchRadius",2);
  
  PeakPickerCWT pp(param);
  PeakPickerCWT pp_copy;
  pp_copy = pp;
  TEST_REAL_EQUAL(pp_copy.getWaveletScale(),0.3)
  TEST_REAL_EQUAL(pp_copy.getPeakCorrBound(),0.8)
  TEST_REAL_EQUAL(pp_copy.getNoiseLevel(),9)
  TEST_EQUAL(pp_copy.getOptimizationValue() == false, true)
  TEST_REAL_EQUAL(pp_copy.getSearchRadius(),2)
  TEST_EQUAL(pp_copy.getParam() == param, true)
  std::cout << "operator=" << std::endl;
RESULT

CHECK((PeakPickerCWT(const PeakPickerCWT& pp)))
  Param param;
  param.setValue("Thresholds:Correlation",0.8);
  param.setValue("Optimization:SkipOptimization","yes");
  param.setValue("WaveletTransform:Scale",0.3);
  param.setValue("Thresholds:NoiseLevel",9);
  param.setValue("Thresholds:SearchRadius",2);
  
  PeakPickerCWT pp(param);
  PeakPickerCWT pp_copy(pp);
  TEST_REAL_EQUAL(pp_copy.getWaveletScale(),0.3)
  TEST_REAL_EQUAL(pp_copy.getPeakCorrBound(),0.8)
  TEST_REAL_EQUAL(pp_copy.getNoiseLevel(),9)
  TEST_EQUAL(pp_copy.getOptimizationValue() == false, true)
  TEST_REAL_EQUAL(pp_copy.getSearchRadius(),2)
  TEST_EQUAL(pp_copy.getParam() == param, true)
  std::cout << "copy constr" << std::endl;
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
  std::cout << "pick(it,it,cont)" << std::endl;
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
  std::cout << "pick(cont, cont)" << std::endl;
RESULT

CHECK((template<typename InputPeakType, typename OutputPeakType > void pickExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_peaks)))
  MSExperiment<DPickedPeak<1> > peaks;
  PeakPickerCWT pp;
   
  pp.pickExperiment(exp_raw,peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 277)
  ExperimentalSettings e = peaks;
  TEST_EQUAL(e == exp_raw, true)
  std::cout << "pickExp(it,it,cont) " << std::endl;
RESULT

CHECK((template<typename InputSpectrumIterator, typename OutputPeakType > void pickExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_peaks)))
  MSExperiment<DPickedPeak<1> > peaks;
  PeakPickerCWT pp;
   
  pp.pickExperiment(exp_raw.begin(),exp_raw.end(),peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)   
  TEST_EQUAL((peaks[0].size() + peaks[1].size()),277)
  std::cout << "pickExperiment(cont,cont)" << std::endl;
RESULT

CHECK((const ContinuousWaveletTransform<1>& getWaveletTransform() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getWaveletTransform().getSpacing(), 0)
RESULT

CHECK((const bool& getOptimizationValue() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getOptimizationValue(),false)
RESULT

CHECK((const float& getNoiseLevel() const))
  PeakPickerCWT pp;
  
  TEST_REAL_EQUAL(pp.getNoiseLevel(),10)
RESULT

CHECK((const float& getPeakBoundCWT() const))
  PeakPickerCWT pp;
  PRECISION(0.001)
  
  TEST_REAL_EQUAL(pp.getPeakBoundCWT(),39.546)
RESULT

CHECK((const float& getPeakBoundMs2LevelCWT() const))
  PeakPickerCWT pp;
  PRECISION(0.001)
  
  TEST_REAL_EQUAL(pp.getPeakBoundMs2LevelCWT(),9.886)
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

CHECK((float& getPeakBoundCWT()))
  PeakPickerCWT pp;
  
  pp.getPeakBoundCWT() = 41;  
  TEST_REAL_EQUAL(pp.getPeakBoundCWT(),41)
RESULT


CHECK((float& getNoiseLevel()))
  PeakPickerCWT pp;
  
  pp.getNoiseLevel() = 12;
  TEST_REAL_EQUAL(pp.getNoiseLevel(),12)
RESULT

CHECK((bool& getOptimizationValue()))
  PeakPickerCWT pp;
  
  pp.getOptimizationValue() = true;
  TEST_REAL_EQUAL(pp.getOptimizationValue(),true)
RESULT

CHECK((float& getPeakBoundMs2LevelCWT()))
  PeakPickerCWT pp;
  
  pp.getPeakBoundMs2LevelCWT() = 13.1;
  TEST_REAL_EQUAL(pp.getPeakBoundMs2LevelCWT(),13.1)
RESULT

CHECK((float& getPeakCorrBound()))
  PeakPickerCWT pp;
  
  pp.getPeakCorrBound() = 1.;
  TEST_REAL_EQUAL(pp.getPeakCorrBound(),1.)
RESULT

CHECK((float& getWaveletScale()))
  PeakPickerCWT pp;
  
  pp.getWaveletScale() = 0.3;
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.3)
RESULT

CHECK((unsigned int& getSearchRadius()))
  PeakPickerCWT pp;
  
  pp.getSearchRadius() = 1;
  TEST_REAL_EQUAL(pp.getSearchRadius(),1)
RESULT

CHECK((std::vector<PeakShape>& getPeakShapes()))
  PeakPickerCWT pp;
  
  std::vector<PeakShape> shapes;
  pp.getPeakShapes() = shapes;
  TEST_EQUAL(shapes.size() == pp.getPeakShapes().size(), true)  
RESULT

CHECK((void setNoiseLevel(const float& noise_level)))
  PeakPickerCWT pp;
  
  pp.setNoiseLevel(12);
  TEST_REAL_EQUAL(pp.getNoiseLevel(),12)
RESULT

CHECK((void setOptimizationValue(const bool& optimization)))
  PeakPickerCWT pp;
  
  pp.setOptimizationValue(true);
  TEST_REAL_EQUAL(pp.getOptimizationValue(),true)
RESULT

CHECK((void setPeakBoundCWT(const float peak_bound_cwt)))
  PeakPickerCWT pp;
  
  pp.setPeakBoundCWT(41);  
  TEST_REAL_EQUAL(pp.getPeakBoundCWT(),41)
RESULT

CHECK((void setPeakBoundMs2LevelCWT(const float& peak_bound_ms2_level_cwt)))
   PeakPickerCWT pp;
   
    pp.setPeakBoundMs2LevelCWT(12);
   TEST_REAL_EQUAL(pp.getPeakBoundMs2LevelCWT(),12)
RESULT

CHECK((void setPeakCorrBound(const float& peak_corr_bound)))
  PeakPickerCWT pp;
  
  pp.setPeakCorrBound(1.);
  TEST_REAL_EQUAL(pp.getPeakCorrBound(),1.)
RESULT

CHECK((void setPeakShapes(const std::vector<PeakShape>& peak_shapes)))
  PeakPickerCWT pp;
  
  std::vector<PeakShape> shapes;
  pp.setPeakShapes(shapes);
  TEST_EQUAL(shapes.size() == pp.getPeakShapes().size(), true) 
RESULT

CHECK((void setSearchRadius(const unsigned int& radius)))
  PeakPickerCWT pp;
  
  pp.setSearchRadius(1);
  TEST_REAL_EQUAL(pp.getSearchRadius(),1)
RESULT

CHECK((void setWaveletScale(const float& scale)))
  PeakPickerCWT pp;
  
  pp.setWaveletScale(0.1);
  TEST_REAL_EQUAL(pp.getWaveletScale(),0.1)
RESULT

CHECK((void setWaveletTransform(const ContinuousWaveletTransform& wt)))
  ContinuousWaveletTransformNumIntegration cwt;
  cwt.setSpacing(0.1);
  
  PeakPickerCWT pp;
  pp.setWaveletTransform(cwt);
  
  TEST_REAL_EQUAL(cwt.getSpacing(),pp.getWaveletTransform().getSpacing())
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



