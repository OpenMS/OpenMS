// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// $Authors: Stephan Aiche, Chris Bielow, Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/RangeUtils.h>

///////////////////////////
#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RawTandemMSSignalSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RawTandemMSSignalSimulation* ptr = 0;
RawTandemMSSignalSimulation* null_ptr = 0;
SimRandomNumberGenerator rng;

START_SECTION((RawTandemMSSignalSimulation(const SimRandomNumberGenerator &rng)))
{
	ptr = new RawTandemMSSignalSimulation(rng);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~RawTandemMSSignalSimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((RawTandemMSSignalSimulation(const RawTandemMSSignalSimulation &source)))
{
    ptr = new RawTandemMSSignalSimulation(rng);
    Param tmp_par = ptr->getParameters();
    tmp_par.setValue("status", "precursor");
    ptr->setParameters(tmp_par);

    RawTandemMSSignalSimulation copy(*ptr);
    TEST_EQUAL(copy.getParameters(), ptr->getParameters())
}
END_SECTION

START_SECTION((RawTandemMSSignalSimulation& operator=(const RawTandemMSSignalSimulation &source)))
{
    RawTandemMSSignalSimulation copy(rng);
    copy = *ptr;
    TEST_EQUAL(copy.getParameters(), ptr->getParameters())
}
END_SECTION

START_SECTION((void generateRawTandemSignals(const FeatureMapSim &, MSSimExperiment &)))
{
    rng.biological_rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng.biological_rng, 0);

    //Load featureXML and MSExperiment from MSSimulator run without MS2 simulation
    String feature_filename, exp_no_ms2_file, exp_with_ms2_file;
    feature_filename = OPENMS_GET_TEST_DATA_PATH("RawTandemMSSignalSimulation_no_ms2.featureXML");
    exp_no_ms2_file = OPENMS_GET_TEST_DATA_PATH("RawTandemMSSignalSimulation_no_ms2.mzML");
    exp_with_ms2_file = OPENMS_GET_TEST_DATA_PATH("RawTandemMSSignalSimulation_with_ms2.mzML");;
    FeatureMapSim features;
    MSSimExperiment exp_no_ms2, exp_with_ms2;
    FeatureXMLFile().load(feature_filename, features);
    MzMLFile().load(exp_no_ms2_file, exp_no_ms2);
    MzMLFile().load(exp_with_ms2_file, exp_with_ms2);

    RawTandemMSSignalSimulation sim(rng);
    Param p;
    p.setValue("status", "precursor");
    p.setValue("tandem_mode", 2);
    p.setValue("TandemSim:SVM:hide_losses", "true");
    p.setValue("Precursor:Exclusion:use_dynamic_exclusion", "true");
    p.setValue("Precursor:Exclusion:exclusion_time", 50.0);
    sim.setParameters(p);

    sim.generateRawTandemSignals(features, exp_no_ms2);
    IntList levels;
    levels.push_back(1);
    exp_no_ms2.erase(remove_if(exp_no_ms2.begin(), exp_no_ms2.end(), InMSLevelRange<MSSimExperiment::SpectrumType>(levels)), exp_no_ms2.end());
    exp_with_ms2.erase(remove_if(exp_with_ms2.begin(), exp_with_ms2.end(), InMSLevelRange<MSSimExperiment::SpectrumType>(levels)), exp_with_ms2.end());

    TEST_EQUAL(exp_with_ms2.size(), exp_no_ms2.size());
    TEST_EQUAL(exp_with_ms2[0].size(), exp_no_ms2[0].size());
    TEST_EQUAL(exp_with_ms2[1].size(), exp_no_ms2[1].size());
    TEST_REAL_SIMILAR(exp_with_ms2[0].getPrecursors()[0].getMZ(), exp_no_ms2[0].getPrecursors()[0].getMZ());
    TEST_REAL_SIMILAR(exp_with_ms2[1].getPrecursors()[0].getMZ(), exp_no_ms2[1].getPrecursors()[0].getMZ());
    TEST_REAL_SIMILAR(exp_with_ms2[0][0].getIntensity(), exp_no_ms2[0][0].getIntensity());
    TEST_REAL_SIMILAR(exp_with_ms2[0][1].getIntensity(), exp_no_ms2[0][1].getIntensity());
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



