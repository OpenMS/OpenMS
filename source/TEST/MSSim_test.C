// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Chris Bielow, Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/MSSim.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSSim, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSim* ptr = 0;
START_SECTION(MSSim())
{
	ptr = new MSSim();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~MSSim())
{
	delete ptr;
}
END_SECTION

START_SECTION((void simulate(gsl_rng * const rnd_gen, const SampleProteins &peptides)))
{
  // TODO
#if 0 // core from old LCMSSim_test

	LCMSSample smp;
	smp.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSim_test.fasta"));
	smp.setPdModelFile(OPENMS_GET_TEST_DATA_PATH("LCMSSim_test_pd.svm"));
	smp.digest();

	Param p;
	p.setValue("ion_model",0);
	p.setValue("total_gradient_time",20.0f);
	p.setValue("random_seed",42);

	LCMSSim sim;
	sim.setParameters(p);
	sim.setRTModelFile(OPENMS_GET_TEST_DATA_PATH("LCMSSim_test.svm"));
	sim.setSample(smp);
	sim.run();	

	std::string tmp_mzdata;
	NEW_TMP_FILE(tmp_mzdata);
	sim.exportMzData(tmp_mzdata);
	
	//test if stored data is equal to created data
	MSExperiment<> exp, exp_orig;
	MzDataFile().load(tmp_mzdata, exp);
	MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("LCMSSim_test_out.mzData"), exp_orig);
	
	TEST_EQUAL(exp.size(), exp_orig.size())
	for (Size s=0; s<exp.size(); ++s)
	{
		TEST_EQUAL(exp[s].size(), exp_orig[s].size())
		for (Size p=0; p<exp[s].size(); ++p)
		{
			TEST_REAL_SIMILAR(exp[s][p].getMZ(), exp_orig[s][p].getMZ())
			TEST_REAL_SIMILAR(exp[s][p].getIntensity(), exp_orig[s][p].getIntensity())
		}
	}

#endif  
}
END_SECTION

START_SECTION((MSSimExperiment const& getExperiment() const ))
{
  MSSimExperiment empty_experiment;
  MSSim mssim;

  TEST_EQUAL(mssim.getExperiment().getSize(), empty_experiment.getSize())

  // TODO we need some more sophisticated testing here
}
END_SECTION

START_SECTION((FeatureMapSim const& getSimulatedFeatures() const ))
{
  // TODO
}
END_SECTION

START_SECTION((ConsensusMap const & getSimulatedConsensus() const))
{
  // TODO
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



