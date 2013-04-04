// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
MSSim* nullPointer = 0;
START_SECTION(MSSim())
{
	ptr = new MSSim();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MSSim())
{
	delete ptr;
}
END_SECTION

START_SECTION((void simulate(const SimRandomNumberGenerator &rnd_gen, SampleChannels &peptides)))
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

START_SECTION((ConsensusMap& getChargeConsensus() ))
{
  // TODO
}
END_SECTION

START_SECTION((ConsensusMap& getLabelingConsensus() ))
{
  // TODO
}
END_SECTION

START_SECTION((FeatureMapSim const& getContaminants() const ))
{
  // TODO
}
END_SECTION

START_SECTION((Param getParameters() const ))
{
  // TODO
}
END_SECTION

START_SECTION((MSSimExperiment const& getPeakMap() const ))
{
  // TODO
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



