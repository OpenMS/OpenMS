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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/ANALYSIS/ID/PILISNeutralLossModel.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
///////////////////////////

#include <algorithm>

using namespace OpenMS;
using namespace std;

START_TEST(PILISNeutralLossModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RichPeakSpectrum spec1, spec2, spec3;

TheoreticalSpectrumGenerator tsg;
Param tsg_param(tsg.getParameters());
tsg_param.setValue("add_metainfo", "true");
tsg_param.setValue("add_losses", "true");
tsg.setParameters(tsg_param);
tsg.getSpectrum(spec1, AASequence("DFPIANGER"), 1);
tsg.getSpectrum(spec2, AASequence("DFPIANGEK"), 1);
tsg.getSpectrum(spec3, AASequence("DFPIANGEREK"), 1);

PILISNeutralLossModel* ptr = 0;
PILISNeutralLossModel* nullPointer = 0;
START_SECTION(PILISNeutralLossModel())
{
	ptr = new PILISNeutralLossModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~PILISNeutralLossModel())
{
	delete ptr;
}
END_SECTION

START_SECTION((PILISNeutralLossModel(const PILISNeutralLossModel &model)))
{
  PILISNeutralLossModel model1;
	Param p(model1.getParameters());
	p.setValue("ion_name", "y");
	model1.setParameters(p);

	PILISNeutralLossModel model2(model1);
	TEST_EQUAL(model1.getParameters() == model2.getParameters(), true)

  HiddenMarkovModel hmm1;
  hmm1.setPseudoCounts(13);
  PILISNeutralLossModel model;
  TEST_REAL_SIMILAR(model.getHMM().getPseudoCounts(), HiddenMarkovModel().getPseudoCounts())
  model.setHMM(hmm1);
  TEST_REAL_SIMILAR(model.getHMM().getPseudoCounts(), 13)

  PILISNeutralLossModel model3(model);
  TEST_REAL_SIMILAR(model3.getHMM().getPseudoCounts(), 13)
}
END_SECTION

START_SECTION((DoubleReal train(const RichPeakSpectrum & spec, const AASequence &peptide, DoubleReal ion_weight, UInt charge, DoubleReal peptide_weight)))
{
  PILISNeutralLossModel model;
	Param p(model.getParameters());
	p.setValue("ion_name", "y");
	model.setParameters(p);

	model.generateModel();

	for (RichPeakSpectrum::ConstIterator it = spec1.begin(); it != spec1.end(); ++it)
	{
		String ion_name = (String)it->getMetaValue("IonName");
		UInt charge = (UInt)count(ion_name.begin(), ion_name.end(), '+');
		if (ion_name.hasSubstring("y"))
		{
			ion_name.remove('+');
			ion_name.remove('y');
			AASequence suffix = AASequence("DFPIANGER").getSuffix(ion_name.toInt());
      model.train(spec1, suffix, suffix.getMonoWeight(Residue::YIon), charge, AASequence("DFPIANGER").getMonoWeight());
		}
	}

  for (RichPeakSpectrum::ConstIterator it = spec2.begin(); it != spec2.end(); ++it)
  {
    String ion_name = (String)it->getMetaValue("IonName");
    UInt charge = (UInt)count(ion_name.begin(), ion_name.end(), '+');
    if (ion_name.hasSubstring("y"))
    {
			ion_name.remove('+');
			ion_name.remove('y');
			AASequence suffix = AASequence("DFPIANGEK").getSuffix(ion_name.toInt());
      model.train(spec1, suffix, suffix.getMonoWeight(Residue::YIon), charge, AASequence("DFPIANGEK").getMonoWeight());
    }
  }

  for (RichPeakSpectrum::ConstIterator it = spec2.begin(); it != spec2.end(); ++it)
  {
    String ion_name = (String)it->getMetaValue("IonName");
    UInt charge = (UInt)count(ion_name.begin(), ion_name.end(), '+');
    if (ion_name.hasSubstring("y"))
    {
			ion_name.remove('+');
			ion_name.remove('y');
			AASequence suffix = AASequence("DFPIANGEREK").getSuffix(ion_name.toInt());
      model.train(spec1, suffix, suffix.getMonoWeight(Residue::YIon), charge, AASequence("DFPIANGEREK").getMonoWeight());
    }
  }

	vector<RichPeak1D> peaks1;
	model.getIons(peaks1, AASequence("ANGER"), 1.0);
	TEST_EQUAL(peaks1.size(), 9)

/*
	for (Size i = 0; i != peaks1.size(); ++i)
	{
		cout << peaks1[i].getMZ() << " " << peaks1[i].getIntensity() << " " << peaks1[i].getMetaValue("IonName") << endl;
	}
*/

	model.evaluate();

	vector<RichPeak1D> peaks2;
	model.getIons(peaks2, AASequence("ANGER"), 1.0);
	TEST_EQUAL(peaks2.size(), 9)

/*
	for (Size i = 0; i != peaks2.size(); ++i)
	{
		cout << peaks2[i].getMZ() << " " << peaks2[i].getIntensity() << " " << peaks2[i].getMetaValue("IonName") << endl;
	}
*/

	TEST_NOT_EQUAL(peaks1 == peaks2, true)
	
}
END_SECTION

START_SECTION((void getIons(std::vector< RichPeak1D > &peaks, const AASequence &peptide, DoubleReal initial_prob)))
{
  NOT_TESTABLE // implicitely tested above
}
END_SECTION

START_SECTION((void setHMM(const HiddenMarkovModel &model)))
{
  HiddenMarkovModel hmm1;
	hmm1.setPseudoCounts(13);
	PILISNeutralLossModel model;
	TEST_REAL_SIMILAR(model.getHMM().getPseudoCounts(), HiddenMarkovModel().getPseudoCounts())
	model.setHMM(hmm1);
	TEST_REAL_SIMILAR(model.getHMM().getPseudoCounts(), 13)
}
END_SECTION

START_SECTION((const HiddenMarkovModel& getHMM() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((void generateModel()))
{
  PILISNeutralLossModel model;
	model.generateModel();
	TEST_NOT_EQUAL(model.getHMM().getNumberOfStates(), 0)
}
END_SECTION

START_SECTION((void evaluate()))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((PILISNeutralLossModel& operator=(const PILISNeutralLossModel &mode)))
{
  HiddenMarkovModel hmm1;
  hmm1.setPseudoCounts(13);
  PILISNeutralLossModel model;
  TEST_REAL_SIMILAR(model.getHMM().getPseudoCounts(), HiddenMarkovModel().getPseudoCounts())
  model.setHMM(hmm1);
  TEST_REAL_SIMILAR(model.getHMM().getPseudoCounts(), 13)

	PILISNeutralLossModel model2;
	model2 = model;
	TEST_REAL_SIMILAR(model2.getHMM().getPseudoCounts(), 13)

	PILISNeutralLossModel model3;
	Param p(model.getParameters());
	p.setValue("ion_name", "y");
	model.setParameters(p);
	model3 = model;
	TEST_EQUAL(model.getParameters() == model3.getParameters(), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



