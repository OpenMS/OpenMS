// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

START_TEST(PILISNeutralLossModel, "$Id: PILISNeutralLossModel_test.C 6512 2010-01-10 18:14:45Z andreas_bertsch $")

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
START_SECTION(PILISNeutralLossModel())
{
	ptr = new PILISNeutralLossModel();
	TEST_NOT_EQUAL(ptr, 0)
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
	model.getIons(peaks1, "ANGER", 1.0);
	TEST_EQUAL(peaks1.size(), 9)

/*
	for (Size i = 0; i != peaks1.size(); ++i)
	{
		cout << peaks1[i].getMZ() << " " << peaks1[i].getIntensity() << " " << peaks1[i].getMetaValue("IonName") << endl;
	}
*/

	model.evaluate();

	vector<RichPeak1D> peaks2;
	model.getIons(peaks2, "ANGER", 1.0);
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



