// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIonScoringCID, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CompNovoIonScoringCID* ptr = 0;
CompNovoIonScoringCID* nullPointer = 0;
START_SECTION(CompNovoIonScoringCID())
{
	ptr = new CompNovoIonScoringCID();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((CompNovoIonScoringCID(const CompNovoIonScoringCID &source)))
{
  CompNovoIonScoringCID cnis;
	Param p(cnis.getParameters());
	p.setValue("fragment_mass_tolerance", 0.6);
	cnis.setParameters(p);
	TEST_EQUAL(CompNovoIonScoringCID(cnis).getParameters() == p, true)
}
END_SECTION

START_SECTION((virtual ~CompNovoIonScoringCID()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void scoreSpectrum(Map<DoubleReal, IonScore>& CID_ion_scores, PeakSpectrum& CID_spec, DoubleReal precursor_weight, Size charge)))
{
  Map<DoubleReal, CompNovoIonScoringBase::IonScore> ion_scores;
	TheoreticalSpectrumGenerator tsg;
	Param tsg_param(tsg.getParameters());
	tsg_param.setValue("add_losses", "true");
	tsg_param.setValue("add_isotopes", "true");
	tsg.setParameters(tsg_param);

	RichPeakSpectrum rspec;
	tsg.getSpectrum(rspec, AASequence("DFPIANGER"));

	PeakSpectrum spec;
	for (Size i = 0; i != rspec.size(); ++i)
	{
		Peak1D p;
		p.setMZ(rspec[i].getMZ());
		p.setIntensity(rspec[i].getIntensity());
		spec.push_back(p);
	}

	CompNovoIonScoringCID cnis;
	cnis.scoreSpectrum(ion_scores, spec, 1018.48, 1);
	
	for (Map<DoubleReal, CompNovoIonScoringBase::IonScore>::ConstIterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
	{
/*
y1 175.118952187571
y2 304.161545285171
y3 361.183009010571
y4 475.225936461371
y5 546.263050250571
y6 659.347114231171
y7 756.399878084171
y8 903.468292000971

b1 117.042044532471
b2 263.102633417371
b3 360.155397270371
b4 473.239461250971
b5 544.276575040171
b6 658.319502490971
b7 715.340966216371
b8 844.383559313971

*/	
		cerr << it->first << " " << it->second.score << endl;
		if (fabs(it->first - 903.468292000971) < 0.01 || 
				fabs(it->first - 756.399878084171) < 0.01 ||
				fabs(it->first - 659.347114231171) < 0.01 ||
				fabs(it->first - 546.263050250571) < 0.01 ||
				fabs(it->first - 475.225936461371) < 0.01 ||
				fabs(it->first - 361.183009010571) < 0.01 ||
				fabs(it->first - 304.161545285171) < 0.01 ||
				fabs(it->first - 175.118952187571) < 0.01 ||
				fabs(it->first - 263.102633417371) < 0.01 ||
				fabs(it->first - 360.155397270371) < 0.01 ||
				/*fabs(it->first - 473.239461250971) < 0.01 ||*/
				/*fabs(it->first - 544.276575040171) < 0.01 ||*/
				fabs(it->first - 658.319502490971) < 0.01 ||
				fabs(it->first - 715.340966216371) < 0.01 ||
				fabs(it->first - 844.383559313971) < 0.01)
		{
			TEST_EQUAL(it->second.score > 1, true)
		}
		else
		{
			TEST_EQUAL(it->second.score <= 1, true)
		}
	}
}
END_SECTION

START_SECTION((CompNovoIonScoringCID& operator=(const CompNovoIonScoringCID &source)))
{
  CompNovoIonScoringCID cnis;
	Param p(cnis.getParameters());
	p.setValue("fragment_mass_tolerance", 0.6);
	cnis.setParameters(p);
	CompNovoIonScoringCID cnis2;
	cnis2 = cnis;
	TEST_EQUAL(cnis2.getParameters() == cnis.getParameters(), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



