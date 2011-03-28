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
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationCID.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIdentificationCID, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CompNovoIdentificationCID* ptr = 0;
CompNovoIdentificationCID* nullPointer = 0;
START_SECTION(CompNovoIdentificationCID())
{
	ptr = new CompNovoIdentificationCID();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CompNovoIdentificationCID())
{
	delete ptr;
}
END_SECTION

START_SECTION((CompNovoIdentificationCID(const CompNovoIdentificationCID& source)))
	CompNovoIdentificationCID cnis;
	Param p(cnis.getParameters());
	p.setValue("precursor_mass_tolerance", 3.0);
	cnis.setParameters(p);
	TEST_EQUAL(CompNovoIdentificationCID(cnis).getParameters() == p, true)
END_SECTION

START_SECTION((CompNovoIdentificationCID& operator = (const CompNovoIdentificationCID& source)))
  CompNovoIdentificationCID cnis;
  Param p(cnis.getParameters());
  p.setValue("precursor_mass_tolerance", 3.0);
  cnis.setParameters(p);
	CompNovoIdentificationCID cnis2;
	cnis2 = cnis;
  TEST_EQUAL(cnis2.getParameters() == p, true)
END_SECTION

START_SECTION((void getIdentifications(std::vector<PeptideIdentification>& ids, const PeakMap& exp)))
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

  Precursor prec;
  prec.setMZ((AASequence("DFPLANGER").getMonoWeight() + 2.0 * Constants::PROTON_MASS_U) / 2.0);
  prec.setCharge(2);
  vector<Precursor> precs;
  precs.push_back(prec);
  spec.setPrecursors(precs);

  vector<PeptideIdentification> ids;
  CompNovoIdentificationCID cni;
  Param cni_param(cni.getParameters());
  cni_param.setValue("precursor_mass_tolerance", 0.3);
  cni.setParameters(cni_param);
	PeakMap exp;
	exp.push_back(spec);
  cni.getIdentifications(ids, exp);
  TEST_EQUAL(ids.size(), 1)
  TEST_EQUAL(ids.begin()->getHits().size() > 0, true)
  TEST_EQUAL(ids.begin()->getHits().begin()->getSequence() == AASequence("DFPLANGER"), true)
END_SECTION

START_SECTION((void getIdentification(PeptideIdentification& id, const PeakSpectrum& CID_spec)))
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

	Precursor prec;
	prec.setMZ((AASequence("DFPLANGER").getMonoWeight() + 2.0 * Constants::PROTON_MASS_U) / 2.0);
	prec.setCharge(2);
	vector<Precursor> precs;
	precs.push_back(prec);
	spec.setPrecursors(precs);	

	PeptideIdentification id;
	CompNovoIdentificationCID cni;
	Param cni_param(cni.getParameters());
	cni_param.setValue("precursor_mass_tolerance", 0.3);
	cni.setParameters(cni_param);
	cni.getIdentification(id, spec);
	TEST_EQUAL(id.getHits().size() > 0, true)
	TEST_EQUAL(id.getHits().begin()->getSequence() == AASequence("DFPLANGER"), true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



