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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentification.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CompNovoIdentification* ptr = 0;
CompNovoIdentification* nullPointer = 0;
START_SECTION(CompNovoIdentification())
{
	ptr = new CompNovoIdentification();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~CompNovoIdentification())
{
	delete ptr;
}
END_SECTION

START_SECTION((CompNovoIdentification(const CompNovoIdentification& source)))
	CompNovoIdentification cni;
	Param p(cni.getParameters());
	p.setValue("fragment_mass_tolerance", 0.5);
	cni.setParameters(p);
	TEST_EQUAL(CompNovoIdentification(cni).getParameters() == p, true)
END_SECTION


START_SECTION((void getIdentifications(std::vector< PeptideIdentification > &ids, const PeakMap &exp)))
{
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

  RichPeakSpectrum rspec_ETD;
  tsg.addPeaks(rspec_ETD, AASequence("DFPIANGER"), Residue::ZIon, 1);
  tsg.addPrecursorPeaks(rspec_ETD, AASequence("DFPIANGER"), 2);
  PeakSpectrum spec_ETD;
  for (Size i = 0; i != rspec_ETD.size(); ++i)
  {
    Peak1D p;
    p.setMZ(rspec_ETD[i].getMZ());
    p.setIntensity(rspec_ETD[i].getIntensity());
    spec_ETD.push_back(p);
  }

  Precursor prec;
  prec.setMZ((AASequence("DFPLANGER").getMonoWeight() + 2.0 * Constants::PROTON_MASS_U) / 2.0);
  prec.setCharge(2);
  vector<Precursor> precs;
  precs.push_back(prec);
  spec.setPrecursors(precs);
  spec_ETD.setPrecursors(precs);

	PeakMap exp;
	exp.push_back(spec);
	exp.push_back(spec_ETD);

  vector<PeptideIdentification> ids;
  CompNovoIdentification cni;
  Param cni_param(cni.getParameters());
  cni.setParameters(cni_param);
  cni.getIdentifications(ids, exp);
  TEST_EQUAL(ids.size(), 1)
  TEST_EQUAL(ids.begin()->getHits().size() > 0, true)
  TEST_EQUAL(ids.begin()->getHits().begin()->getSequence() == AASequence("DFPLANGER"), true)  
}
END_SECTION

START_SECTION((void getIdentification(PeptideIdentification &id, const PeakSpectrum &CID_spec, const PeakSpectrum &ETD_spec)))
{
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

	RichPeakSpectrum rspec_ETD;
	tsg.addPeaks(rspec_ETD, AASequence("DFPIANGER"), Residue::ZIon, 1);
	tsg.addPrecursorPeaks(rspec_ETD, AASequence("DFPIANGER"), 2);
	PeakSpectrum spec_ETD;
	for (Size i = 0; i != rspec_ETD.size(); ++i)
  {
    Peak1D p;
    p.setMZ(rspec_ETD[i].getMZ());
    p.setIntensity(rspec_ETD[i].getIntensity());
    spec_ETD.push_back(p);
  }

  Precursor prec;
  prec.setMZ((AASequence("DFPLANGER").getMonoWeight() + 2.0 * Constants::PROTON_MASS_U) / 2.0);
  prec.setCharge(2);
  vector<Precursor> precs;
  precs.push_back(prec);
  spec.setPrecursors(precs);
	spec_ETD.setPrecursors(precs);

  PeptideIdentification id;
  CompNovoIdentification cni;
  Param cni_param(cni.getParameters());
  cni.setParameters(cni_param);
  cni.getIdentification(id, spec, spec_ETD);
  TEST_EQUAL(id.getHits().size() > 0, true)
  TEST_EQUAL(id.getHits().begin()->getSequence() == AASequence("DFPLANGER"), true)

}
END_SECTION

START_SECTION((CompNovoIdentification& operator=(const CompNovoIdentification &source)))
{
  CompNovoIdentification cni;
  Param p(cni.getParameters());
  p.setValue("fragment_mass_tolerance", 0.5);
  cni.setParameters(p);
	CompNovoIdentification cni2;
	cni2 = cni;
  TEST_EQUAL(cni2.getParameters() == p, true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



