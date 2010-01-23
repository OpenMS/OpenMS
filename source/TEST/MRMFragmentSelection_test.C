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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MRMFragmentSelection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFragmentSelection* ptr = 0;
START_SECTION(MRMFragmentSelection())
{
	ptr = new MRMFragmentSelection();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~MRMFragmentSelection())
{
	delete ptr;
}
END_SECTION

START_SECTION((MRMFragmentSelection(const MRMFragmentSelection &rhs)))
{
  MRMFragmentSelection mrmfs;
	Param p = mrmfs.getParameters();
	p.setValue("num_top_peaks", 18);
	mrmfs.setParameters(p);
	TEST_EQUAL(MRMFragmentSelection(mrmfs).getParameters() == p, true)
}
END_SECTION

START_SECTION((MRMFragmentSelection& operator=(const MRMFragmentSelection &rhs)))
{
  MRMFragmentSelection mrmfs;
	Param p = mrmfs.getParameters();
	p.setValue("num_top_peaks", 18);
	mrmfs.setParameters(p);
	MRMFragmentSelection mrmfs2;
	mrmfs2 = mrmfs;
	TEST_EQUAL(mrmfs2.getParameters() == p, true)
}
END_SECTION

START_SECTION((void selectFragments(std::vector< RichPeak1D > &selected_peaks, const RichPeakSpectrum &spec)))
{
	RichPeakSpectrum spec;
	TheoreticalSpectrumGenerator tsg;
	Param tsg_param(tsg.getParameters());
	tsg_param.setValue("add_metainfo", "true");
	tsg.setParameters(tsg_param);
	tsg.addPeaks(spec, AASequence("DFPIANGER"), Residue::YIon, 1);
	tsg.addPeaks(spec, AASequence("DFPIANGER"), Residue::BIon, 1);

	spec.sortByPosition();
	Precursor prec;
	prec.setMZ(1019.1);
	vector<Precursor> precursors;
	precursors.push_back(prec);
	spec.setPrecursors(precursors);
	
	PeptideHit hit;
	hit.setCharge(1);
	hit.setSequence("DFPIANGER");
	vector<PeptideHit> hits;
	hits.push_back(hit);
	PeptideIdentification id;
	id.setHits(hits);
	vector<PeptideIdentification> ids;
	ids.push_back(id);
	spec.setPeptideIdentifications(ids);

	MRMFragmentSelection mrmfs;
	Param p(mrmfs.getParameters());
	p.setValue("num_top_peaks", 1);
	p.setValue("allowed_ion_types", StringList::create("y"));
	mrmfs.setParameters(p);

	vector<RichPeak1D> selected_peaks;
	mrmfs.selectFragments(selected_peaks, spec);
	TEST_EQUAL(selected_peaks.size(), 1)

	p.setValue("num_top_peaks", 3);
	p.setValue("min_pos_precursor_percentage", 10.0);
	mrmfs.setParameters(p);
	selected_peaks.clear();
	mrmfs.selectFragments(selected_peaks, spec);
	TEST_EQUAL(selected_peaks.size(), 3)

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



