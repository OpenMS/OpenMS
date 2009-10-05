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

#include <iostream>

#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(ProtonDistributionModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ProtonDistributionModel* ptr = 0;

START_SECTION(ProtonDistributionModel())
	ptr = new ProtonDistributionModel();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~ProtonDistributionModel())
	delete ptr;
END_SECTION

ptr = new ProtonDistributionModel();

START_SECTION(ProtonDistributionModel(const ProtonDistributionModel& model))
	ProtonDistributionModel copy(*ptr);
	NOT_TESTABLE
END_SECTION

START_SECTION(ProtonDistributionModel& operator = (const ProtonDistributionModel& pdm))
	ProtonDistributionModel copy;
	copy = *ptr;
	NOT_TESTABLE
END_SECTION

START_SECTION(void getProtonDistribution(vector<DoubleReal>& bb_charges, vector<DoubleReal>& sc_charges, const AASequence& peptide, Int charge, Residue::ResidueType res_type = Residue::YIon))
	vector<DoubleReal> bb_charges, sc_charges;
	DoubleReal bb_tmp[] = {1.76496e-09, 2.9459e-13, 6.3724e-12, 2.96724e-13, 0.69332e-13, 6.56286e-13, 4.82365e-13, 3.51139e-13, 5.82514e-23, 1.35049e-12};
	AASequence peptide("DFPIANGER");
	ptr->getProtonDistribution(bb_charges, sc_charges, peptide, 1);
	for (Size i = 0; i <= peptide.size(); ++i)
	{
		TEST_REAL_SIMILAR(bb_charges[i], bb_tmp[i])
	}

	DoubleReal sc_tmp[] = {2.7239e-23, 0, 0, 0, 0, 7.77547e-15, 0, 1.15343e-22, 1};
	for (Size i = 0; i != peptide.size(); ++i)
	{
		TEST_REAL_SIMILAR(sc_charges[i], sc_tmp[i])
	}

END_SECTION

START_SECTION((void setPeptideProtonDistribution(const std::vector< DoubleReal > &bb_charge, const std::vector< DoubleReal > &sc_charge)))
	vector<DoubleReal> bb_charges, sc_charges;
	AASequence peptide("DFPIANGER");
	ptr->getProtonDistribution(bb_charges, sc_charges, peptide, 1);

	ptr->setPeptideProtonDistribution(bb_charges, sc_charges);
	NOT_TESTABLE
END_SECTION

START_SECTION((void getChargeStateIntensities(const AASequence &peptide, const AASequence &n_term_ion, const AASequence &c_term_ion, Int charge, Residue::ResidueType n_term_type, std::vector< DoubleReal > &n_term_intensities, std::vector< DoubleReal > &c_term_intensities, FragmentationType type)))
	vector<DoubleReal> bb_charges, sc_charges;
	AASequence peptide("DFPIANGER");
	ptr->getProtonDistribution(bb_charges, sc_charges, peptide, 1);
	
	// set the full proton distribution
	ptr->setPeptideProtonDistribution(bb_charges, sc_charges);

	//DoubleReal n_term1(0), n_term2(0), c_term1(0), c_term2(0);
	AASequence pre1("DFP"), suf1("IANGER");
	vector<DoubleReal> pre_ints, suf_ints;
	ptr->getChargeStateIntensities(peptide, pre1, suf1, 1, Residue::YIon, pre_ints, suf_ints, ProtonDistributionModel::ChargeDirected);

	//TEST_REAL_SIMILAR(n_term1, 0.0);
	//TEST_REAL_SIMILAR(n_term2, 0.0);
	//TEST_REAL_SIMILAR(c_term1, 1.0);
	//TEST_REAL_SIMILAR(c_term2, 0.0);

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
