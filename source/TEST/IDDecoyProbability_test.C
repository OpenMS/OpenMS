// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch, Sven Nahnsen$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IDDecoyProbability, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IDDecoyProbability* ptr = 0;
CHECK(IDDecoyProbability())
{
	ptr = new IDDecoyProbability();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(virtual ~IDDecoyProbability())
{
	delete ptr;
}
RESULT

CHECK((IDDecoyProbability(const IDDecoyProbability &rhs)))
{
 	NOT_TESTABLE 
}
RESULT

CHECK((IDDecoyProbability& operator=(const IDDecoyProbability &rhs)))
{
  NOT_TESTABLE
}
RESULT

CHECK((void apply(std::vector<PeptideIdentification>& prob_ids, const std::vector< PeptideIdentification > &fwd_ids, const std::vector< PeptideIdentification > &rev_ids)))
{
  IDDecoyProbability decoy;
	vector<ProteinIdentification> prot_ids_fwd, prot_ids_rev;
	vector<PeptideIdentification> pep_ids_fwd, pep_ids_rev, prob_ids;
	IdXMLFile().load("data/XTandem_fwd_ids.idXML", prot_ids_fwd, pep_ids_fwd);
	IdXMLFile().load("data/XTandem_rev_ids.idXML", prot_ids_rev, pep_ids_rev);

	decoy.apply(prob_ids, pep_ids_fwd, pep_ids_rev);

	for (vector<PeptideIdentification>::const_iterator it = prob_ids.begin(); it != prob_ids.end(); ++it)
	{
		if (it->getHits().size() > 0)
		{
			for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
			{
				double prob(pit->getScore());
				double orig_score((double)pit->getMetaValue("XTandem_score"));

				if (orig_score > 40.0)
				{
					TEST_EQUAL(prob > 0.97, true)
				}
				if (orig_score < 25)
				{
					TEST_EQUAL(prob < 0.05, true)
				}
			}
		}
	}
}
RESULT


CHECK(void generateDistributionImage(const Map< double, double > &ids, const String &formula, const String &filename))
	NOT_TESTABLE
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



