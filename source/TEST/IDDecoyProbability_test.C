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
// $Maintainer: Andreas Bertsch, Sven Nahnsen$
// $Authors: Andreas Bertsch $
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
IDDecoyProbability* nullPointer = 0;
START_SECTION(IDDecoyProbability())
{
	ptr = new IDDecoyProbability();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~IDDecoyProbability())
{
	delete ptr;
}
END_SECTION

START_SECTION((IDDecoyProbability(const IDDecoyProbability &rhs)))
{
 	NOT_TESTABLE 
}
END_SECTION

START_SECTION((IDDecoyProbability& operator=(const IDDecoyProbability &rhs)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void apply(std::vector<PeptideIdentification>& prob_ids, const std::vector< PeptideIdentification > &fwd_ids, const std::vector< PeptideIdentification > &rev_ids)))
{
  IDDecoyProbability decoy;
	vector<ProteinIdentification> prot_ids_fwd, prot_ids_rev;
	vector<PeptideIdentification> pep_ids_fwd, pep_ids_rev, prob_ids;
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), prot_ids_fwd, pep_ids_fwd, document_id);
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), prot_ids_rev, pep_ids_rev, document_id);

	decoy.apply(prob_ids, pep_ids_fwd, pep_ids_rev);

	for (vector<PeptideIdentification>::const_iterator it = prob_ids.begin(); it != prob_ids.end(); ++it)
	{
		if (it->getHits().size() > 0)
		{
			for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
			{
				DoubleReal prob(pit->getScore());
				DoubleReal orig_score((DoubleReal)pit->getMetaValue("XTandem_score"));
				if (orig_score > 40.0)
				{
					TEST_EQUAL(prob > 0.9, true)
				}
				if (orig_score < 20)
				{
					TEST_EQUAL(prob < 0.05, true)
				}
			}
		}
	}
}
END_SECTION

START_SECTION((void apply(std::vector< PeptideIdentification > &ids)))
{
	IDDecoyProbability decoy;
	vector<ProteinIdentification> prot_ids_fwd, prot_ids_rev, prot_ids;
  vector<PeptideIdentification> pep_ids_fwd, pep_ids_rev, prob_ids, pep_ids;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), prot_ids_fwd, pep_ids_fwd, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), prot_ids_rev, pep_ids_rev, document_id);

  for (vector<PeptideIdentification>::iterator it = pep_ids_fwd.begin(); it != pep_ids_fwd.end(); ++it)
  {
    vector<PeptideHit> hits = it->getHits();
    for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
    {
      pit->setMetaValue("target_decoy", "target");
    }
    it->setHits(hits);
    pep_ids.push_back(*it);
  }
  for (vector<PeptideIdentification>::iterator it = pep_ids_rev.begin(); it != pep_ids_rev.end(); ++it)
  {
    vector<PeptideHit> hits = it->getHits();
    for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
    {
      pit->setMetaValue("target_decoy", "decoy");
    }
    it->setHits(hits);
    pep_ids.push_back(*it);
  }

  decoy.apply(pep_ids);

  for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
  {
    if (it->getHits().size() > 0)
    {
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        DoubleReal prob(pit->getScore());
        DoubleReal orig_score((DoubleReal)pit->getMetaValue("XTandem_score"));
        if (orig_score > 40.0)
        {
          TEST_EQUAL(prob > 0.9, true)
        }
        if (orig_score < 20)
        {
          TEST_EQUAL(prob < 0.05, true)
        }
      }
    }
  }
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



