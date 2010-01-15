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
#include <OpenMS/FORMAT/IdXMLFile.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FalseDiscoveryRate, "$Id: FalseDiscoveryRate_test.C 6446 2009-11-20 16:21:41Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FalseDiscoveryRate* ptr = 0;
START_SECTION(FalseDiscoveryRate())
{
	ptr = new FalseDiscoveryRate();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~FalseDiscoveryRate())
{
	delete ptr;
}
END_SECTION

START_SECTION((void apply(std::vector< PeptideIdentification > &fwd_ids, std::vector< PeptideIdentification > &rev_ids)))
{
  ptr = new FalseDiscoveryRate();
	vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids;
	vector<PeptideIdentification> fwd_pep_ids, rev_pep_ids;
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);
	ptr->apply(fwd_pep_ids, rev_pep_ids);
	TOLERANCE_ABSOLUTE(0.0001)
	for (vector<PeptideIdentification>::const_iterator it = fwd_pep_ids.begin(); it != fwd_pep_ids.end(); ++it)
	{
		if (it->getHits().size() > 0)
		{
			PeptideHit hit(*it->getHits().begin());
			DoubleReal fdr(hit.getScore());
			DoubleReal orig_score((DoubleReal)hit.getMetaValue("XTandem_score"));

			if (orig_score >= 39.4)
			{
				TEST_REAL_SIMILAR(fdr, 0)
			}
			if (orig_score <= 37.9 + 0.0001 && orig_score >= 37.9 - 0.0001)
			{
				TEST_REAL_SIMILAR(fdr, 0.08)
			}
		}
	}
}
END_SECTION

START_SECTION((void apply(std::vector< ProteinIdentification > &fwd_ids, std::vector< ProteinIdentification > &rev_ids)))
{
  vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids;
  vector<PeptideIdentification> fwd_pep_ids, rev_pep_ids;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);
  ptr->apply(fwd_prot_ids, rev_prot_ids);
  TOLERANCE_ABSOLUTE(0.001)
	
	for (vector<ProteinIdentification>::const_iterator prot_it = fwd_prot_ids.begin(); prot_it != fwd_prot_ids.end(); ++prot_it)
	{
		if (prot_it->getHits().size() > 0)
		{
			for (vector<ProteinHit>::const_iterator it = prot_it->getHits().begin(); it != prot_it->getHits().end(); ++it)
			{
				ProteinHit hit(*it);
				DoubleReal fdr(hit.getScore());
				DoubleReal orig_score((DoubleReal)hit.getMetaValue("XTandem_score"));

				if (orig_score < -1.8)
				{
					TEST_REAL_SIMILAR(fdr, 0)
				}
				if ((orig_score == -1.7))
				{
					TEST_REAL_SIMILAR(fdr, 0.0617284)
				}
			}
		}
	}
}
END_SECTION

START_SECTION((void apply(std::vector<PeptideIdentification>& ids)))
	vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids, prot_ids;
  vector<PeptideIdentification> fwd_pep_ids, rev_pep_ids, pep_ids;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);

	for (vector<PeptideIdentification>::iterator it = fwd_pep_ids.begin(); it != fwd_pep_ids.end(); ++it)
	{
		vector<PeptideHit> hits = it->getHits();
		for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
		{
			pit->setMetaValue("target_decoy", "target");
		}
		it->setHits(hits);
		pep_ids.push_back(*it);
	}
	for (vector<PeptideIdentification>::iterator it = rev_pep_ids.begin(); it != rev_pep_ids.end(); ++it)
	{
		vector<PeptideHit> hits = it->getHits();
		for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
		{
			pit->setMetaValue("target_decoy", "decoy");
		}
		it->setHits(hits);
		pep_ids.push_back(*it);
	}


  ptr->apply(pep_ids);
  TOLERANCE_ABSOLUTE(0.001)
	for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
  {
    if (it->getHits().size() > 0)
    {
      PeptideHit hit(*it->getHits().begin());
      DoubleReal fdr(hit.getScore());
      DoubleReal orig_score((DoubleReal)hit.getMetaValue("XTandem_score"));

      if (orig_score >= 39.4)
      {
        TEST_REAL_SIMILAR(fdr, 0)
      }
      if (orig_score <= 37.9 + 0.0001 && orig_score >= 37.9 - 0.0001)
      {
        TEST_REAL_SIMILAR(fdr, 0.1)
      }
    }
  }
END_SECTION

START_SECTION((void apply(std::vector<ProteinIdentification>& ids)))
	vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids, prot_ids;
  vector<PeptideIdentification> fwd_pep_ids, rev_pep_ids, pep_ids;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);

	for (vector<ProteinIdentification>::const_iterator it = fwd_prot_ids.begin(); it != fwd_prot_ids.end(); ++it)
	{
		prot_ids.push_back(*it);
	}

	for (vector<ProteinIdentification>::const_iterator it = rev_prot_ids.begin(); it != rev_prot_ids.end(); ++it)
	{
		prot_ids.push_back(*it);
	}

	ptr->apply(prot_ids);

  TOLERANCE_ABSOLUTE(0.001)

  for (vector<ProteinIdentification>::const_iterator prot_it = prot_ids.begin(); prot_it != prot_ids.end(); ++prot_it)
  {
    if (prot_it->getHits().size() > 0)
    {
      for (vector<ProteinHit>::const_iterator it = prot_it->getHits().begin(); it != prot_it->getHits().end(); ++it)
      {
        ProteinHit hit(*it);
        DoubleReal fdr(hit.getScore());
        DoubleReal orig_score((DoubleReal)hit.getMetaValue("XTandem_score"));

        if (orig_score < -1.8)
        {
          TEST_REAL_SIMILAR(fdr, 0)
        }
        if (orig_score == -1.7)
        {
          TEST_REAL_SIMILAR(fdr, 0.0617284)
        }
				if (orig_score > -1.2)
				{
					cerr << fdr << " " << orig_score << endl;
					TEST_EQUAL(fdr > 0.1, true)
				}
      }
    }
  }
	

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



