// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FalseDiscoveryRate, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FalseDiscoveryRate* ptr = nullptr;
FalseDiscoveryRate* nullPointer = nullptr;
START_SECTION(FalseDiscoveryRate())
{
  ptr = new FalseDiscoveryRate();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FalseDiscoveryRate())
{
  delete ptr;
}
END_SECTION

START_SECTION((void apply(std::vector<PeptideIdentification> &fwd_ids, std::vector<PeptideIdentification> &rev_ids)))
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
    if (!it->getHits().empty())
    {
      PeptideHit hit(*it->getHits().begin());
      double fdr(hit.getScore());
      double orig_score((double)hit.getMetaValue("XTandem_score"));
      
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

START_SECTION((void apply(std::vector<ProteinIdentification> &fwd_ids, std::vector<ProteinIdentification> &rev_ids)))
{
  vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids;
  vector<PeptideIdentification> fwd_pep_ids, rev_pep_ids;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids_withProtScores.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids_withProtScores.idXML"), rev_prot_ids, rev_pep_ids, document_id);
  ptr->apply(fwd_prot_ids, rev_prot_ids);
  TOLERANCE_ABSOLUTE(0.001)
	
  for (vector<ProteinIdentification>::const_iterator prot_it = fwd_prot_ids.begin(); prot_it != fwd_prot_ids.end(); ++prot_it)
  {
    if (!prot_it->getHits().empty())
    {
      for (vector<ProteinHit>::const_iterator it = prot_it->getHits().begin(); it != prot_it->getHits().end(); ++it)
      {
        ProteinHit hit(*it);
        double fdr(hit.getScore());
        double orig_score((double)hit.getMetaValue("XTandem_score"));
        
        // it gets here, but neither of the conditions below are ever satisfied
        if (orig_score < -1.8)
        {
          TEST_REAL_SIMILAR(fdr, 0)
        }
        if (orig_score == 0.0)
        {
          TEST_REAL_SIMILAR(fdr, 0.897384)
        }
      }
    }
  }
}
END_SECTION

START_SECTION((void apply(std::vector<PeptideIdentification> &id)))
{
  vector<ProteinIdentification> prot_ids;
  vector<PeptideIdentification> pep_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FalseDiscoveryRate_OMSSA.idXML"), prot_ids, pep_ids);

  ptr->apply(pep_ids);
  TOLERANCE_ABSOLUTE(0.001)
	
  for (Size z = 1; z <= 4; ++z)
  {
    for (vector<PeptideIdentification>::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
  	{
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
      	double fdr(pit->getScore());
      	double orig_score((double)pit->getMetaValue("OMSSA_score"));

      	if (orig_score <= 10e-4)
        {
          TEST_REAL_SIMILAR(fdr, 0)
        }
      	if (orig_score >= 1000 && pit->getCharge() != 1)
        {
          TEST_EQUAL(fdr > 0.1, true)
        }
      }
    }

    // target hit
    PeptideIdentification pep_id = pep_ids[0];
    PeptideHit pit = pep_id.getHits()[0];
    double fdr(pit.getScore());
    TEST_REAL_SIMILAR(fdr, 0.0730478589420655);

    // target+decoy hit considered as target
    pep_id = pep_ids[5];
    pit = pep_id.getHits()[0];
    fdr = pit.getScore();
    TEST_REAL_SIMILAR(fdr, 0.409926470588235);

    // decoy hit removed
    pep_id = pep_ids[9];
    TEST_EQUAL(pep_id.getHits().size(), 0)
  }
}
END_SECTION

START_SECTION((void apply(std::vector<ProteinIdentification>& ids)))
{
  vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids, prot_ids;
  vector<PeptideIdentification> fwd_pep_ids, rev_pep_ids, pep_ids;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);

  for (vector<ProteinIdentification>::const_iterator it = fwd_prot_ids.begin(); it != fwd_prot_ids.end(); ++it)
  {
    prot_ids.push_back(*it);
    for (vector<ProteinHit>::iterator hit_it = prot_ids.back().getHits().begin(); hit_it != prot_ids.back().getHits().end(); ++hit_it)
    {
      hit_it->setMetaValue("target_decoy", "target");
    }
  }
  
  for (vector<ProteinIdentification>::const_iterator it = rev_prot_ids.begin(); it != rev_prot_ids.end(); ++it)
  {
    prot_ids.push_back(*it);
    for (vector<ProteinHit>::iterator hit_it = prot_ids.back().getHits().begin(); hit_it != prot_ids.back().getHits().end(); ++hit_it)
    {
      hit_it->setMetaValue("target_decoy", "decoy");
    }
  }

  ptr->apply(prot_ids);

  TOLERANCE_ABSOLUTE(0.001)

  for (vector<ProteinIdentification>::const_iterator prot_it = prot_ids.begin(); prot_it != prot_ids.end(); ++prot_it)
  {
    if (!prot_it->getHits().empty())
    {
      for (vector<ProteinHit>::const_iterator it = prot_it->getHits().begin(); it != prot_it->getHits().end(); ++it)
      {
        ProteinHit hit(*it);
        double fdr(hit.getScore());
        double orig_score((double)hit.getMetaValue("XTandem_score"));

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
          TEST_EQUAL(fdr > 0.1, true)
        }
      }
    }
  }
}
END_SECTION

START_SECTION((void applyPicked(std::vector<ProteinIdentification>& ids)))
{
  vector<ProteinIdentification> prot_ids;
  vector<PeptideIdentification> pep_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FalseDiscoveryRate_picked_in.idXML"), prot_ids, pep_ids);

  ptr->applyPickedProteinFDR(prot_ids[0],"decoy_");

  TOLERANCE_ABSOLUTE(0.001)
  const auto& hits = prot_ids[0].getHits();
  TEST_REAL_SIMILAR(hits[0].getScore(),0.25)
  TEST_REAL_SIMILAR(hits[1].getScore(),0.25)
  TEST_REAL_SIMILAR(hits[2].getScore(),0.25)
  TEST_REAL_SIMILAR(hits[3].getScore(),0.4)
  TEST_REAL_SIMILAR(hits[4].getScore(),0.4)
  TEST_REAL_SIMILAR(hits[5].getScore(),0.5)
}
END_SECTION


START_SECTION((void apply(std::vector<ProteinIdentification>& ids)))
{
  vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids, prot_ids;
  vector<PeptideIdentification> fwd_pep_ids, rev_pep_ids, pep_ids;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);

  for (vector<ProteinIdentification>::const_iterator it = fwd_prot_ids.begin(); it != fwd_prot_ids.end(); ++it)
  {
    prot_ids.push_back(*it);
    for (vector<ProteinHit>::iterator hit_it = prot_ids.back().getHits().begin(); hit_it != prot_ids.back().getHits().end(); ++hit_it)
    {
      hit_it->setMetaValue("target_decoy", "target");
    }
  }

  for (vector<ProteinIdentification>::const_iterator it = rev_prot_ids.begin(); it != rev_prot_ids.end(); ++it)
  {
    prot_ids.push_back(*it);
    for (vector<ProteinHit>::iterator hit_it = prot_ids.back().getHits().begin(); hit_it != prot_ids.back().getHits().end(); ++hit_it)
    {
      hit_it->setMetaValue("target_decoy", "decoy");
    }
  }

  ptr->apply(prot_ids);

  TOLERANCE_ABSOLUTE(0.001)

  for (vector<ProteinIdentification>::const_iterator prot_it = prot_ids.begin(); prot_it != prot_ids.end(); ++prot_it)
  {
    if (!prot_it->getHits().empty())
    {
      for (vector<ProteinHit>::const_iterator it = prot_it->getHits().begin(); it != prot_it->getHits().end(); ++it)
      {
        ProteinHit hit(*it);
        double fdr(hit.getScore());
        double orig_score((double)hit.getMetaValue("XTandem_score"));

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
          TEST_EQUAL(fdr > 0.1, true)
        }
      }
    }
  }
}
END_SECTION
delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



