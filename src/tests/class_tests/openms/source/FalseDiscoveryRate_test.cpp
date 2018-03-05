// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    if (it->getHits().size() > 0)
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
    if (prot_it->getHits().size() > 0)
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
    if (prot_it->getHits().size() > 0)
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



