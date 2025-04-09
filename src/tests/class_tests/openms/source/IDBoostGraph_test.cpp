// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/test_config.h>

using namespace OpenMS;
using namespace std;
using Internal::IDBoostGraph;

static void runIBGResolve(vector<ProteinIdentification>& inferred_protein_ids,
                          vector<PeptideIdentification>& inferred_peptide_ids)
{
  IDBoostGraph ibg(inferred_protein_ids[0], inferred_peptide_ids, 1, false, false);
  ibg.computeConnectedComponents();
  ibg.clusterIndistProteinsAndPeptides(); //TODO check in resolve or do it there if not done yet!
  //Note: the above does not add singleton groups to graph
  ibg.resolveGraphPeptideCentric(true);
  inferred_protein_ids[0].getIndistinguishableProteins().clear();
  inferred_protein_ids[0].getProteinGroups().clear();
  ibg.annotateIndistProteins(true); // this does not really add singletons since they are not in the graph
  IDFilter::removeUnreferencedProteins(inferred_protein_ids[0], inferred_peptide_ids);
  IDFilter::updateProteinGroups(inferred_protein_ids[0].getIndistinguishableProteins(),inferred_protein_ids[0].getHits());
  inferred_protein_ids[0].fillIndistinguishableGroupsWithSingletons();
  auto & ipg = inferred_protein_ids[0].getIndistinguishableProteins();
  std::sort(std::begin(ipg), std::end(ipg));
}

START_TEST(IDBoostGraph, "$Id$")

    START_SECTION(IDBoostGraph only best PSMs)
    {
      vector<ProteinIdentification> prots;
      vector<PeptideIdentification> peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      IDBoostGraph idb{prots[0], peps, 1, false, false};
      TEST_EQUAL(idb.getNrConnectedComponents(), 0)
      // 6 proteins (1 unmatched and omitted since we build the graph psm-centric) plus 4 peptides (top per psm).
      TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 9)
      idb.computeConnectedComponents();
      TEST_EQUAL(idb.getNrConnectedComponents(), 3)
      TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 3)
      TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 4)
      TEST_EQUAL(boost::num_vertices(idb.getComponent(2)), 2)
      TEST_EXCEPTION(Exception::MissingInformation, idb.clusterIndistProteinsAndPeptidesAndExtendGraph());
      idb.clusterIndistProteinsAndPeptides();
      TEST_EQUAL(idb.getNrConnectedComponents(), 3)
      // Only cc 0 and 1 have indist prot group
      TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 4)
      TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 5)
      TEST_EQUAL(boost::num_vertices(idb.getComponent(2)), 2)
    }
    END_SECTION

    /* TODO test graph-based resolution
    START_SECTION(IDBoostGraph graph-based group resolution)
        {
          vector<ProteinIdentification> prots;
          vector<PeptideIdentification> peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
          IDBoostGraph idb{prots[0], peps, 1, false};
          TEST_EQUAL(idb.getNrConnectedComponents(), 0)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 9)
          idb.computeConnectedComponents();
          TEST_EQUAL(idb.getNrConnectedComponents(), 3)
          // The next lines do not sum up to 9 because protein PH2 is unmatched
          // If you want to avoid that, filter unmatched proteins first!
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 3)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 4)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(2)), 2)
          TEST_EXCEPTION(Exception::MissingInformation, idb.clusterIndistProteinsAndPeptidesAndExtendGraph());
          idb.clusterIndistProteinsAndPeptides();
          TEST_EQUAL(idb.getNrConnectedComponents(), 3)
          // Only cc 0 and 1 have indist prot group
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 4)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 5)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(2)), 2)
        }
    END_SECTION
     */

    START_SECTION(IDBoostGraph all PSMs)
        {
          vector<ProteinIdentification> prots;
          vector<PeptideIdentification> peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
          IDBoostGraph idb{prots[0], peps, 0, false, false};
          TEST_EQUAL(idb.getNrConnectedComponents(), 0)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 14)
          idb.computeConnectedComponents();
          // Now it is 5 ccs because there is an unmatched peptide and a new PSM that only matches to
          // previously uncovered protein PH2.
          TEST_EQUAL(idb.getNrConnectedComponents(), 5)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 4)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 2)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(2)), 5)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(3)), 1)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(4)), 2)
        }
    END_SECTION

    START_SECTION(IDBoostGraph only best PSMs with runinfo)
        {
          vector<ProteinIdentification> prots;
          vector<PeptideIdentification> peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("IDBoostGraph_test_in.idXML"),prots,peps);

          IDBoostGraph idb{prots[0], peps, 1, true, false};
          TEST_EQUAL(idb.getNrConnectedComponents(), 0)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 8)
          idb.computeConnectedComponents();
          TEST_EQUAL(idb.getNrConnectedComponents(), 2)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 6)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 2)

          idb.clusterIndistProteinsAndPeptidesAndExtendGraph();

          TEST_EQUAL(idb.getNrConnectedComponents(), 2)
          // Only cc 0 and 1 have indist prot group
          //TODO we could reduce the number of nodes by removing ones without evidence
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 25)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 11)

        }
    END_SECTION

    START_SECTION(IDBoostGraph graph-based group resolution)
        {
          vector<ProteinIdentification> prots;
          vector<PeptideIdentification> peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
          IDBoostGraph idb{prots[0], peps, 1, false, false};
          TEST_EQUAL(idb.getNrConnectedComponents(), 0)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 9)
          idb.computeConnectedComponents();
          TEST_EQUAL(idb.getNrConnectedComponents(), 3)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 3)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 4)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(2)), 2)
          TEST_EXCEPTION(Exception::MissingInformation, idb.clusterIndistProteinsAndPeptidesAndExtendGraph());
          idb.clusterIndistProteinsAndPeptides();
          // Only cc 0 and 1 have indist prot group
          TEST_EQUAL(boost::num_edges(idb.getComponent(0)), 3)
          TEST_EQUAL(boost::num_edges(idb.getComponent(1)), 4)
          TEST_EQUAL(boost::num_edges(idb.getComponent(2)), 1)
          idb.resolveGraphPeptideCentric();
          TEST_EQUAL(idb.getNrConnectedComponents(), 3)
          // Only cc 0 and 1 have indist prot group
          TEST_EQUAL(boost::num_edges(idb.getComponent(0)), 3)
          // There is one shared peptide in the second component whose edge will be resolved
          TEST_EQUAL(boost::num_edges(idb.getComponent(1)), 3)
          TEST_EQUAL(boost::num_edges(idb.getComponent(2)), 1)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 4)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(1)), 5)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(2)), 2)
        }
    END_SECTION


    START_SECTION(Resolution)
    {
      // TODO problem is that there is no way to build the graph using existing groups.
      //  therefore resolution on the graph will redo groups and assign new scores.
      //  Therefore we need slightly different test files.
      vector<ProteinIdentification> prots;
      vector<PeptideIdentification> peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_in.idXML"), prots, peps);
      runIBGResolve(prots, peps);
      std::string tmp_filename;
      NEW_TMP_FILE(tmp_filename);
      IdXMLFile().store(tmp_filename, prots, peps);
      TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_out_ibg.idXML"), tmp_filename);

      prots.clear();
      peps.clear();
      tmp_filename.clear();
      NEW_TMP_FILE(tmp_filename);
      idf.load(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_in2.idXML"), prots, peps);
      runIBGResolve(prots, peps);
      IdXMLFile().store(tmp_filename, prots, peps);
      TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_out2_ibg.idXML"), tmp_filename);
    }
    END_SECTION

    START_SECTION(IDBoostGraph on consensusXML TODO)
    {

    }
    END_SECTION

END_TEST
