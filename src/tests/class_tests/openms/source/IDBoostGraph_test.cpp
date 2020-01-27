// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2019.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/test_config.h>

using namespace OpenMS;
using namespace std;
using Internal::IDBoostGraph;

START_TEST(IDBoostGraph, "$Id$")

    START_SECTION(IDBoostGraph only best PSMs)
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

    START_SECTION(IDBoostGraph all PSMs)
        {
          vector<ProteinIdentification> prots;
          vector<PeptideIdentification> peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
          IDBoostGraph idb{prots[0], peps, 0, false};
          TEST_EQUAL(idb.getNrConnectedComponents(), 0)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 14)
          idb.computeConnectedComponents();
          // Now it is 5 ccs because there is an unmatched peptide and a new PSM that only matches to
          // previously uncovered protein PH2.
          TEST_EQUAL(idb.getNrConnectedComponents(), 5)
          // The next lines do not sum up to 9 because protein PH2 is unmatched
          // If you want to avoid that, filter unmatched proteins first!
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

          IDBoostGraph idb{prots[0], peps, 1, true};
          TEST_EQUAL(idb.getNrConnectedComponents(), 0)
          TEST_EQUAL(boost::num_vertices(idb.getComponent(0)), 8)
          idb.computeConnectedComponents();
          TEST_EQUAL(idb.getNrConnectedComponents(), 2)
          // The next lines do not sum up to 9 because protein PH2 is unmatched
          // If you want to avoid that, filter unmatched proteins first!
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

    START_SECTION(IDBoostGraph on consensusXML TODO)
    {

    }
    END_SECTION

END_TEST
