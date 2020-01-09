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
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/test_config.h>

using namespace OpenMS;
using namespace std;

START_TEST(BasicProteinInferenceAlgorithm, "$Id$")

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID)
    {
      vector<ProteinIdentification> prots;
      vector<PeptideIdentification> peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BasicProteinInferenceAlgorithm bpia;
      Param p = bpia.getParameters();
      p.setValue("min_peptides_per_protein", 0);
      bpia.setParameters(p);
      bpia.run(peps, prots);
      TEST_EQUAL(prots[0].getHits()[0].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[1].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[2].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[3].getScore(), 0.8)
      TEST_EQUAL(prots[0].getHits()[4].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[5].getScore(), 0.9)

      TEST_EQUAL(prots[0].getHits()[0].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[1].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[2].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[3].getMetaValue("nr_found_peptides"), 2)
      TEST_EQUAL(prots[0].getHits()[4].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[5].getMetaValue("nr_found_peptides"), 1)
    }
    END_SECTION

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID without shared peps)
    {
      vector<ProteinIdentification> prots;
      vector<PeptideIdentification> peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BasicProteinInferenceAlgorithm bpia;
      Param p = bpia.getParameters();
      p.setValue("use_shared_peptides","false");
      p.setValue("min_peptides_per_protein", 0);
      bpia.setParameters(p);
      bpia.run(peps, prots);
      TEST_EQUAL(prots[0].getHits()[0].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[1].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[2].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[3].getScore(), 0.8)
      TEST_EQUAL(prots[0].getHits()[4].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[5].getScore(), 0.9)

      TEST_EQUAL(prots[0].getHits()[0].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[1].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[2].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[3].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[4].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[5].getMetaValue("nr_found_peptides"), 1)
    }
    END_SECTION

END_TEST
