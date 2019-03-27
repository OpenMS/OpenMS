// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow$
// $Authors: Dominik Schmitz, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/QC/Contaminants.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

//FASTAFile fasta_file;
//vector<FASTAFile::FASTAEntry> contaminants;
//fasta_file.load("/buffer/ag_bsc/pmsb/data/databases/crab.fasta",contaminants);
//FeatureXMLFile fxml_file;
//FeatureMap fmap;
//fxml_file.load("/buffer/ag_bsc/pmsb/data/Example_Data/lfq_spikein_dilution_1.featureXML", fmap);



START_TEST(Contaminants, "$Id$")

FASTAFile::FASTAEntry contaminantsProtein("test_protein", "protein consists of only Alanine or Cytosine", "AAAAAAAAAAKRAAAAAAAAAAKRCCCCCCCCCCKRCCCCCCCCCC");
vector<FASTAFile::FASTAEntry> contaminantsFile;
contaminantsFile.push_back(contaminantsProtein);


FeatureMap fmap;
Feature f;
fmap.getProteinIdentifications().resize(1);
fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("trypsin");
fmap.getProteinIdentifications()[0].getSearchParameters().missed_cleavages = 1;
PeptideIdentification id;
PeptideHit hit;

hit.setSequence(AASequence::fromString("AAAAAAAAAAK"));
id.setHits(std::vector<PeptideHit>(1, hit));
f.setPeptideIdentifications({id});
f.setIntensity(12.0);
fmap.push_back(f);

hit.setSequence(AASequence::fromString("RAAAAAAAAAAK"));
id.setHits(std::vector<PeptideHit>(1, hit));
f.setPeptideIdentifications({id});
f.setIntensity(8.0);
fmap.push_back(f);

hit.setSequence(AASequence::fromString("RAAAAAAAAAAK"));
id.setHits(std::vector<PeptideHit>(1, hit));
f.setPeptideIdentifications({id});
f.setIntensity(10.0);
fmap.push_back(f);

hit.setSequence(AASequence::fromString("RARKAAAAAAAK"));
id.setHits(std::vector<PeptideHit>(1, hit));
f.setPeptideIdentifications({id});
f.setIntensity(20.0);
fmap.push_back(f);

hit.setSequence(AASequence::fromString("QQQQQQQQQQ"));
id.setHits(std::vector<PeptideHit>(1, hit));
f.setPeptideIdentifications({id});
f.setIntensity(10.0);
fmap.push_back(f);

std::vector<PeptideIdentification> ids2(3);
PeptideHit hit2;
hit2.setSequence(AASequence::fromString("AAAAAAAAAAK"));
ids2[0].setHits(std::vector<PeptideHit>(1, hit));
hit2.setSequence(AASequence::fromString("RCCCCCCCCCCK"));
ids2[1].setHits(std::vector<PeptideHit>(1, hit));
hit2.setSequence(AASequence::fromString("DDDDDDDDDD"));
ids2[2].setHits(std::vector<PeptideHit>(1, hit));
fmap.setUnassignedPeptideIdentifications(ids2);


    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////

START_SECTION((void compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants)))
{
  Contaminants conts1, conts2;
  //conts1.compute(fmap, contaminantsFile);

  conts2.compute(fmap, contaminantsFile);
  //std::vector<std::pair<double, double>> result1 = conts1.getResults();
  std::vector<std::pair<double, double>> result2 = conts2.getResults();
  //ABORT_IF(result1.size() != 1);
  //TEST_REAL_SIMILAR(result1[0].first, 3/8.0);
  //TEST_REAL_SIMILAR(result1[0].second, 1/2.0);
  ABORT_IF(result2.size() != 1);
  TEST_REAL_SIMILAR(result2[0].first, 5/8.0);
  TEST_REAL_SIMILAR(result2[0].second, 1/2.0);
}
END_SECTION


    /////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////
END_TEST
