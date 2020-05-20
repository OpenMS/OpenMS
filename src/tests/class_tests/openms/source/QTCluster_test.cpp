// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
///////////////////////////

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/BaseFeature.h>

START_TEST(QTCluster, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

QTCluster* qtc_ptr = nullptr;
QTCluster* qtc_nullPointer = nullptr;
QTCluster::BulkData* qtc_data_ptr = nullptr;

BaseFeature bf;
bf.setRT(1.1);
bf.setMZ(2.2);
bf.setCharge(3);
bf.getPeptideIdentifications().resize(2);
PeptideHit hit;
hit.setSequence(AASequence::fromString("AAA"));
bf.getPeptideIdentifications()[0].insertHit(hit);
hit.setSequence(AASequence::fromString("CCC"));
bf.getPeptideIdentifications()[1].insertHit(hit);
GridFeature gf(bf, 123, 456);

START_SECTION((QTCluster::BulkData(const OpenMS::GridFeature* const center_point, Size num_maps, double max_distance, Int x_coord, Int y_coord, Size id)))
{
  qtc_data_ptr = new QTCluster::BulkData(&gf, 2, 11.1, 0, 0, 0);
  TEST_NOT_EQUAL(qtc_data_ptr, qtc_nullPointer);
}
END_SECTION

START_SECTION((QTCluster(BulkData* const data, bool use_IDs)))
{
  qtc_ptr = new QTCluster(qtc_data_ptr, false);
  TEST_NOT_EQUAL(qtc_ptr, qtc_nullPointer);
}
END_SECTION

START_SECTION((~QTCluster()))
{
  delete qtc_ptr;
}
END_SECTION

START_SECTION((QTCluster::~BulkData()))
{
  delete qtc_data_ptr;
}
END_SECTION

QTCluster::BulkData qtc_data(&gf, 2, 11.1, 7, 9, 1);
QTCluster cluster(&qtc_data, true);

START_SECTION((double getCenterRT() const))
{
  TEST_EQUAL(cluster.getCenterRT(), 1.1);
}
END_SECTION

START_SECTION((double getCenterMZ() const))
{
  TEST_EQUAL(cluster.getCenterMZ(), 2.2);
}
END_SECTION

START_SECTION((Int getXCoord() const))
{
  TEST_EQUAL(cluster.getXCoord(), 7);
}
END_SECTION

START_SECTION((Int getYCoord() const))
{
  TEST_EQUAL(cluster.getYCoord(), 9);
}
END_SECTION

START_SECTION((Size getId() const))
{
  TEST_EQUAL(cluster.getId(), 1);
}
END_SECTION

START_SECTION((Size size() const))
{
  TEST_EQUAL(cluster.size(), 1);
}
END_SECTION

GridFeature gf2(bf, 789, 1011);

START_SECTION((void add(const GridFeature* const element, double distance)))
{
  cluster.initializeCluster();
  cluster.add(&gf2, 3.3);
  cluster.finalizeCluster();
  TEST_EQUAL(cluster.size(), 2);
}
END_SECTION

START_SECTION((bool operator<(QTCluster& cluster)))
{
  QTCluster::BulkData data(&gf, 2, 11.1, 0, 0, 2);
  QTCluster cluster2(&data, false);
  TEST_EQUAL(cluster2 < cluster, true);
}
END_SECTION

START_SECTION((QTCuster::Elements getElements() const))
{
  QTCluster::Elements elements = cluster.getElements();
  TEST_EQUAL(elements.size(), 2);

  if (elements[0].feature != &gf)
  {
    TEST_EQUAL(elements[0].feature, &gf2);
    TEST_EQUAL(elements[1].feature, &gf);
  }
  else
  {
    TEST_EQUAL(elements[0].feature, &gf);
    TEST_EQUAL(elements[1].feature, &gf2);
  }
}
END_SECTION

START_SECTION((QTCluster::Elements getAllNeighbors() const))
{
  GridFeature gf3(bf, 789, 1012);
  GridFeature gf4(bf, 222, 1011);

  QTCluster::BulkData data(&gf, 2, 11.1, 0, 0, 2);
  QTCluster cluster2(&data, false);
  TEST_EQUAL(cluster2.getAllNeighbors().size(), 0)
  cluster2.initializeCluster();
  cluster2.add(&gf2, 3.3);
  cluster2.finalizeCluster();
  TEST_EQUAL(cluster2.getAllNeighbors().size(), 1)
  TEST_EQUAL(cluster2.getAllNeighbors()[0].feature, &gf2)

  // adding a better feature from the same map does not increase neighbor size
  cluster2.initializeCluster();
  cluster2.add(&gf3, 3.0);
  cluster2.finalizeCluster();
  TEST_EQUAL(cluster2.getAllNeighbors().size(), 1)
  TEST_EQUAL(cluster2.getAllNeighbors()[0].feature, &gf3)

  // adding features from a new map will increase neighbor size
  cluster2.initializeCluster();
  cluster2.add(&gf4, 3.9);
  cluster2.add(&gf4, 3.2);
  cluster2.add(&gf4, 3.1);
  cluster2.add(&gf4, 3.8);
  cluster2.finalizeCluster();
  
  QTCluster::Elements neighbors = cluster2.getAllNeighbors();

  TEST_EQUAL(neighbors.size(), 2)
  if (neighbors[0].feature != &gf3)
  {
    TEST_EQUAL(neighbors[0].feature, &gf4);
    TEST_EQUAL(neighbors[1].feature, &gf3);
  }
  else
  {
    TEST_EQUAL(neighbors[0].feature, &gf3);
    TEST_EQUAL(neighbors[1].feature, &gf4);
  }
}
END_SECTION

START_SECTION(bool update(const QTCluster::Elements& removed))
{
  QTCluster::Elements removed;
  removed.push_back({789, &gf2});
  TEST_EQUAL(cluster.update(removed), true);
  TEST_EQUAL(cluster.size(), 1);
  removed.push_back({123, &gf});

  // removing the center invalidates the cluster:
  TEST_EQUAL(cluster.update(removed), false);
  TEST_EQUAL(cluster.isInvalid(), true);
}
END_SECTION

QTCluster::BulkData qtc_data2(&gf, 2, 11.1, 7, 9, 3);

START_SECTION((double getQuality()))
{
  // cluster is invalid, we shouldnt use it any more -> create a new one
  TEST_EQUAL(cluster.isInvalid(), true);

  cluster = QTCluster(&qtc_data2, true);

  cluster.initializeCluster();
  cluster.add(&gf2, 3.3);
  cluster.finalizeCluster();
  TEST_EQUAL(cluster.getQuality(), (11.1 - 3.3) / 11.1);
  TEST_EQUAL(cluster.isInvalid(), false);
}
END_SECTION

START_SECTION(double getCurrentQuality() const)
{
  TEST_EQUAL(cluster.getCurrentQuality(), cluster.getQuality());
}
END_SECTION

START_SECTION((const set<AASequence>& getAnnotations()))
{
  TEST_EQUAL(cluster.getAnnotations().size(), 2);
  TEST_EQUAL(*(cluster.getAnnotations().begin()), AASequence::fromString("AAA"));
  TEST_EQUAL(*(cluster.getAnnotations().rbegin()), AASequence::fromString("CCC"));
  QTCluster::BulkData data(&gf, 2, 11.1, 0, 0, 2);
  QTCluster cluster2(&data, false);
  TEST_EQUAL(cluster2.getAnnotations().empty(), true);
}
END_SECTION

START_SECTION((inline bool isInvalid() const))
{
  TEST_EQUAL(cluster.isInvalid(), false);
}
END_SECTION

START_SECTION((void setInvalid()))
{
  cluster.setInvalid();
  TEST_EQUAL(cluster.isInvalid(), true);
}
END_SECTION

START_SECTION((void finalizeCluster() ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((void initializeCluster() ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
