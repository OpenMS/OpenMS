// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

QTCluster* qtc_ptr = 0;
QTCluster* qtc_nullPointer = 0;

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

START_SECTION((QTCluster(GridFeature* center_point, Size num_maps, double max_distance, bool use_IDs)))
{
  qtc_ptr = new QTCluster(&gf, 2, 11.1, false);
  TEST_NOT_EQUAL(qtc_ptr, qtc_nullPointer);
}
END_SECTION

START_SECTION((~QTCluster()))
{
  delete qtc_ptr;
}
END_SECTION

QTCluster cluster(&gf, 2, 11.1, true);

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

START_SECTION((Size size() const))
{
  TEST_EQUAL(cluster.size(), 1);
}
END_SECTION

GridFeature gf2(bf, 789, 1011);

START_SECTION((void add(GridFeature* element, double distance)))
{
  cluster.add(&gf2, 3.3);
  TEST_EQUAL(cluster.size(), 2);
}
END_SECTION

START_SECTION((bool operator<(QTCluster& cluster)))
{
  QTCluster cluster2(&gf, 2, 11.1, false);
  TEST_EQUAL(cluster2 < cluster, true);
}
END_SECTION

START_SECTION((void getElements(boost::unordered::unordered_map<Size, GridFeature*>& elements)))
{
  boost::unordered::unordered_map<Size, GridFeature*> elements;
  cluster.getElements(elements);
  TEST_EQUAL(elements.size(), 2);
  TEST_EQUAL(elements[123], &gf);
  TEST_EQUAL(elements[789], &gf2);
}
END_SECTION

START_SECTION((bool update(const boost::unordered::unordered_map<Size, GridFeature*>& removed)))
{
  boost::unordered::unordered_map<Size, GridFeature*> removed;
  removed[789] = &gf2;
  TEST_EQUAL(cluster.update(removed), true);
  TEST_EQUAL(cluster.size(), 1);
  removed[123] = &gf;
  // removing the center invalidates the cluster:
  TEST_EQUAL(cluster.update(removed), false);
}
END_SECTION

START_SECTION((double getQuality()))
{
  TEST_EQUAL(cluster.getQuality(), 0.0);
  cluster.add(&gf2, 3.3);
  TEST_EQUAL(cluster.getQuality(), (11.1 - 3.3) / 11.1);
}
END_SECTION

START_SECTION((const set<AASequence>& getAnnotations()))
{
  TEST_EQUAL(cluster.getAnnotations().size(), 2);
  TEST_EQUAL(*(cluster.getAnnotations().begin()), AASequence::fromString("AAA"));
  TEST_EQUAL(*(cluster.getAnnotations().rbegin()), AASequence::fromString("CCC"));
  QTCluster cluster2(&gf, 2, 11.1, false);
  TEST_EQUAL(cluster2.getAnnotations().empty(), true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
