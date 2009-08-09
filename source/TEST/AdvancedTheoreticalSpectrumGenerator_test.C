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
// $Maintainer: Sandro Andreotti $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/CHEMISTRY/AdvancedTheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(AdvancedTheoreticalSpectrumGenerator, "$Id: TheoreticalSpectrumGenerator_test.C 4776 2009-03-05 14:14:35Z groepl $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

AdvancedTheoreticalSpectrumGenerator* ptr = 0;

START_SECTION(AdvancedTheoreticalSpectrumGenerator())
  ptr = new AdvancedTheoreticalSpectrumGenerator();
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(AdvancedTheoreticalSpectrumGenerator(const AdvancedTheoreticalSpectrumGenerator& source))
  AdvancedTheoreticalSpectrumGenerator copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~AdvancedTheoreticalSpectrumGenerator())
  delete ptr;
END_SECTION

ptr = new AdvancedTheoreticalSpectrumGenerator();
AASequence peptide("IFSQVGK");

START_SECTION(AdvancedTheoreticalSpectrumGenerator& operator = (const AdvancedTheoreticalSpectrumGenerator& tsg))
  AdvancedTheoreticalSpectrumGenerator copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1))
//TODO when one probabilistic model is fix
//This also tests the method: void loadProbabilisticModel();

/*RichPeakSpectrum spec;
  ptr->getSpectrum(spec, peptide, 1);
  TEST_EQUAL(spec.size(), 12)

  TOLERANCE_ABSOLUTE(0.001)

  double result[] = {115.1, 147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear();
  ptr->getSpectrum(spec, peptide, 2);
  TEST_EQUAL(spec.size(), 24)
  */
END_SECTION

START_SECTION(UInt IndexConverter::operator(const UInt &type_id_a, const UInt &intensity_level_a, const UInt &intensity_level_parent, const UInt &number_intensity_levels))
  AdvancedTheoreticalSpectrumGenerator::IndexConverter ind_conv;
  TEST_EQUAL(ind_conv(10,3,2,5), 263)
END_SECTION

AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork * tan_ptr = 0;

START_SECTION(AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork())
tan_ptr = new AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork();
TEST_NOT_EQUAL(tan_ptr, 0)
END_SECTION

START_SECTION(AdvancedTheoreticalSpectrumGenerator::~TreeAugmentedNetwork())
  delete tan_ptr;
END_SECTION

typedef AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork::TanEdge TanEdge;
typedef std::vector<TanEdge>EdgeVector;

START_SECTION(AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork(EdgeVector & edges))
EdgeVector edges;
TanEdge e1={1,2,-2.0};
TanEdge e2={1,3,-5.0};
TanEdge e3={1,4,-6.0};
TanEdge e4={2,4,-3.0};
TanEdge e5={2,3,-7.0};
TanEdge e6={3,4,-4.0};
edges.push_back(e1);
edges.push_back(e2);
edges.push_back(e3);
edges.push_back(e4);
edges.push_back(e5);
edges.push_back(e6);

tan_ptr = new AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork(edges);
END_SECTION

std::vector<Int>has_parent;
std::vector<UInt>dfs_order;

START_SECTION(AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork void generateTree(std::vector<Int> &tree_structure))
tan_ptr->generateTree(has_parent);
TEST_EQUAL(has_parent.size(), 5)

TEST_EQUAL(has_parent[0],-1)
TEST_EQUAL(has_parent[1],-1)
TEST_EQUAL(has_parent[2], 3)
TEST_EQUAL(has_parent[3], 1)
TEST_EQUAL(has_parent[4], 1)

tan_ptr->getDFSOrder(dfs_order);

TEST_EQUAL(dfs_order.size(), 4)
TEST_EQUAL(dfs_order[0], 1)
TEST_EQUAL(dfs_order[1], 4)
TEST_EQUAL(dfs_order[2], 3)
TEST_EQUAL(dfs_order[3], 2)
END_SECTION

START_SECTION(AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork(AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork & rhs))
AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork copy(*tan_ptr);
std::vector<Int>copy_has_parent;
std::vector<UInt>copy_dfs_order;
copy.generateTree(copy_has_parent);
copy.getDFSOrder(copy_dfs_order);

TEST_EQUAL(copy_has_parent.size(), has_parent.size());
TEST_EQUAL(copy_dfs_order.size(), dfs_order.size());

for(Size i=0; i<copy_has_parent.size();++i)
  TEST_EQUAL(copy_has_parent[i], has_parent[i]);

for(Size i=0; i<copy_dfs_order.size();++i)
  TEST_EQUAL(copy_dfs_order[i], dfs_order[i]);
END_SECTION

START_SECTION(AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork operator =(const AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork & rhs))
AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork copy;
copy=*tan_ptr;

std::vector<Int>copy_has_parent;
std::vector<UInt>copy_dfs_order;
copy.generateTree(copy_has_parent);
copy.getDFSOrder(copy_dfs_order);

TEST_EQUAL(copy_has_parent.size(), has_parent.size());
TEST_EQUAL(copy_dfs_order.size(), dfs_order.size());

for(Size i=0; i<copy_has_parent.size();++i)
  TEST_EQUAL(copy_has_parent[i], has_parent[i]);

for(Size i=0; i<copy_dfs_order.size();++i)
  TEST_EQUAL(dfs_order[i], copy_dfs_order[i]);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
