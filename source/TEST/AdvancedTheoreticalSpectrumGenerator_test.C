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
#include <OpenMS/FORMAT/MzDataFile.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

///////////////////////////

START_TEST(AdvancedTheoreticalSpectrumGenerator, "$Id$")

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

START_SECTION(void simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Int charge=1))
  // init rng
  gsl_rng* rnd_gen = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rnd_gen, 0);
  RichPeakSpectrum spec;
  ptr->loadProbabilisticModel();
  ptr->simulate(spec, peptide,rnd_gen);
  gsl_rng_free(rnd_gen);

  MSExperiment<RichPeak1D>exp;
  //exp.push_back(spec);
  MzDataFile mz_data_file;
  //mz_data_file.store("tmp_sandro.mzData", exp);

  mz_data_file.load(OPENMS_GET_TEST_DATA_PATH("AdvancedTheoreticalSpectrumGenerator_test.mzData"),exp);

  TEST_EQUAL(exp.size(), 1);
  if(exp.size())
  {
    TEST_EQUAL(spec.size(), exp[0].size());
    Size min_size = min(spec.size(), exp[0].size());

    for(Size i = 0; i<min_size; ++i)
    {
      TEST_REAL_SIMILAR(spec[i].getPosition()[0],(exp[0][i]).getPosition()[0]);
      TEST_EQUAL(spec[i].getIntensity(),(exp[0][i]).getIntensity());
      }
  }
END_SECTION

START_SECTION(void loadProbabilisticModel())
NOT_TESTABLE
END_SECTION



START_SECTION([EXTRA]UInt IndexConverter::operator(const UInt &type_id_a, const UInt &intensity_level_a, const UInt &intensity_level_parent, const UInt &number_intensity_levels))
  AdvancedTheoreticalSpectrumGenerator::IndexConverter ind_conv;
  TEST_EQUAL(ind_conv(10,3,2,5), 263)
END_SECTION

AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork * tan_ptr = 0;

START_SECTION([EXTRA]AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork())
tan_ptr = new AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork();
TEST_NOT_EQUAL(tan_ptr, 0)
END_SECTION

START_SECTION([EXTRA]AdvancedTheoreticalSpectrumGenerator::~TreeAugmentedNetwork())
  delete tan_ptr;
END_SECTION

typedef AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork::TanEdge TanEdge;
typedef std::vector<TanEdge>EdgeVector;

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

std::vector<Int>has_parent;
std::vector<UInt>dfs_order;

START_SECTION([EXTRA]AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork void generateTree(std::vector<Int> &tree_structure))
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

START_SECTION([EXTRA]AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork(AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork & rhs))
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

START_SECTION([EXTRA]AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork operator =(const AdvancedTheoreticalSpectrumGenerator::TreeAugmentedNetwork & rhs))
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
