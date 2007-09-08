// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/HierarchicalClustering.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(HierarchicalClustering, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HierarchicalClustering<>* ptr = 0;
CHECK((HierarchicalClustering()))
        ptr = new HierarchicalClustering<>;
        TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~HierarchicalClustering()))
        delete ptr;
RESULT

CHECK((HierarchicalClustering(const HierarchicalClustering &source)))
  HierarchicalClustering<> hc;

  Param p = hc.getParameters();
  p.setValue("cluster_cutoff", 41.5);
  hc.setParameters(p);
  
  HierarchicalClustering<> hc_copy(hc);
  
  Param p_copy = hc_copy.getParameters();
  
  TEST_EQUAL((double) p_copy.getValue("cluster_cutoff"), 41.5);
  
RESULT

CHECK((HierarchicalClustering& operator=(const HierarchicalClustering &source)))
  HierarchicalClustering<> hc;

  Param p = hc.getParameters();
  p.setValue("cluster_cutoff", 41.5);
  hc.setParameters(p);
  
  HierarchicalClustering<> hc_copy = hc;
  
  Param p_copy = hc_copy.getParameters();
  
  TEST_EQUAL((double) p_copy.getValue("cluster_cutoff"), 41.5);
RESULT


CHECK((const ClusterIdxVectorType& getClusters() const))
  //create some data ...
  
  std::vector< HierarchicalClustering<>::ClusterPointType> points(4);
  points[0][0] = 1;
  points[0][1] = 1;
  points[1][0] = 3;
  points[1][1] = 1;
  points[2][0] = 2;
  points[2][1] = 2;
  points[3][0] = 4;
  points[3][1] = 4;
  
  HierarchicalClustering<> hc;
  Param p = hc.getParameters();
  p.setValue("cluster_cutoff", 1.5);  
  hc.setParameters(p);
  
  hc.compute(points);
  
  HierarchicalClustering<>::ClusterIdxVectorType cl = hc.getClusters();

  TEST_EQUAL(cl.size(),2);  // we should get 2 clusters
  TEST_EQUAL(cl[0].size(), 3); // ... the first one having three points
  TEST_EQUAL(cl[1].size(), 1); // ... the second with one point
  TEST_EQUAL(cl[1][0], 3); // ... , namely point #3
      
  hc.printStatistics(std::cout);
  
RESULT

CHECK((void compute(const std::vector< ClusterPointType > &points) throw (Exception::NotImplemented)))
  // all done above... no need to do it again; its hard to split those three
RESULT

CHECK((void printStatistics(std::ostream &os)))
  // all done above... no need to do it again; its hard to split those three
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
