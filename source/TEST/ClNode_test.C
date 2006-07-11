// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  $
// --------------------------------------------------------------------------
//



#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>

///////////////////////////

///////////////////////////
START_TEST(ClusterNode, "$Id: ClNode_test.C,v 1.3 2006/03/29 12:30:29 andreas_bertsch Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ClusterNode* clp;

CHECK(ClusterNode::ClusterNode())
  clp = new ClusterNode();
  TEST_NOT_EQUAL(clp, 0)
RESULT

CHECK(ClusterNode::insert())
  ClusterNode* min = new ClusterNode(1,0.2);
  ClusterNode* max = new ClusterNode(2,2.2);
  clp->insert(min);
  clp->insert(max);
  TEST_REAL_EQUAL(clp->getMinParentMass(),0.2)
  TEST_REAL_EQUAL(clp->getMaxParentMass(),2.2)
  TEST_EQUAL(2,clp->size())
RESULT

CHECK(ClusterNode::ClusterNode(ClusterNode*,ClusterNode*))
  ClusterNode* clp2 = new ClusterNode(new ClusterNode(3,0.1),new ClusterNode(4,3));
  TEST_NOT_EQUAL(clp2, 0)
  TEST_REAL_EQUAL(clp2->getMinParentMass(),0.1)
  TEST_REAL_EQUAL(clp2->getMaxParentMass(),3)
  TEST_EQUAL(clp2->size(),2)
RESULT

CHECK(ClusterNode::ClusterNode(const ClusterNode&))
  ClusterNode cl(*clp);
  TEST_REAL_EQUAL(cl.getMinParentMass(),0.2)
  TEST_REAL_EQUAL(cl.getMaxParentMass(),2.2)
  TEST_EQUAL(2,cl.size())
RESULT

CHECK(ClusterNode::operator=(const ClusterNode&))
  ClusterNode cl = *clp;
  TEST_REAL_EQUAL(cl.getMinParentMass(),0.2)
  TEST_REAL_EQUAL(cl.getMaxParentMass(),2.2)
  TEST_EQUAL(2,cl.size())
RESULT

CHECK(ClusterNode::children())
  list<int>::const_iterator cit = clp->children().begin();
  TEST_EQUAL(*cit,1)
  ++cit;
  TEST_EQUAL(*cit,2)
  ++cit;
  bool end = (cit == clp->children().end());
  TEST_EQUAL(end,1)
RESULT

CHECK(ClusterNode::~ClusterNode())
  delete clp;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
