// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BinaryTreeNode, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BinaryTreeNode* ptr = 0;
BinaryTreeNode* null_ptr = 0;
START_SECTION(BinaryTreeNode())
{
	ptr = new BinaryTreeNode();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~BinaryTreeNode())
{
	delete ptr;
}
END_SECTION

START_SECTION((BinaryTreeNode(const Size i, const Size j, const float x)))
{
  // TODO
}
END_SECTION

START_SECTION((~BinaryTreeNode()))
{
  // TODO
}
END_SECTION

START_SECTION((BinaryTreeNode(const BinaryTreeNode &source)))
{
  // TODO
}
END_SECTION

START_SECTION((BinaryTreeNode& operator=(const BinaryTreeNode &source)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



