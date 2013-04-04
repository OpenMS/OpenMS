// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace OpenMS;
using namespace std;


class LowLevelComparator
{
public:
  double operator()(const double first, const double second) const
  {
    double x, y;
    x = min(second, first);
    y = max(first, second);
    if ((y - x) > 1)
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    return 1 - (y - x);
  }

}; // end of LowLevelComparator

Int main()
{
  // data
  vector<double> data; // must be filled

  srand(333);
  for (int i = 0; i < 12; ++i)
  {
    data.push_back((double)rand() / RAND_MAX);
  }

  LowLevelComparator llc;
  CompleteLinkage sl;
  vector<BinaryTreeNode> tree;
  DistanceMatrix<Real> dist; // will be filled
  ClusterHierarchical ch;
  ch.setThreshold(0.15);

  // clustering
  ch.cluster<double, LowLevelComparator>(data, llc, sl, tree, dist);

  ClusterAnalyzer ca;
  std::cout << ca.newickTree(tree) << std::endl;

  return 0;
} //end of main
