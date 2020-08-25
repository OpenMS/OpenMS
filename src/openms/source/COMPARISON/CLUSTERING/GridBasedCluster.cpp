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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>

namespace OpenMS
{

GridBasedCluster::GridBasedCluster(const Point &centre, const Rectangle &bounding_box, const std::vector<int> &point_indices, const int &property_A, const std::vector<int> &properties_B)
: centre_(centre), bounding_box_(bounding_box), point_indices_(point_indices), property_A_(property_A), properties_B_(properties_B)
{
}

GridBasedCluster::GridBasedCluster(const Point &centre, const Rectangle &bounding_box, const std::vector<int> &point_indices)
: centre_(centre), bounding_box_(bounding_box), point_indices_(point_indices), property_A_(-1), properties_B_(point_indices.size(),-1)
{
}

const GridBasedCluster::Point& GridBasedCluster::getCentre() const
{
  return centre_;
}

const GridBasedCluster::Rectangle& GridBasedCluster::getBoundingBox() const
{
  return bounding_box_;
}

const std::vector<int>& GridBasedCluster::getPoints() const
{
  return point_indices_;
}

int GridBasedCluster::getPropertyA() const
{
  return property_A_;
}

const std::vector<int>& GridBasedCluster::getPropertiesB() const
{
  return properties_B_;
}

bool GridBasedCluster::operator<(const GridBasedCluster& other) const
{
  return centre_.getY() < other.centre_.getY();
}

bool GridBasedCluster::operator>(const GridBasedCluster& other) const
{
  return centre_.getY() > other.centre_.getY();
}

bool GridBasedCluster::operator==(const GridBasedCluster& other) const
{
  return centre_.getY() == other.centre_.getY();
}

}
