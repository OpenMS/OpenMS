// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedClustering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>

using namespace OpenMS;

START_TEST(GridBasedClustering, "$Id$")

MultiplexClustering::MultiplexDistance metric(1);

std::vector<double> grid_spacing_x;
std::vector<double> grid_spacing_y;
for (double i = 0; i <= 10; ++i)
{
    grid_spacing_x.push_back(i);
    grid_spacing_y.push_back(i);
}

std::vector<double> data_x;
std::vector<double> data_y;
std::vector<int> properties_A;
std::vector<int> properties_B;
for (int i = 0; i < 1000; ++i)
{
    data_x.push_back(5*(sin(static_cast<double>(i))+1));
    data_y.push_back(5*(sin(static_cast<double>(i+18))+1));
    properties_A.push_back(1);    // Should be the same within each cluster.
    properties_B.push_back(i);    // Should be different within each cluster.
}

GridBasedClustering<MultiplexClustering::MultiplexDistance>* nullPointer = nullptr;
GridBasedClustering<MultiplexClustering::MultiplexDistance>* ptr;

START_SECTION(GridBasedClustering(Metric metric, const std::vector<double> &data_x, const std::vector<double> &data_y, const std::vector<int> &properties_A, const std::vector<int> &properties_B, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y))
    GridBasedClustering<MultiplexClustering::MultiplexDistance> clustering(metric, data_x, data_y, properties_A, properties_B, grid_spacing_x, grid_spacing_y);
    clustering.cluster();
    TEST_EQUAL(clustering.getResults().size(), 12);
    ptr = new GridBasedClustering<MultiplexClustering::MultiplexDistance>(metric, data_x, data_y, properties_A, properties_B, grid_spacing_x, grid_spacing_y);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(GridBasedClustering(Metric metric, const std::vector<double> &data_x, const std::vector<double> &data_y, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y))
    GridBasedClustering<MultiplexClustering::MultiplexDistance> clustering(metric, data_x, data_y, grid_spacing_x, grid_spacing_y);
    clustering.cluster();
    TEST_EQUAL(clustering.getResults().size(), 12);
    ptr = new GridBasedClustering<MultiplexClustering::MultiplexDistance>(metric, data_x, data_y, grid_spacing_x, grid_spacing_y);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

GridBasedClustering<MultiplexClustering::MultiplexDistance> clustering(metric, data_x, data_y, grid_spacing_x, grid_spacing_y);

START_SECTION(void cluster())
    clustering.cluster();
    TEST_EQUAL(clustering.getResults().size(), 12);
END_SECTION

START_SECTION(std::map<int Cluster> getResults() const)
    clustering.cluster();
    TEST_EQUAL(clustering.getResults().size(), 12);
END_SECTION

START_SECTION(void extendClustersY())
    clustering.cluster();
    clustering.extendClustersY();
    TEST_EQUAL(clustering.getResults().size(), 11);
END_SECTION

START_SECTION(void removeSmallClustersY(double threshold_y))
    clustering.cluster();
    clustering.removeSmallClustersY(0.9);
    TEST_EQUAL(clustering.getResults().size(), 8);
END_SECTION

END_TEST
