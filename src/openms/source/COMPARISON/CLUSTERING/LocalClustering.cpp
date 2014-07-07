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
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/HashGrid2.h>
#include <OpenMS/COMPARISON/CLUSTERING/LocalClustering.h>
#include <OpenMS/COMPARISON/CLUSTERING/Cluster.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>

namespace OpenMS
{
    
LocalClustering::LocalClustering(const std::vector<double> &data_x, const std::vector<double> &data_y, const std::vector<int> &properties_A, const std::vector<int> &properties_B, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y, double scaling_y)
: grid_(grid_spacing_x,grid_spacing_y), scaling_y_(scaling_y)
{
    init(data_x, data_y, properties_A, properties_B);
}

LocalClustering::LocalClustering(const std::vector<double> &data_x, const std::vector<double> &data_y, std::vector<double> grid_spacing_x, std::vector<double> grid_spacing_y, double scaling_y)
: grid_(grid_spacing_x,grid_spacing_y), scaling_y_(scaling_y)
{
    // set properties A and B to -1, i.e. ignore properties when clustering
    init(data_x, data_y, minusOnes(data_x.size()), minusOnes(data_x.size()));
}

void LocalClustering::init(const std::vector<double> &data_x, const std::vector<double> &data_y, const std::vector<int> &properties_A, const std::vector<int> &properties_B)
{
    // fill the grid with points to be clustered (initially each cluster contains a single point)
    for (unsigned i = 0; i < data_x.size(); ++i)
    {
        //std::cout << "i = " << i << "    cell count = " << grid_.getCellCount() << "\n";
        
        Point position(data_x[i],data_y[i]);
        Rectangle box(position,position);
        
        std::vector<int> pi;    // point indicies
        pi.push_back(i);
        std::vector<int> pb;    // properties B
        pb.push_back(properties_B[i]);
        
        // add to cluster list
        Cluster cluster(position, box, pi, properties_A[i], pb);
        clusters_.insert(std::make_pair(i,cluster));
        
        // register on hash grid
        //std::cout << "(x,y) = (" << position.getX() << ", " << position.getY() << ")    (i,j) = (" << grid_.getIndex(position).first << ", " << grid_.getIndex(position).second << ")\n";
        grid_.addCluster(grid_.getIndex(position), i);
    }
    //std::cout << "\n";
    
    // fill list of minimum distances
    std::map<int, Cluster>::iterator iterator = clusters_.begin();
    while (iterator != clusters_.end())
    {
        int cluster_index = iterator->first;
        Cluster cluster = iterator->second;
        
        if (findNearestNeighbour(cluster, cluster_index))
        {
            // remove from grid
            grid_.removeCluster(grid_.getIndex(cluster.getCentre()), cluster_index);
            // remove from cluster list
            clusters_.erase(iterator++);
        }
        else
        {
            ++iterator;
        }
    }
    
}

void LocalClustering::cluster()
{
    MinimumDistance zero_distance(-1, -1, 0);
    typedef std::multiset<MinimumDistance, std::less<vector<MinimumDistance>::value_type> >::iterator MultisetIterator;
    
    // combine clusters until all have been moved to the final list
    while (clusters_.size() > 0)
    {
        //std::cout << "size of cluster list = " << clusters_.size() << "\n";
        
        MultisetIterator smallest_distance_it = distances_.lower_bound(zero_distance);
        MinimumDistance smallest_distance(*smallest_distance_it);
        distances_.erase(smallest_distance_it);
        int cluster_index1 = smallest_distance.getClusterIndex();
        int cluster_index2 = smallest_distance.getNearestNeighbourIndex();
        
        // update cluster list
        Cluster cluster1 = clusters_.find(cluster_index1)->second;
        Cluster cluster2 = clusters_.find(cluster_index2)->second;
        std::vector<int> points1 = cluster1.getPoints();
        std::vector<int> points2 = cluster2.getPoints();
        std::vector<int> new_points;
        new_points.reserve(points1.size() + points2.size());
        new_points.insert(new_points.end(), points1.begin(), points1.end());
        new_points.insert(new_points.end(), points2.begin(), points2.end());
        
        double new_x = (cluster1.getCentre().getX() * points1.size() + cluster2.getCentre().getX() * points2.size()) / (points1.size() + points2.size());
        double new_y = (cluster1.getCentre().getY() * points1.size() + cluster2.getCentre().getY() * points2.size()) / (points1.size() + points2.size());
    
        Rectangle box1 = cluster1.getBoundingBox();
        Rectangle box2 = cluster2.getBoundingBox();
        Rectangle new_box(box1);
        new_box.enlarge(box2.minPosition());
        new_box.enlarge(box2.maxPosition());
        
        // Properties A of both clusters should by now be the same. The merge veto has been checked
        // when a new entry to the minimum distance list was added, @see findNearestNeighbour.
        if (cluster1.getPropertyA() != cluster2.getPropertyA())
        {
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Property A of both clusters not the same. ", "A");
        }
        int new_A = cluster1.getPropertyA();
        
        std::vector<int> B1 = cluster1.getPropertiesB();
        std::vector<int> B2 = cluster2.getPropertiesB();
        std::vector<int> new_B;
        new_B.reserve(B1.size() + B2.size());
        new_B.insert(new_B.end(), B1.begin(), B1.end());
        new_B.insert(new_B.end(), B2.begin(), B2.end());
        
        Cluster new_cluster(DPosition<2>(new_x,new_y), new_box, new_points, new_A, new_B);
        
        clusters_.erase(clusters_.find(cluster_index1));
        clusters_.erase(clusters_.find(cluster_index2));
        clusters_.insert(std::make_pair(cluster_index1, new_cluster));
        
        // update hash grid
        CellIndex cell_for_cluster1 = grid_.getIndex(cluster1.getCentre());
        CellIndex cell_for_cluster2 = grid_.getIndex(cluster2.getCentre());
        CellIndex cell_for_new_cluster = grid_.getIndex(DPosition<2>(new_x,new_y));
        
        grid_.removeCluster(cell_for_cluster1, cluster_index1);
        grid_.removeCluster(cell_for_cluster2, cluster_index2);
        grid_.addCluster(cell_for_new_cluster, cluster_index1);
    
        // update minimum distance list
        std::vector<int> clusters_to_be_updated;
        clusters_to_be_updated.push_back(cluster_index1);    // index of the new cluster
        
        MultisetIterator it = distances_.begin();
        while (it != distances_.end())
        {
            MinimumDistance distance(*it);
            if (distance.getClusterIndex() == cluster_index2)
            {
                // distance.getClusterIndex() == cluster_index1 should never happen since this entry has already been erased.
                // Cluster 2 no longer exists. Hence only remove the entry from the list of minimum distances.
                distances_.erase(it++);
            }
            else if (distance.getNearestNeighbourIndex() == cluster_index1 || distance.getNearestNeighbourIndex() == cluster_index2)
            {
                clusters_to_be_updated.push_back(distance.getClusterIndex());
                distances_.erase(it++);
            }
            else
            {
                ++it;
            }
        }
        
        for (std::vector<int>::const_iterator cluster_index = clusters_to_be_updated.begin(); cluster_index != clusters_to_be_updated.end(); ++cluster_index)
        {
            if (findNearestNeighbour(clusters_.find(*cluster_index)->second,*cluster_index))
            {
                Cluster c = clusters_.find(*cluster_index)->second;
                grid_.removeCluster(grid_.getIndex(c.getCentre()), *cluster_index);    // remove from grid
                clusters_.erase(clusters_.find(*cluster_index));    // remove from cluster list
           }
        }
                 
    }
}

void LocalClustering::extendClustersY()
{
    std::cout << "Extending clusters in y-direction.\n";
    
    // construct new grid (grid only in x-direction, single cell in y-direction)
    std::vector<double> grid_spacing_x = grid_.getGridSpacingX();
    std::vector<double> grid_spacing_y = grid_.getGridSpacingY();
    std::vector<double> grid_spacing_y_new;
    grid_spacing_y_new.push_back(grid_spacing_y.front());
    grid_spacing_y_new.push_back(grid_spacing_y.back());
    HashGrid2 grid_x_only(grid_spacing_x, grid_spacing_y_new);
    
    std::cout << "size x = " << grid_spacing_x.size() << "    size y = " << grid_spacing_y_new.size() << "\n";
    
    // register final clusters on the new hash grid
    for (std::map<int, Cluster>::iterator it = clusters_final_.begin(); it != clusters_final_.end(); ++it)
    {
        int cluster_index = it->first;
        Cluster cluster = it->second;
        grid_x_only.addCluster(grid_x_only.getIndex(cluster.getCentre()), cluster_index);
        
        std::cout << "    cell index = " << grid_x_only.getIndex(cluster.getCentre()).first << "  " << grid_x_only.getIndex(cluster.getCentre()).second << "\n";
        std::cout << "    cell non-empty = " << grid_x_only.isNonEmptyCell(grid_x_only.getIndex(cluster.getCentre())) << "\n";
    }
    
    // scan through x on the grid
    for (unsigned cell = 0; cell < grid_spacing_x.size(); ++cell)
    {
        CellIndex grid_index(cell,0);
        if (grid_x_only.isNonEmptyCell(grid_index))
        {
            std::cout << "    non-empty cell = " << cell << "\n";
        
            std::list<int> cluster_indices = grid_x_only.getClusters(grid_index);    // indices of clusters in this x-range
            if (cluster_indices.size() > 1)
            {
                // First, order the clusters in ascending y.
                std::list<Cluster> cluster_list;    // list to order clusters in y
                std::map<Cluster,int> index_list;    // allows us to keep track of cluster indices after sorting
                for (std::list<int>::iterator it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
                {
                    cluster_list.push_back(clusters_final_.find(*it)->second);
                    index_list.insert(std::make_pair(clusters_final_.find(*it)->second,*it));
                }
                cluster_list.sort();
                
                // Now check if two adjacent clusters c1 and c2 can be merged.
                std::list<Cluster>::iterator c1 = cluster_list.begin();
                std::list<Cluster>::iterator c2 = cluster_list.begin();
                ++c1;
                while (c1 != cluster_list.end())
                {
                    double centre1x = (*c1).getCentre().getX();
                    double centre1y = (*c1).getCentre().getY();
                    double centre2x = (*c2).getCentre().getX();
                    double centre2y = (*c2).getCentre().getY();
                    
                    double box1x_min = (*c1).getBoundingBox().minX();
                    double box1x_max = (*c1).getBoundingBox().maxX();
                    double box1y_min = (*c1).getBoundingBox().minY();
                    double box1y_max = (*c1).getBoundingBox().maxY();
                    double box2x_min = (*c2).getBoundingBox().minX();
                    double box2x_max = (*c2).getBoundingBox().maxX();
                    double box2y_min = (*c2).getBoundingBox().minY();
                    double box2y_max = (*c2).getBoundingBox().maxY();
                    
                    double y_range1 = box1y_max - box1y_min;
                    double y_range2 = box2y_max - box2y_min;
                    double y_gap = box1y_min - box2y_max;
                   
                    // Is the x-centre of one cluster in the x-range of the other?
                    bool centre_in_range1 = (box2x_min <= centre1x && centre1x <= box2x_max);
                    bool centre_in_range2 = (box1x_min <= centre2x && centre2x <= box1x_max);
                        
                    // Is the y-gap between the two clusters smaller than 1/s of their average y-range?
                    double s = 6;    // scaling factor
                    bool clusters_close = (y_gap * s <= (y_range1 - y_range2)/2);
                       
                    // Shall we merge the two adjacent clusters?
                    if ((centre_in_range1 || centre_in_range2) && clusters_close)
                    {
                        std::vector<int> points1 = (*c1).getPoints();
                        std::vector<int> points2 = (*c2).getPoints();
                        std::vector<int> new_points;
                        new_points.reserve(points1.size() + points2.size());
                        new_points.insert(new_points.end(), points1.begin(), points1.end());
                        new_points.insert(new_points.end(), points2.begin(), points2.end());
                        
                        double new_x = (centre1x * points1.size() + centre2x * points2.size()) / (points1.size() + points2.size());
                        double new_y = (centre1y * points1.size() + centre2y * points2.size()) / (points1.size() + points2.size());
                        
                        double min_x = std::min(box1x_min, box2x_min); 
                        double min_y = std::min(box1y_min, box2y_min); 
                        double max_x = std::max(box1x_max, box2x_max); 
                        double max_y = std::max(box1y_max, box2y_max);
                        
                        Point new_centre(new_x, new_y);
                        Point position_min(min_x,min_y);
                        Point position_max(max_x,max_y);
                        Rectangle new_bounding_box(position_min,position_max);
                        
                        Cluster new_cluster(new_centre, new_bounding_box, new_points);
                        
                        // update final cluster list
                        clusters_final_.erase(clusters_final_.find(index_list.find(*c1)->second));
                        clusters_final_.erase(clusters_final_.find(index_list.find(*c2)->second));
                        clusters_final_.insert(std::make_pair(index_list.find(*c1)->second, new_cluster));
                        
                        // update hash grid
                        CellIndex cell_for_cluster1 = grid_x_only.getIndex((*c1).getCentre());
                        CellIndex cell_for_cluster2 = grid_x_only.getIndex((*c2).getCentre());
                        CellIndex cell_for_new_cluster = grid_x_only.getIndex(new_centre);
                        
                        grid_x_only.removeCluster(cell_for_cluster1, index_list.find(*c1)->second);
                        grid_x_only.removeCluster(cell_for_cluster2, index_list.find(*c2)->second);
                        grid_x_only.addCluster(cell_for_new_cluster, index_list.find(*c1)->second);
                    }
                    else
                    {
                        ++c1;
                        ++c2;
                    }

                } 
            }
        }
    }
        
}

void LocalClustering::removeSmallClustersY(double threshold_y)
{
    std::map<int, Cluster>::iterator it = clusters_final_.begin();
    while (it != clusters_final_.end())
    {
        Rectangle box = it->second.getBoundingBox();
        if (box.maxY() - box.minY() < threshold_y)
        {
            clusters_final_.erase(it++);
        }
        else
        {
            ++it;
        }
    }    
}

std::map<int, Cluster> LocalClustering::getResults() const
{
    //std::cout << "number of final clusters = " << clusters_final_.size() << "\n";
    return clusters_final_;
}

LocalClustering::MinimumDistance::MinimumDistance(const int &cluster_index, const int &nearest_neighbour_index, const double &distance)
:cluster_index_(cluster_index), nearest_neighbour_index_(nearest_neighbour_index), distance_(distance)
{
}

int LocalClustering::MinimumDistance::getClusterIndex() const
{
    return cluster_index_;
}

int LocalClustering::MinimumDistance::getNearestNeighbourIndex() const
{
    return  nearest_neighbour_index_;
}

bool LocalClustering::MinimumDistance::operator<(MinimumDistance other) const
{
    return distance_ < other.distance_;
}

bool LocalClustering::MinimumDistance::operator>(MinimumDistance other) const
{
    return distance_ > other.distance_;
}

bool LocalClustering::MinimumDistance::operator==(MinimumDistance other) const
{
    return distance_ == other.distance_;
}

std::vector<int> LocalClustering::minusOnes(int l)
{
    std::vector<int> v;
    for (int i=0; i<l; ++i)
    {
        v.push_back(-1);
    }
    
    return v;
}

double LocalClustering::clusterDistance(Point c1, Point c2, double scaling) const
{
    return sqrt((c1.getX() - c2.getX())*(c1.getX() - c2.getX()) + scaling * scaling * (c1.getY() - c2.getY())*(c1.getY() - c2.getY()));
}

bool LocalClustering::mergeVeto(Cluster c1, Cluster c2) const
{
    int A1 = c1.getPropertyA();
    int A2 = c2.getPropertyA();
    std::vector<int> B1 = c1.getPropertiesB();
    std::vector<int> B2 = c2.getPropertiesB();
    
    // check if any of the properties A and B is not set i.e. =-1
    if (A1 == -1 || A2 == -1 || std::find(B1.begin(), B1.end(), -1) != B1.end() || std::find(B2.begin(), B2.end(), -1) != B2.end())
    {
        return false;
    }
    
    // Will the merged cluster have the same properties A?
    bool vetoA = !(A1 == A2);
    
    // Will the merged cluster have different properties B?
    std::vector<int> B;
    sort(B1.begin(), B1.end());
    sort(B2.begin(), B2.end());
    set_intersection(B1.begin(),B1.end(),B2.begin(),B2.end(),back_inserter(B));
    bool vetoB = !B.empty();
    
    return (vetoA || vetoB);
}

bool LocalClustering::findNearestNeighbour(Cluster cluster, int cluster_index)
{
    Point centre = cluster.getCentre();
    CellIndex cell_index = grid_.getIndex(centre);
    
    double min_dist = 0;
    int nearest_neighbour = -1;
    
    // search in the grid cell and its 8 neighbouring cells for the nearest neighbouring cluster
    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            CellIndex cell_index2(cell_index);
            cell_index2.first += i;
            cell_index2.second += j;
            if (grid_.isNonEmptyCell(cell_index2))
            {
                std::list<int> cluster_indices = grid_.getClusters(cell_index2);
                for (std::list<int>::const_iterator cluster_index2 = cluster_indices.begin(); cluster_index2 != cluster_indices.end(); ++cluster_index2) {
                    if (*cluster_index2 != cluster_index)
                    {
                        Cluster cluster2 = clusters_.find(*cluster_index2)->second;
                        Point centre2 = cluster2.getCentre();
                        double distance = clusterDistance(centre, centre2, scaling_y_);
                        bool veto = mergeVeto(cluster, cluster2);    // If clusters cannot be merged anyhow, they are no nearest neighbours.
                        if (!veto && (distance < min_dist || nearest_neighbour == -1))
                        {
                            min_dist = distance;
                            nearest_neighbour = *cluster_index2;
                        }
                    }
                }
            }
        }
    }
    
    if (nearest_neighbour == -1)
    {
        // no other cluster nearby, hence move the cluster to the final results
        clusters_final_.insert(std::make_pair(cluster_index, clusters_.find(cluster_index)->second));
        return true;
    }
    
    // add to the list of minimal distances
    distances_.insert(MinimumDistance(cluster_index, nearest_neighbour, min_dist));
    return false;
}

}
