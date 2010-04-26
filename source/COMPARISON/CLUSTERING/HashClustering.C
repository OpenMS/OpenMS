/*
 * HashClustering.C
 *
 *  Created on: 26.04.2010
 *      Author: steffen
 */

#include <OpenMS/COMPARISON/CLUSTERING/HashClustering.h>

typedef std::map<std::pair<int,int>, std::list<DataSubset*> > GridElements;

namespace OpenMS {
void HashClustering::init()
{
	min_distance=numeric_limits<double>::max();
	 hashing_map=grid.elements;

		for (GridElements::iterator it=hashing_map.begin();it!=hashing_map.end();++it)
		{
			std::pair<int,int> act_coords=it->first;
			std::list<GridElement*>& elements=it->second;
			int x=act_coords.first;
			int y=act_coords.second;
			for (int i=x-1;i<=x+1;++i)
			{
				if (i<0 || i>max_x)
					continue;
				for (int j=y;j<=y+1;++j)
				{
					if (j<0 || j>max_y || (i==x-1 && j==y))
						continue;
	//				cout << x <<"/"<< y << endl;
					GridElements::iterator actPos=hashing_map.find(make_pair(i,j));
					if (actPos==hashMap.end())
						continue;
					std::list<GridElement*>& neighbor_elements=actPos->second;
					for (std::list<GridElement*>::iterator lit=elements.begin();lit!=elements.end();++lit)
					{
						for (std::list<GridElement*>::iterator nit=neighbor_elements.begin();nit!=neighbor_elements.end();++nit)
						{
								if ((i==x && j==y && (*lit)->getID()>(*nit)->getID()) || **lit==**nit)
									continue;
								DataSubset* element_ptr = dynamic_cast<DataSubset*> (*lit);
								DataSubset* neighbor_ptr = dynamic_cast<DataSubset*> (*nit);
								double act_distance=method.getDistance(*element_ptr,*neighbor_ptr);
								std::cout << (*lit)->getID() << " " << (*nit)->getID() << " " << act_distance << std::endl;
								std::pair<DistanceSet::iterator,bool> position=distances.insert(DistanceEntry(element_ptr,neighbor_ptr,act_distance));
								if (position.second)
								{
									element_ptr->distance_iterators.insert(make_pair(neighbor_ptr,position.first));
									if (act_distance<min_distance)
									{
										min_distance=act_distance;
										min_distance_subsets=make_pair(element_ptr,neighbor_ptr);
									}
								}

						}
					}
				}
			}
}
void HashClustering::merge(DataSubset& subset1,DataSubset& subset2)
{

}
void HashClustering::updateMinElements()
{

}
vector<DataSubset*> HashClustering::performClustering()
{

}
}


