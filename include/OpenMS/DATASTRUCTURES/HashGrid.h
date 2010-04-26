/*
 * GeometricHashing.h
 *
 *  Created on: 24.03.2010
 *      Author: steffen
 */

#include<map>
#include<list>
#include<set>
#include<vector>
#include<iostream>

#include <OpenMS/DATASTRUCTURES/GridElement.h>



#ifndef HASHGRID_H_
#define HASHGRID_H_

namespace OpenMS
{

class HashGrid {

private:
	DistanceSet distances;
	double rt_threshold;
	double mz_threshold;
	double hyp_threshold;
	double rt_scaling;
	int grid_size_x;
	int grid_size_y;
	double min_distance;
	std::pair<DataSubset*,DataSubset*> min_distance_subsets;



public:
	std::map<std::pair<int,int>, std::list<GridElement*> > elements;
	HashGrid()
		{
		}

	HashGrid(double rt_threshold_,double mz_threshold_)
	{
		rt_threshold=rt_threshold_;
		mz_threshold=mz_threshold_;
		grid_size_x=-1;
		grid_size_y=-1;
	}
	~HashGrid() {
		for (std::map<std::pair<int,int>, std::list<GridElement*> >::iterator it=elements.begin();it!=elements.end();++it)
		{
			std::list<GridElement*>& elements=it->second;
			for (std::list<GridElement*>::iterator lit=elements.begin();lit!=elements.end();++lit)
			{
				delete(*lit);
			}
		}

	}
	void setRTThreshold(double threshold_)
	{
		rt_threshold=threshold_;
	}

	void setMZThreshold(double threshold_)
	{
		mz_threshold=threshold_;
	}

	double getDistance(DataSubset& element1,DataSubset& element2)
	{
		DistanceSet::iterator pos=distances.find(boost::make_tuple(&element1,&element2));
		if (pos!=distances.end())
		{
			DistanceEntry el=*pos;
			return el.distance;
		}
		pos=distances.find(boost::make_tuple(&element2,&element1));
		if (pos!=distances.end())
		{
			DistanceEntry el=*pos;
			return el.distance;
		}
		return -1;
	}

	void removeElement(GridElement & element_,int x,int y)
		{
		std::list<GridElement*>& subsets = elements[std::make_pair(x,y)];
			subsets.remove(&element_);
			if (subsets.empty())
				elements.erase(std::make_pair(x,y));
		}

	void removeElement(GridElement & element_)
	{
		int x = element_.mz / mz_threshold;
		int y = element_.rt / rt_threshold;
		removeElement(element_,x,y);
	}



	void insert(GridElement* element_)
	{
			int x = element_->mz / mz_threshold;
			if (x>grid_size_x)
				grid_size_x=x;
			int y = element_->rt / rt_threshold;
			if (y>grid_size_y)
				grid_size_y=y;

			elements[std::make_pair(x,y)].push_back(element_);
	}

	void consoleOut()
	{
		for (std::map<std::pair<int,int>, std::list<GridElement*> >::iterator it=elements.begin();it!=elements.end();++it)
		{
			std::pair<int,int> coords=it->first;
			std::list<GridElement*> act_elements= it->second;
			if (it->second.size()>0)
				std::cout << coords.first << "/" << coords.second<< ": ";
			for (std::list<GridElement*>::iterator list_it=act_elements.begin();list_it!=act_elements.end();++list_it)
			{
				std::cout << (*list_it)->getID() << " | ";
			}
			std::cout << std::endl;

		}
		std::cout << std::endl;
	}

	int size()
	{
		return elements.size();
	}

    double getRT_threshold() const
    {
        return rt_threshold;
    }

    double getMZ_threshold() const
    {
        return mz_threshold;
    }

};
}




#endif /* HASHGRID_H_ */
