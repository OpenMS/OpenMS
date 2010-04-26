/*
 * SILACData.h
 *
 *  Created on: 25.03.2010
 *      Author: steffen
 */

#ifndef GRIDELEMENT_H_
#define GRIDELEMENT_H_

#include<map>
#include<list>
#include<set>
#include<vector>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>


namespace OpenMS
{
class OPENMS_DLLAPI GridElement {
public :
	double mz;
	double rt;
	virtual int getID();
};

class OPENMS_DLLAPI DataPoint : public GridElement {

public:
	double intensity; // intensity at RT and m/z
	int cluster_id; // ID number of the cluster the data point belongs to
	int cluster_size; // number of points in cluster 'cluster_id'
	int feature_id;
	DataPoint()
	{
		rt = 0;
		mz = 0;
		intensity = 0;
		feature_id= 0;
		cluster_id = 0;
		cluster_size = 0;
	}
	DataPoint(const DataPoint &copyin) : GridElement(copyin)
	{
		rt = copyin.rt;
		mz = copyin.mz;
		feature_id=copyin.feature_id;
		intensity = copyin.intensity;
		cluster_id = copyin.cluster_id;
		cluster_size = copyin.cluster_size;
	}
	DataPoint(double rt_, double mz_, double intensity_, int feature_id_)
	{
		rt = rt_;
		mz = mz_;
		feature_id=feature_id_;
		intensity = intensity_;
		cluster_id = 0;
		cluster_size = 0;
	}
	~DataPoint(){};
	DataPoint& operator=(const DataPoint &rhs)
	{
		this->rt = rhs.rt;
		this->mz = rhs.mz;
		this->intensity = rhs.intensity;
		this->feature_id=rhs.feature_id;
		this->cluster_id = rhs.cluster_id;
		this->cluster_size = rhs.cluster_size;
		return *this;
	}

	bool operator==(const DataPoint &rhs) const
	{
		if( this->feature_id != rhs.feature_id) return false;
		if( this->rt != rhs.rt) return false;
		if( this->mz != rhs.mz) return false;
		if( this->intensity != rhs.intensity) return false;
		//if( this->cluster_id != rhs.cluster_id) return false;
		//if( this->cluster_size != rhs.cluster_size) return false;
		return true;
	}


	bool operator!=(const DataPoint &rhs) const
	{
		return !(*this==rhs);
	}
	bool operator<(const DataPoint &rhs) const // required for built-in STL functions like sort
	{
		if ( *this == rhs) return false;

		return this->feature_id < rhs.feature_id;
	}

	int getID()
	{
		return feature_id;
	}
};


class OPENMS_DLLAPI BinaryTreeNode {
public:
	DataPoint* data1;
		DataPoint* data2;
		double distance;
	BinaryTreeNode(DataPoint* data1_,DataPoint* data2_,double distance_) {
		data1=data1_;
		data2=data2_;
		distance=distance_;
	}

	~BinaryTreeNode() {
		// TODO Auto-generated destructor stub
	}

	BinaryTreeNode& operator = (const BinaryTreeNode& source)
		{
			if (this != &source)
			{
				data1 = source.data1;
				data2 = source.data2;
				distance = source.distance;
			}
			return *this;
		}

	bool operator==(const BinaryTreeNode &cp) const
	{
		if( this->data1 != cp.data1) return false;
		if( this->data2 != cp.data2) return false;
		if( this->distance != cp.distance) return false;
		return true;
	}

	bool operator<(const BinaryTreeNode &cp) const
	{
		if (std::abs(this->distance - cp.distance) <= 0.00000001)
				return *(this->data1) < *(cp.data1);
			else
				return this->distance < cp.distance;
	}

};

struct Dist{};

class DataSubset;

struct DistanceEntry
{
	DataSubset* data_point;
	DataSubset* owner;
	double distance;
	DistanceEntry(DataSubset* owner_,DataSubset* data_point_,double distance_)
	{
		owner=owner_;
		data_point=data_point_;
		distance=distance_;
	}
	int operator<(const DistanceEntry &i) const
	{
		return distance < i.distance;
	}

};


typedef boost::multi_index::multi_index_container<
DistanceEntry,
boost::multi_index::indexed_by<
boost::multi_index::hashed_unique<
boost::multi_index::composite_key<
			DistanceEntry,
			boost::multi_index::member<DistanceEntry,DataSubset*,&DistanceEntry::owner>,
			boost::multi_index::member<DistanceEntry,DataSubset*,&DistanceEntry::data_point> > >,
			boost::multi_index::ordered_non_unique< boost::multi_index::tag<Dist>,
			boost::multi_index::member<DistanceEntry,double,&DistanceEntry::distance> > >
> DistanceSet;


class OPENMS_DLLAPI DataSubset : public GridElement
{

public:
	std::map<DataSubset*,DistanceSet::iterator> distance_iterators;
	std::list<DataPoint*> data_points;
	std::vector<BinaryTreeNode> tree;

	DataSubset(DataPoint& data_point)
	{
		data_points.push_back(&data_point);
		rt = data_point.rt;
		mz = data_point.mz;
	}
	DataSubset(const DataSubset& copy) :GridElement(copy)
	{
		data_points = copy.data_points;
		tree = copy.tree;
		distance_iterators = copy.distance_iterators;
		mz = copy.mz;
		rt = copy.rt;
	}
	DataSubset(const DataSubset* copy_ptr)
	{
		data_points=copy_ptr->data_points;
		tree=copy_ptr->tree;
		distance_iterators=copy_ptr->distance_iterators;
		mz=copy_ptr->mz;
		rt=copy_ptr->rt;
	}
	int operator<(const DataSubset &el) const
	{
		std::list<DataPoint*> data1=this->data_points;
		std::list<DataPoint*> data2=el.data_points;
		return data1.size() < data2.size();
	}


	int size()
	{
		return data_points.size();
	}

	int getID()
	{
		return data_points.front()->feature_id;
	}

	bool operator !=(const DataSubset &el) const
	{
		if (this->data_points.size()!=el.data_points.size())
		{
			return true;
		}
		else
		{
			std::list<DataPoint*> data1=this->data_points;
			std::list<DataPoint*>::iterator it1=data1.begin();
			std::list<DataPoint*> data2=el.data_points;
			std::list<DataPoint*>::iterator it2=data2.begin();
			while((*it1)->mz==(*it2)->mz && (*it1)->rt==(*it2)->rt)
			{
				++it1;
				++it2;
				if (it1==data1.end() && it2==data2.end())
					return false;
			}
			return true;
		}
	}

	bool operator ==(const DataSubset & el) const
	{
		return !(*this!=el);
	}
	void removeDistanceIterator(std::map<DataSubset*,double>::iterator pos)
	{
	}
};
}



#endif /* GRIDELEMENT_H_ */
