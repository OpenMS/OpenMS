#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
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
		double x,y;
		x = min(second,first);
		y = max(first,second);
		if((y-x)>1)
		{
			throw Exception::InvalidRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}
		return 1-(y-x);
	}
}; // end of LowLevelComparator

Int main()
{
	// data
	vector<double> data; // must be filled

	srand(333);
	for(int i=0;i<12;++i)
	{
		data.push_back((double)rand()/RAND_MAX);
	}

	LowLevelComparator llc;
	SingleLinkage sl;
	vector< vector<UInt> > result; // will be filled
	ClusterHierarchical ch;
	ch.setThreshold(0.15);

	// clustering
	ch.clusterForVector<double,LowLevelComparator>(data,llc,sl,result);
	for(vector< vector<UInt> >::iterator outer_it=result.begin();outer_it!= result.end();++outer_it)
	{
		for(vector<UInt>::iterator inner_it=outer_it->begin();inner_it!= outer_it->end();++inner_it)
		{
			cout << " | " << *inner_it ;
		}
		cout << endl;
	}

	result.clear();

	ch.setThreshold(1.0);
	ch.clusterForDendrogramm<double,LowLevelComparator>(data,llc,sl,result,"output/Tutorial_Clustering.den");


	return 0;
} //end of main
