// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer_impl.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder_impl.h>


using namespace std;

//TODO merge the two algorithms as soon as possible

namespace OpenMS
{

	MapAlignmentAlgorithmPoseClustering::MapAlignmentAlgorithmPoseClustering()
		: MapAlignmentAlgorithm()
	{
		setName("MapAlignmentAlgorithmPoseClustering");
	
    defaults_.setValue("superimposer_type", "poseclustering_affine","The superimposer used ");
		defaults_.setValidStrings("superimposer_type",Factory<BaseSuperimposer< PointMapType > >::registeredProducts());
		defaults_.setValue("pairfinder_type", "DelaunayPairFinder","The pair finder used ");
		defaults_.setValidStrings("pairfinder_type",Factory<BasePairFinder< PointMapType > >::registeredProducts());
		
		subsections_.push_back("superimposer");
		subsections_.push_back("pairfinder");
		
		defaultsToParam_();
	}

	MapAlignmentAlgorithmPoseClustering::~MapAlignmentAlgorithmPoseClustering()
	{
	}

	void MapAlignmentAlgorithmPoseClustering::alignPeakMaps(vector< MSExperiment<> >& maps)
	{		
		//define reference map (the one with most peaks)
		UInt reference_map_index = 0;
		UInt max_count = 0;		
		for (UInt m=0; m<maps.size(); ++m)
		{
			//init getSize() by calling updateRanges
			maps[m].updateRanges(1);
			if (maps[m].getSize()>max_count)
			{
				max_count = maps[m].getSize();
				reference_map_index = m;
			}
		}
		
    // build a consensus map of the elements of the reference map (take the 400 highest peaks)
    PointMapType cons_ref_map;
    { //new scope to get rid of the tmp variables
			DPeakArray<RawDataPoint2D> tmp;
			tmp.sortByIntensity(true);
			if (tmp.size()>400) tmp.resize(400); //TODO make this a parameter
			maps[reference_map_index].get2DData(tmp);
	    for (UInt i=0; i < tmp.size(); ++i)
	    {
	      ConsensusFeature< DPeakArray<RawDataPoint2D> > c(reference_map_index,i,tmp[i]);
	      cons_ref_map.push_back(c);
	    }
	  }
		
		//init superimposer and pairfinder with model and parameters
		BaseSuperimposer<PointMapType>* superimposer = Factory<BaseSuperimposer<PointMapType> >::create(param_.getValue("superimposer_type"));
    superimposer->setElementMap(0, cons_ref_map);
    superimposer->setParameters(param_.copy("superimposer:",true));
    
    BasePairFinder<PointMapType>* pairfinder = Factory<BasePairFinder<PointMapType> >::create(param_.getValue("pairfinder_type"));
		pairfinder->setElementMap(0, cons_ref_map);
    pairfinder->setParameters(param_.copy("pairfinder:",true));
		
		for (UInt i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				//build a consensus map of map i
				PointMapType map;
		    DPeakArray<RawDataPoint2D> tmp;
		    maps[i].get2DData(tmp);
		    for (UInt i2=0; i2 < tmp.size(); ++i2)
		    {
		      ConsensusFeature< DPeakArray<RawDataPoint2D> >  c(tmp[i2].getPosition(),tmp[i2].getIntensity());
		      map.push_back(c);
		    }
				
				//run superimposer    
	      superimposer->setElementMap(1, map);
	      superimposer->run();
	      //run pairfinder
	      pairfinder->setTransformation(0, superimposer->getTransformation(0));        
	      pairfinder->setElementMap(1, map);
	      vector< ElementPair < ConsensusFeature< DPeakArray <RawDataPoint2D> > > > element_pairs;
	      pairfinder->setElementPairs(element_pairs);
	      pairfinder->findElementPairs();

				// calculate the transformation
				LinearMapping trafo;
				if (element_pairs.size() > 2) // estimate for each grid cell a better transformation using the element pairs
				{
					trafo = calculateRegression_(element_pairs);
				}
				else // otherwise take the estimated transformation of the superimposer
				{
					trafo = superimposer->getTransformation(0);
				}
				
				// apply transformation
				for (UInt j=0; j< maps[i].size(); ++j)
				{
					DoubleReal rt = maps[i][j].getRT();
					trafo.apply(rt);
					maps[i][j].setRT(rt);
				}
			}
		}
	}


	void MapAlignmentAlgorithmPoseClustering::alignFeatureMaps(vector< FeatureMap<> >& maps)
	{
		//define reference map (the one with most peaks)
		UInt reference_map_index = 0;
		UInt max_count = 0;		
		for (UInt m=0; m<maps.size(); ++m)
		{
			if (maps[m].size()>max_count)
			{
				max_count = maps[m].size();
				reference_map_index = m;
			}
		}
		
    // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
    FeatureMapType cons_ref_map;
    const FeatureMap<>& ref_map = maps[reference_map_index];
    for (UInt i=0; i < ref_map.size(); ++i)
    {
      cons_ref_map.push_back(ConsensusFeature< FeatureMap<> > (reference_map_index,i,ref_map[i]));
    }
   
		//init superimposer and pairfinder with model and parameters
		BaseSuperimposer<FeatureMapType>* superimposer = Factory<BaseSuperimposer<FeatureMapType> >::create(param_.getValue("superimposer_type"));
    superimposer->setElementMap(0, cons_ref_map);
    superimposer->setParameters(param_.copy("superimposer:",true));
    
    BasePairFinder<FeatureMapType>* pairfinder = Factory<BasePairFinder<FeatureMapType> >::create(param_.getValue("pairfinder_type"));
		pairfinder->setElementMap(0, cons_ref_map);
    pairfinder->setParameters(param_.copy("pairfinder:",true));

    for (UInt i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				//build a consensus map of map i
				FeatureMapType map;
		    const FeatureMap<>& map2 = maps[i];
		    map.reserve(map2.size());
		    for (UInt i2=0; i2 < map2.size(); ++i2)
		    {
		      map.push_back(ConsensusFeature< FeatureMap<> >(map2[i2].getPosition(),map2[i2].getIntensity()));
		    }

				//run superimposer    
	      superimposer->setElementMap(1, map);
	      superimposer->run();
	      //run pairfinder
	      pairfinder->setTransformation(0, superimposer->getTransformation(0));        
	      pairfinder->setElementMap(1, map);
	      vector< ElementPair < ConsensusFeature< FeatureMap<> > > > element_pairs;
	      pairfinder->setElementPairs(element_pairs);
	      pairfinder->findElementPairs();

				// calculate the transformation
				LinearMapping trafo;
				if (element_pairs.size() > 2) // estimate for each grid cell a better transformation using the element pairs
				{
					trafo = calculateRegression_(element_pairs);
				}
				else // otherwise take the estimated transformation of the superimposer
				{
					trafo =  superimposer->getTransformation(0);
				}

				// apply transformation
				for (UInt j = 0; j < map.size(); ++j)
				{
					DoubleReal rt = maps[i][j].getRT();
					trafo.apply(rt);
					maps[i][j].setRT(rt);
				}
			}
		}
	}

} //namespace 
