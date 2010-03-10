// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

//Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmApplyGivenTrafo.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

// debugging yes/no
// #define V_MapAlignmentAlgorithm(a) std::cout << a << std::endl;
#define V_MapAlignmentAlgorithm(a)

namespace OpenMS
{
	//register products here
	void MapAlignmentAlgorithm::registerChildren()
	{
		Factory<MapAlignmentAlgorithm>::registerProduct(
			MapAlignmentAlgorithmIdentification::getProductName(),
			&MapAlignmentAlgorithmIdentification::create);		

		Factory<MapAlignmentAlgorithm>::registerProduct ( MapAlignmentAlgorithmPoseClustering::    getProductName(), &MapAlignmentAlgorithmPoseClustering::    create );
		Factory<MapAlignmentAlgorithm>::registerProduct ( MapAlignmentAlgorithmSpectrumAlignment:: getProductName(), &MapAlignmentAlgorithmSpectrumAlignment:: create );
		Factory<MapAlignmentAlgorithm>::registerProduct ( MapAlignmentAlgorithmApplyGivenTrafo::   getProductName(), &MapAlignmentAlgorithmApplyGivenTrafo::   create );
	}

	MapAlignmentAlgorithm::MapAlignmentAlgorithm()
		: DefaultParamHandler("MapAlignmentAlgorithm"),
			ProgressLogger()
	{
	}

	MapAlignmentAlgorithm::~MapAlignmentAlgorithm()
	{
	}

	void MapAlignmentAlgorithm::alignPeakMaps(std::vector< MSExperiment<> >&, std::vector<TransformationDescription>&)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}

	void MapAlignmentAlgorithm::alignFeatureMaps(std::vector< FeatureMap<> >&, std::vector<TransformationDescription>&)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);				
	}

	void MapAlignmentAlgorithm::alignPeptideIdentifications(std::vector< std::vector< PeptideIdentification > >&, std::vector<TransformationDescription>&)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);				
	}



  void MapAlignmentAlgorithm::transformPeakMaps( std::vector< MSExperiment<> >& maps,
                                                                const  std::vector<TransformationDescription>& given_trafos
                                                              )
  {
    V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::transformPeakMaps()");

    if ( given_trafos.size() != maps.size() )
    {
      throw Exception::IllegalArgument
        (__FILE__, __LINE__, __PRETTY_FUNCTION__,
         String("MapAlignmentAlgorithm expects one given transformation (got: ")
         + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal"
        );
    }

    for ( UInt index = 0; index < maps.size(); ++index )
    {
      MSExperiment<> & mse = maps[index];
      const TransformationDescription& td = given_trafos[index];
      transformSinglePeakMap(mse,td);
    }
    return;
  }

  void MapAlignmentAlgorithm::transformSinglePeakMap( MSExperiment<>& msexp, const TransformationDescription& trafo )
  {
    msexp.clearRanges();
    trafo.init_();
    for ( MSExperiment<>::iterator mse_iter = msexp.begin(); mse_iter != msexp.end(); ++mse_iter )
    {
      DoubleReal rt = mse_iter->getRT();
      (*trafo.trafo_)(rt);
      mse_iter->setRT(rt);
    }
    msexp.updateRanges();
    return;
  }


  void MapAlignmentAlgorithm::transformFeatureMaps( std::vector< FeatureMap<> >& maps,
                                                                    const std::vector<TransformationDescription>& given_trafos
                                                                  )
   {
     V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::transformFeatureMaps()");

     if ( given_trafos.size() != maps.size() )
     {
       throw Exception::IllegalArgument
         (__FILE__, __LINE__, __PRETTY_FUNCTION__,
          String("MapAlignmentAlgorithm expects one given transformation (got: ")
          + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal"
         );
     }

     for ( UInt index = 0; index < maps.size(); ++index )
     {
       FeatureMap<> & fm = maps[index];
       const TransformationDescription& td = given_trafos[index];
       transformSingleFeatureMap( fm, td );
     }

     return;
   }

   void MapAlignmentAlgorithm::transformSingleFeatureMap( FeatureMap<>& fmap, const TransformationDescription& trafo )
   {
     trafo.init_();
     for ( std::vector<Feature>::iterator fmit = fmap.begin(); fmit != fmap.end(); ++fmit )
     {
       applyToFeature_(fmit, trafo);
     }
		 
		 // adapt RT values of unassigned peptides:
		 if (!fmap.getUnassignedPeptideIdentifications().empty())
		 {
			 transformSinglePeptideIdentification(
				 fmap.getUnassignedPeptideIdentifications(), trafo);
		 }

     return;
   }


   void MapAlignmentAlgorithm::applyToFeature_( const std::vector<Feature>::iterator &iter, const TransformationDescription& trafo )
   {
     // V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::applyToFeature_()");

     // transform feature position
     DoubleReal rt = iter->getRT();
     (*trafo.trafo_)(rt);
     iter->setRT(rt);

     // loop over all convex hulls
     std::vector<ConvexHull2D> & convex_hulls = iter->getConvexHulls();
     for ( std::vector<ConvexHull2D>::iterator chiter = convex_hulls.begin();
           chiter!= convex_hulls.end();
           ++chiter
         )
     {
       // transform all hull point positions within convex hull
       ConvexHull2D::PointArrayType & points = const_cast<ConvexHull2D::PointArrayType&>(chiter->getPoints());
       for ( ConvexHull2D::PointArrayType::iterator points_iter = points.begin();
             points_iter != points.end();
             ++points_iter
           )
       {
         DoubleReal rt = (*points_iter)[Feature::RT];
         (*trafo.trafo_)(rt);
         (*points_iter)[Feature::RT] = rt;
       }
     }

		 // adapt RT values of annotated peptides:
		 if (!iter->getPeptideIdentifications().empty())
		 {
			 transformSinglePeptideIdentification(iter->getPeptideIdentifications(),
																						trafo);
		 }
		 
     // recurse into subordinates
     for ( std::vector<Feature>::iterator subiter = iter->getSubordinates().begin();
           subiter != iter->getSubordinates().end();
           ++subiter )
     {
       applyToFeature_(subiter,trafo);
     }

     return;
   }



   void MapAlignmentAlgorithm::transformPeptideIdentifications( std::vector< std::vector< PeptideIdentification > >& maps,
                                                                                const std::vector<TransformationDescription>& given_trafos
                                                                              )
    {
      V_MapAlignmentAlgorithm("Hi out there.  This is MapAlignmentAlgorithm::transformPeptideIdentifications()");

      if ( given_trafos.size() != maps.size() )
      {
        throw Exception::IllegalArgument
          (__FILE__, __LINE__, __PRETTY_FUNCTION__,
           String("MapAlignmentAlgorithm expects one given transformation (got: ")
           + given_trafos.size() + ") per input map (got: " + maps.size() + "), these numbers are not equal"
          );
      }

      for ( UInt map_index = 0; map_index < maps.size(); ++map_index )
      {
        V_MapAlignmentAlgorithm("map_index: " << map_index);
        const TransformationDescription& td = given_trafos[map_index];
        std::vector< PeptideIdentification >& pepids = maps[map_index];
        transformSinglePeptideIdentification(pepids,td);
      }
      return;
    }

    void MapAlignmentAlgorithm::transformSinglePeptideIdentification( std::vector< PeptideIdentification >& pepids, const TransformationDescription& trafo )
    {
      const UInt meta_index_RT = MetaInfo::registry().getIndex("RT");
      trafo.init_();
      for ( UInt pepid_index = 0; pepid_index < pepids.size(); ++pepid_index )
      {
        V_MapAlignmentAlgorithm("pepid_index: " << pepid_index);
        PeptideIdentification & pepid = pepids[pepid_index];
        DataValue dv = pepid.getMetaValue(meta_index_RT);
        if (dv!=DataValue::EMPTY)
        {
          DoubleReal rt(dv);
          V_MapAlignmentAlgorithm("RT before: " << rt);
          (*trafo.trafo_)(rt);
          pepid.setMetaValue(meta_index_RT,rt);
          V_MapAlignmentAlgorithm("RT after: " << rt);
        }
      }
      return;
    }

} 
