// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_STARALIGNMENT_H
#define OPENMS_ANALYSIS_MAPMATCHING_STARALIGNMENT_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/GeomHashPairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DMapDewarper.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DMapMatcherRegression.h>

#define DEBUG_ALIGNMENT
//#undef DEBUG_ALIGNMENT

namespace OpenMS
{

  /**
     @brief 
   

  **/
  template < typename ConsensusElementT = ConsensusFeature<> >
  class StarAlignment : public BaseAlignment< ConsensusElementT >
  {
    public:
      /// Defines the coordinates of peaks / features.
      enum DimensionId
      {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

      /** @name Type definitions
       */
      //@{
      /// Base class
      typedef BaseAlignment< ConsensusElementT > Base;
      typedef typename Base::ConsensusElementType ConsensusElementType;
      typedef typename Base::ElementType ElementType;
      typedef typename Base::ElementContainerType ElementContainerType;
      typedef typename Base::ConsensusMapType ConsensusMapType;

      ///
      typedef DPeakConstReferenceArray< ConsensusMapType > PeakConstReferenceMapType;

      /// Traits type
      typedef typename ConsensusElementType::TraitsType TraitsType;

      /// Position
      typedef DPosition < 2, TraitsType > PositionType;

      /// Quality
      typedef typename TraitsType::QualityType QualityType;

      //// Intensity
      typedef typename TraitsType::IntensityType IntensityType;

      /// Type of feature pairs
      typedef DFeaturePair < 2, ConsensusElementType > ElementPairType;

      /// Container for generated feature pairs
      typedef DFeaturePairVector < 2, ConsensusElementType > ElementPairVectorType;

      /// Type of estimated transformation
      typedef DGrid< 2 > GridType;
      //@}

      using Base::element_map_vector_;
      using Base::reference_map_index_;
      using Base::param_;
      using Base::final_consensus_map_;
      using Base::file_names_;
      using Base::assignReferenceMap_;
      using Base::pairwise_matcher_;
      using Base::pair_finder_;


      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      StarAlignment()
          : Base()
      {}

      /// Copy constructor
      StarAlignment(const StarAlignment& source)
          : Base(source)
      {}

      ///  Assignment operator
      virtual StarAlignment& operator = (StarAlignment source)
      {
        if (&source==this)
          return *this;

        return *this;
      }

      /// Destructor
      virtual ~StarAlignment()
    {}
      //@}

      /** @name Accesssor methods
       */
      //@{
      /// Mutable access to the transformations
      void setTransformationVector(const std::vector< GridType >& transformations)
      {
        transformations_ = transformations;
      }
      /// Mutable access to the transformations
      std::vector< GridType >& getTransformationVector()
      {
        return transformations_;
      }
      /// Non-mutable access to the transformations
      const std::vector< GridType >& getTransformationVector() const
      {
        return transformations_;
      }
      //@}


      /// Estimates the transformation for each grid cell
      virtual void run() throw (Exception::InvalidValue)
      {
        if (element_map_vector_.size() < 2)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Take at least 2 maps for alignment.","") ;
        }
        else
        {
          assignReferenceMap_();
        }
#ifdef DEBUG_ALIGNMENT
        std::cout << "*** Reference Map is " << file_names_[reference_map_index_] << " ***" <<  std::endl;
#endif

        alignMultipleMaps_();
      }

      /// Return the Alignment tree in Newick Tree format
      virtual String getAlignmentTree() const
      {
        String tree;
        UnsignedInt n = element_map_vector_.size();
        tree = '(';
        for (UnsignedInt i = 0; i < n; ++i)
        {
          if (i != reference_map_index_)
          {
            tree = tree + '(' + (String)reference_map_index_ + ":0," + (String)i + ':' + (String)(i+1) + "):0";

            if (i < (n-1))
            {
              tree = tree + ',';
            }
          }
        }
        tree = tree + ')';

        return tree;
      }

    protected:

      /** @name Data members
       */
      //@{
      /// transformations
      std::vector< GridType > transformations_;
      //@}


      // align all maps to the reference map
      void alignMultipleMaps_()
      {

#ifdef DEBUG_ALIGNMENT
        std::cout << "*** Build a consensus map of the elements of the reference map (contains only singleton consensus elements) ***" << std::endl;
#endif
        // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
        ConsensusMapType cons_ref_map;
        buildConsensusMapTypeInsertInGroup_(reference_map_index_,cons_ref_map);
        final_consensus_map_ = cons_ref_map;

#ifdef DEBUG_ALIGNMENT
        std::ofstream out("reference_map.dat", std::ios::out);
        for (UnsignedInt i = 0; i < cons_ref_map.size(); ++i)
        {
          out << cons_ref_map[i].getPosition()[0] << ' ' << cons_ref_map[i].getPosition()[1] << '\n';
        }
        out.flush();

        std::cout << "*** Compute the consensus map of all pairwise alignment ***" << std::endl;


#endif
        // compute the consensus map of all pairwise alignment
        Param param_matcher = param_.copy("matching:",true);
        pairwise_matcher_->setParam(param_matcher);
        pairwise_matcher_->setFeatureMap(0,cons_ref_map);

        DMapMatcherRegression<ConsensusElementType> lin_regression;
        UnsignedInt number_maps = element_map_vector_.size();
        transformations_.resize(number_maps);
#ifdef DEBUG_ALIGNMENT

        UnsignedInt number_alignments = 0;
#endif

        for (UnsignedInt i = 0; i < number_maps; ++i)
        {
          std::cout.precision(10);
          if (i != reference_map_index_)
          {
#ifdef DEBUG_ALIGNMENT
            std::cout << "*** Build a consensus map of map " << i << " *** " << std::endl;
#endif
            //build a consensus map of map i
            ConsensusMapType map;
            buildConsensusMapType_(i,map);

#ifdef DEBUG_ALIGNMENT

            std::cout << "*** Compute a transformation for each grid cell and find pairs in the reference_map_ and map " << i << " ***" << std::endl;
#endif
            // compute a transformation for each grid cell and find pairs in the reference_map_ and map_i
            pairwise_matcher_->setFeatureMap(1, map);
            ElementPairVectorType pairs;
            pairwise_matcher_->setFeaturePairs(pairs);
            pairwise_matcher_->run();

#ifdef DEBUG_ALIGNMENT

            std::cout << "*** Estimate for each grid cell a better transformation using the element pairs. number of pairs: " << pairs.size() << " ***" << std::endl;
#endif
            // estimate for each grid cell a better transformation using the element pairs
            lin_regression.setFeaturePairs(pairs);
            lin_regression.setGrid(pairwise_matcher_->getGrid());
            lin_regression.setMinQuality(-1.);
            lin_regression.estimateTransform();

            transformations_[i] = lin_regression.getGrid();
#ifdef DEBUG_ALIGNMENT

            String name = "map_" + (String)number_alignments + ".dat";
            std::ofstream out(name.c_str(), std::ios::out);
#endif
            // iterate over all features...
            UnsignedInt n = map.size();
            for (UnsignedInt j = 0; j < n; ++j)
            {
              // Test in which cell this element is included
              // and apply the corresponding transformation
              typename GridType::iterator grid_it = (lin_regression.getGrid()).begin();
              while (grid_it != (lin_regression.getGrid()).end() )
              {
                if (grid_it->encloses(map[j].getPosition()) )
                {
                  DLinearMapping<1>* mapping_rt = dynamic_cast<DLinearMapping<1>* >(grid_it->getMappings()[RT]);
                  DLinearMapping<1>* mapping_mz = dynamic_cast<DLinearMapping<1>* >(grid_it->getMappings()[MZ]);

                  // apply transform for the singleton group element
                  IndexTuple< ElementContainerType > index_tuple(i,j,(*(element_map_vector_[i]))[j]);
                  PositionType pos = (*(element_map_vector_[i]))[j].getPosition();

                  mapping_rt->apply(pos[RT]);
                  mapping_mz->apply(pos[MZ]);
                  index_tuple.setTransformedPosition(pos);
#ifdef DEBUG_ALIGNMENT

                  out << map[j].getPosition()[0] << ' ' << map[j].getPosition()[1] << ' ' << pos[0] << ' ' << pos[1] << '\n';
#endif

                  map[j].getPosition() = pos;
                  map[j].insert(index_tuple);
                }
                grid_it++;

              } // end while (grid)
            } // end while (features)
#ifdef DEBUG_ALIGNMENT
            out.flush();

            std::cout << "*** Compute the consensus of the reference map and map " << i << " ***" << std::endl;
#endif
            // compute the consensus of the reference map and map i
            DelaunayPairFinder<ConsensusMapType> pair_finder;
            pair_finder.setParam(param_matcher);
            pair_finder.computeConsensusMap(map,final_consensus_map_);

#ifdef DEBUG_ALIGNMENT
            std::cout << "*** DONE!! number of consensus elements " << final_consensus_map_.size() << " ***"<< std::endl;
            ++number_alignments;
#endif

          }
        }
#ifdef DEBUG_ALIGNMENT
        std::cout << "=========== Final Consensus Map =========" << std::endl;
        std::ofstream out_cons("Consensus.dat",std::ios::out);
        out_cons << "cons_rt cons_mz cons_int rt_map1 rt_transf_map1 mz_map1 mz_transf_map1 int_map1 rt_map2 rt_transf_map1 mz_map2 mz_transf_map2 int_map2 ... rt_mapn rt_transf_mapn mz_mapn mz_transf_mapn int_mapn\n";
        for (UnsignedInt i = 0; i < final_consensus_map_.size(); ++i)
        {
          ConsensusElementType* c = &(final_consensus_map_[i]);
          out_cons << c->getPosition()[0] << ' '
          << c->getPosition()[1] << ' '
          << c->getIntensity() << ' ';

          for (typename ConsensusElementType::Group::const_iterator it = c->begin(); it != c->end(); ++it)
          {
            out_cons << it->getElement().getPosition()[0] << ' '
            << it->getTransformedPosition()[0] << ' '
            << it->getElement().getPosition()[1] << ' '
            << it->getTransformedPosition()[1] << ' '
            << it->getElement().getIntensity() << ' ';
          }
          out_cons << std::endl;
        }

        std::ofstream out_gp("Consensus.gp",std::ios::out);
        UnsignedInt first=5;
        UnsignedInt second=7;
        out_gp << "plot \"reference_map.dat\" using 1:2 title \"reference_map\"  w points pointtype 20 lt 1\n"
               << "replot \"Consensus.dat\" using 1:2 title \"consensus\"  w points pointtype 20 lt 2\n"
               << "replot \"Consensus.dat\" using " << first << ':' << second << " title \"\" w points pointtype 20 lt 1\n";
        UnsignedInt n=element_map_vector_.size();
        for (UnsignedInt i=0; i < (n-1); ++i)
        {
          String map = "map_" + (String)i + ".dat";
          out_gp << "replot \"Consensus.dat\" using 1:2:($" << first << "-$1):($" << second << "-$2)  w vectors lt 3 nohead title \"pairs\"\n"
          << "replot \"" << map << "\" using 1:2 title \"original positions map " << i << "\" pointtype 3 lt " << i+3 << '\n'
          << "replot \"" << map << "\" using 3:4 title \"transformed positions map " << i << "\" pointtype 20 lt " << i+3 << '\n'
          << "replot \"" << map << "\" using 1:2:($3-$1):($4-$2) w vectors lt 7 nohead title \"transformed\"\n";
          first +=5;
          second +=5;
        }
        std::cout << "The consensus elements are written to Consensus.dat.\n"
        << "You can visualize the result using the gnuplot script \"Consensus.gp\" (Type \"gnuplot Consensus.gp -\")" << std::endl;
#endif

      }
  }
  ; // StarAlignment
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_STARALIGNMENT_H
