// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/MAPMATCHING/DMapMatcherRegression.h>
#include <OpenMS/KERNEL/DPeakConstReferenceArray.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>


#define DEBUG_ALIGNMENT
#undef DEBUG_ALIGNMENT

namespace OpenMS
{
  /**
     @brief A star alignment class.
     
     This class aligns elements of multiple element maps. 
     An element can be a DPeak, a DFeature, a ConsensusPeak or a ConsensusFeature.
     Corresponding elements are grouped together and stored as a consensus element.
     This class computes a star-alignment, that means a reference map is chosen 
     and in a first step all maps are mapped onto the reference map using the 
     PoseClusteringPairwiseMapMatcher. 
     In the second the final_consensus_map_ is determined. 
     At the beginning of the alignment the final_consensus_map_ contains all elements 
     of the reference_map as single consensus elements.
     Each transformed map is successive aligned to the final_consensus_map_ and corresponding
     elements are grouped together. 
     At the end of the alignment the resulting final_consensus_map_ covers the elements 
     of all maps, whereas corresponding elements are arranged together into ConsensusFeature or ConsensusPeak.
     
     @Note If you use consensus maps, the consensus elements are used as normal elements and you will
     loose the former consensus information.
     
     @ingroup MapAlignment
  */
  template < typename ConsensusElementT = ConsensusFeature< FeatureMap > >
  class StarAlignment : public BaseAlignment< ConsensusElementT >
  {
  public:
    /// Defines the coordinates of peaks / features.
    enum DimensionId
    {
      RT = DimensionDescription < LCMS_Tag >::RT,
      MZ = DimensionDescription < LCMS_Tag >::MZ
    };

    /** Symbolic names for indices of feature maps etc.
      This should make things more understandable and maintainable. 
      */
    enum Maps
    {
      MODEL = 0,
      SCENE = 1
    };

    /// Base class
    typedef BaseAlignment< ConsensusElementT > Base;
    typedef typename Base::ConsensusElementType ConsensusElementType;
    typedef typename Base::ElementType ElementType;
    typedef typename Base::ElementContainerType ElementContainerType;
    typedef typename Base::ConsensusMapType ConsensusMapType;
    typedef typename Base::GridType GridType;

    /// Pointer vector
    typedef DPeakConstReferenceArray< ElementContainerType > PeakConstReferenceMapType;

    /// Traits type
    typedef typename ConsensusElementType::TraitsType TraitsType;

    /// Position
    typedef DPosition < 2, TraitsType > PositionType;

    /// Quality
    typedef typename TraitsType::QualityType QualityType;

    //// Intensity
    typedef typename TraitsType::IntensityType IntensityType;

    /// Type of element pairs
    typedef DFeaturePair < 2, ConsensusElementType > ElementPairType;

    /// Container for generated consensus element pairs
    typedef DFeaturePairVector < 2, ConsensusElementType > ConsensusElementPairVectorType;

    /// Container for generated element pairs
    typedef DFeaturePairVector < 2, ElementType > ElementPairVectorType;

    using Base::element_map_vector_;
    using Base::param_;
    using Base::final_consensus_map_;
    using Base::transformations_;
    using Base::file_names_;
    using Base::map_type_;
    
    /// Constructor
    StarAlignment()
        : Base(),
        reference_map_index_(0)
    {}

    /// Copy constructor
    StarAlignment(const StarAlignment& source)
        : Base(source),
        reference_map_index_(source.reference_map_index_)
    {}

    ///  Assignment operator
    virtual StarAlignment& operator = (StarAlignment source)
    {
      if (&source==this)
        return *this;

      Base::operator = (source);
     
      reference_map_index_ = source.reference_map_index_;

      return *this;
    }

    /// Destructor
    virtual ~StarAlignment()
  {}

    /// Mutable access to the index of the reference map
    void setReferenceMapIndex(UnsignedInt index) throw (Exception::InvalidValue)
    {
      if (index > element_map_vector_.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The index is not contained in the vector of element containers.","") ;
      }
      else
      {
        reference_map_index_ = index;
      }
    }
    /// Non-mutable access to the index of the reference map
    const UnsignedInt& getReferenceMapIndex() const
    {
      return reference_map_index_;
    }

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

      if (map_type_ == "feature_map")
      {
        alignMultipleFeatureMaps_();
      }
      else
      {
        alignMultiplePeakMaps_();
      }
    }

    /// Return the Alignment tree in Newick Tree format
    virtual String getAlignmentTree() const
    {
      String tree;
      UnsignedInt n = element_map_vector_.size();
      tree = '(';
      UnsignedInt j = 0;
      for (UnsignedInt i = 0; i < n; ++i)
      {
        if (i != reference_map_index_)
        {
          tree = tree + '(' + reference_map_index_ + ":0," + i + ':' + (i+1) + "):0";
          ++j;

          if (j < (n-1))
          {
            tree = tree + ',';
          }
        }
      }
      tree = tree + ')';

      return tree;
    }

  protected:
    /// Index of the reference map
    UnsignedInt reference_map_index_;


    /// Define the map with the most elements as the reference map
    void assignReferenceMap_()
    {
      UnsignedInt n = element_map_vector_.size();
      UnsignedInt ref_index = 0;
      UnsignedInt max_number = element_map_vector_[ref_index]->size();

      for (UnsignedInt i = 1; i < n; ++i)
      {
        if (n > max_number)
        {
          max_number = n;
          ref_index = i;
        }
      }
      reference_map_index_ = ref_index;
    }

    /// Align all feature maps to the reference map
    void alignMultipleFeatureMaps_()
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
        out << cons_ref_map[i].getPosition()[RT] << ' ' << cons_ref_map[i].getPosition()[MZ] << '\n';
      }
      out.flush();

      std::cout << "*** Compute the consensus map of all pairwise alignment ***" << std::endl;
#endif
      // compute the consensus map of all pairwise alignment
      Param param_matcher = param_.copy("matching_algorithm:",true);
      /// Pairwise map matcher

      BasePairwiseMapMatcher< ConsensusMapType >* pairwise_matcher_;
      DataValue data_value = param_matcher.getValue("type");
      if (data_value != DataValue::EMPTY)
      {
        pairwise_matcher_ = Factory<BasePairwiseMapMatcher< ConsensusMapType > >::create(data_value);
        pairwise_matcher_->setParameters(param_matcher);
      }
      else
      {
        pairwise_matcher_ = Factory<BasePairwiseMapMatcher< ConsensusMapType > >::create(data_value);
        param_.setValue("matching_algorithm:type","poseclustering_pairwise");
        pairwise_matcher_->setParameters(param_matcher);
      }
      
      pairwise_matcher_->setElementMap(MODEL,cons_ref_map);

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
          std::cout << "*** Build a consensus map of map " << i << " *** " << file_names_[i] << " ***" <<  std::endl;
#endif

          //build a consensus map of map i
          ConsensusMapType map;
          buildConsensusMapType_(i,map);

#ifdef DEBUG_ALIGNMENT

          std::cout << "*** Compute a transformation for each grid cell and find pairs in the reference_map_ and map " << i << " ***" << std::endl;
#endif
          // compute a transformation for each grid cell and find pairs in the reference_map_ and map_i
          pairwise_matcher_->setElementMap(SCENE, map);
          pairwise_matcher_->clearGrid();
          pairwise_matcher_->run();

#ifdef DEBUG_ALIGNMENT

          std::cout << "*** Estimate for each grid cell a better transformation using the element pairs. number of pairs: " << pairwise_matcher_->getElementPairs().size() << " ***" << std::endl;
#endif

          // use the linear regression only if there are more than 2 pairs
          if (pairwise_matcher_->getElementPairs().size() > 2)
          {
            // estimate for each grid cell a better transformation using the element pairs
            lin_regression.setFeaturePairs(pairwise_matcher_->getElementPairs());
            lin_regression.setGrid(pairwise_matcher_->getGrid());
            lin_regression.setMinQuality(-1.);
            lin_regression.estimateTransform();
            transformations_[i] = lin_regression.getGrid();
          }
          // otherwise take the estimated transformation of the superimposer
          else
          {
            transformations_[i] = pairwise_matcher_->getGrid();
          }
#ifdef DEBUG_ALIGNMENT

          String name = "map_" + (String)number_alignments + ".dat";
          std::ofstream out(name.c_str(), std::ios::out);
#endif
          // iterate over all Elements...
          UnsignedInt n = map.size();
          for (UnsignedInt j = 0; j < n; ++j)
          {
            //             std::cout << "insert " << map[j] << std::endl;
            // Test in which cell this element is included
            // and apply the corresponding transformation
            typename GridType::iterator grid_it = transformations_[i].begin();
            while ((grid_it != (transformations_[i]).end()))
            {
              IndexTuple< ElementContainerType > index_tuple(i,j,(*(element_map_vector_[i]))[j]);
              PositionType pos = (*(element_map_vector_[i]))[j].getPosition();
              if (grid_it->encloses(map[j].getPosition()))
              {
                // apply transform for the singleton group element
                if (grid_it->getMappings().size() != 0)
                {
                  DLinearMapping<1>* mapping_rt = dynamic_cast<DLinearMapping<1>* >(grid_it->getMappings()[RT]);
                  DLinearMapping<1>* mapping_mz = dynamic_cast<DLinearMapping<1>* >(grid_it->getMappings()[MZ]);

                  mapping_rt->apply(pos[RT]);
                  mapping_mz->apply(pos[MZ]);
                }

                index_tuple.setTransformedPosition(pos);
#ifdef DEBUG_ALIGNMENT

                out << map[j].getPosition()[RT] << ' ' << map[j].getPosition()[MZ] << ' ' << pos[RT] << ' ' << pos[MZ] << '\n';
#endif

                map[j].getPosition() = pos;
                map[j].insert(index_tuple);
                break;
              }
              grid_it++;
            } // end while (grid)
          } // end for (Elements)

#ifdef DEBUG_ALIGNMENT
          out.flush();
          std::cout << "*** Compute the consensus of the reference map and map " << i << " ***" << std::endl;
#endif
          // compute the consensus of the reference map and map i
          DelaunayPairFinder<ConsensusMapType, ElementContainerType> pair_finder;
          pair_finder.setParameters(param_.copy("consensus_algorithm:",true));
          pair_finder.computeConsensusMap(map,final_consensus_map_);

#ifdef DEBUG_ALIGNMENT

          std::cout << "*** DONE!! number of consensus elements " << final_consensus_map_.size() << " ***"<< std::endl;
          ++number_alignments;
          std::ofstream out_cons("ConsensusMap",std::ios::out);
          for (UnsignedInt i = 0; i < final_consensus_map_.size(); ++i)
          {
            out_cons << final_consensus_map_[i] << std::endl;
          }
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
        out_cons << c->getPosition()[RT] << ' '
        << c->getPosition()[MZ] << ' '
        << c->getIntensity() << ' ';

        for (typename ConsensusElementType::Group::const_iterator it = c->begin(); it != c->end(); ++it)
        {
          out_cons << it->getElement().getPosition()[RT] << ' '
          << it->getTransformedPosition()[RT] << ' '
          << it->getElement().getPosition()[MZ] << ' '
          << it->getTransformedPosition()[MZ] << ' '
          << it->getElement().getIntensity() << ' ';
        }
        out_cons << std::endl;
      }

      std::ofstream out_gp("Consensus.gp",std::ios::out);
      UnsignedInt first=5;
      UnsignedInt second=7;
      out_gp << "plot \"reference_map.dat\" using 1:2 title \"reference_map\"  w points pointtype 20 lt 1\n"
      << "replot \"Consensus.dat\" using 1:2:($" << first << "-$1):($" << second << "-$2)  w vectors lt 3 nohead title \"pairs\"\n"
      << "replot \"Consensus.dat\" using 1:2 title \"consensus\"  w points pointtype 20 lt 2\n"
      << "replot \"Consensus.dat\" using " << first << ':' << second << " title \"\" w points pointtype 20 lt 1\n";
      UnsignedInt n=element_map_vector_.size();
      first +=5;
      second +=5;
      for (UnsignedInt i=0; i < (n-1); ++i)
      {
        String map = "map_" + (String)i + ".dat";
        out_gp << "replot \"Consensus.dat\" using 1:2:($" << first << "-$1):($" << second << "-$2)  w vectors lt 3 nohead title \"\"\n"
        << "replot \"" << map << "\" using 1:2 title \"original positions map " << i << "\" pointtype 3 lt " << i+3 << '\n'
        << "replot \"" << map << "\" using 3:4 title \"transformed positions map " << i << "\" pointtype 20 lt " << i+3 << '\n'
        << "replot \"" << map << "\" using 1:2:($3-$1):($4-$2) w vectors lt 7 nohead title \"transformed\"\n";
        first +=5;
        second +=5;
      }
      std::cout << "The consensus elements are written to Consensus.dat.\n"
      << "You can visualize the result using the gnuplot script \"Consensus.gp\" (Type \"gnuplot Consensus.gp -\")" << std::endl;
#endif

      delete pairwise_matcher_;

    }

    /// Align all peak maps to the reference map
    void alignMultiplePeakMaps_()
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
        out << cons_ref_map[i].getPosition()[RT] << ' ' << cons_ref_map[i].getPosition()[MZ] << '\n';
      }
      out.flush();

      std::cout << "*** Compute the consensus map of all pairwise alignment ***" << std::endl;
#endif
      // compute the consensus map of all pairwise alignment
      Param param_matcher = param_.copy("matching_algorithm:",true);
      
      // take the n-th most intensive Peaks of the reference map
      Size n = 50;
      PeakConstReferenceMapType reference_pointer_map((element_map_vector_[reference_map_index_])->begin(), (element_map_vector_[reference_map_index_])->end());
      reference_pointer_map.sortByIntensity();
      Size number = (reference_pointer_map.size() > n) ? n : reference_pointer_map.size();
      PeakConstReferenceMapType reference_most_intense(reference_pointer_map.end() - number, reference_pointer_map.end());

      BasePairwiseMapMatcher< PeakConstReferenceMapType >* pairwise_matcher_;
      DataValue data_value = param_matcher.getValue("type");
      if (data_value != DataValue::EMPTY)
      {
        pairwise_matcher_ = Factory<BasePairwiseMapMatcher< PeakConstReferenceMapType > >::create(data_value);
        pairwise_matcher_->setParameters(param_matcher);
      }
      else
      {
        pairwise_matcher_ = Factory<BasePairwiseMapMatcher< PeakConstReferenceMapType > >::create(data_value);
        param_.setValue("matching_algorithm:type","poseclustering_pairwise");
        pairwise_matcher_->setParameters(param_matcher);
      }
      pairwise_matcher_->setElementMap(MODEL,reference_most_intense);

      DMapMatcherRegression< ElementType > lin_regression;
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

          PeakConstReferenceMapType pointer_map((element_map_vector_[i])->begin(), (element_map_vector_[i])->end());
          pairwise_matcher_->clearGrid();
          pairwise_matcher_->initGridTransformation(pointer_map);
          /* pointer_map.sortByIntensity();
           Size number = (pointer_map.size() > n) ? n : pointer_map.size();
           PeakConstReferenceMapType most_intense(pointer_map.end() - number, pointer_map.end());
           */

#ifdef DEBUG_ALIGNMENT

          std::cout << "*** Compute a transformation for each grid cell and find pairs in the reference_map_ and map " << i << " ***" << std::endl;
#endif
          // compute a transformation for each grid cell and find pairs in the reference_map_ and map_i
          pairwise_matcher_->setElementMap(SCENE, pointer_map);
          pairwise_matcher_->run();

#ifdef DEBUG_ALIGNMENT

          std::cout << "*** Estimate for each grid cell a better transformation using the element pairs. number of pairs: " << pairwise_matcher_->getElementPairs().size() << " ***" << std::endl;
#endif
          // estimate for each grid cell a better transformation using the element pairs
          lin_regression.setFeaturePairs(pairwise_matcher_->getElementPairs());
          lin_regression.setGrid(pairwise_matcher_->getGrid());
          lin_regression.setMinQuality(-1.);
          lin_regression.estimateTransform();

          transformations_[i] = lin_regression.getGrid();

#ifdef DEBUG_ALIGNMENT

          String name = "map_" + (String)number_alignments + ".dat";
          std::ofstream out(name.c_str(), std::ios::out);
#endif
          // iterate over all Elements...
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

                out << map[j].getPosition()[RT] << ' ' << map[j].getPosition()[MZ] << ' ' << pos[RT] << ' ' << pos[MZ] << '\n';
#endif

                map[j].getPosition() = pos;
                map[j].insert(index_tuple);
              }
              grid_it++;

            } // end while (grid)
          } // end for

#ifdef DEBUG_ALIGNMENT
          out.flush();

          std::cout << "*** Compute the consensus of the reference map and map " << i << " ***" << std::endl;
#endif
          // compute the consensus of the reference map and map i
          DelaunayPairFinder<ConsensusMapType, ElementContainerType> pair_finder;
          pair_finder.setParameters(param_.copy("consensus_algorithm:",true));
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
        out_cons << c->getPosition()[RT] << ' '
        << c->getPosition()[MZ] << ' '
        << c->getIntensity() << ' ';

        for (typename ConsensusElementType::Group::const_iterator it = c->begin(); it != c->end(); ++it)
        {
          out_cons << it->getElement().getPosition()[RT] << ' '
          << it->getTransformedPosition()[RT] << ' '
          << it->getElement().getPosition()[MZ] << ' '
          << it->getTransformedPosition()[MZ] << ' '
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
      n = element_map_vector_.size();
      for (UnsignedInt i=0; i < n; ++i)
      {
        if (i != reference_map_index_)
        {
          String map = "map_" + (String)i + ".dat";
          out_gp << "replot \"Consensus.dat\" using 1:2:($" << first << "-$1):($" << second << "-$2)  w vectors lt 3 nohead title \"pairs\"\n"
          << "replot \"" << map << "\" using 1:2 title \"original positions map " << i << "\" pointtype 3 lt " << i+3 << '\n'
          << "replot \"" << map << "\" using 3:4 title \"transformed positions map " << i << "\" pointtype 20 lt " << i+3 << '\n'
          << "replot \"" << map << "\" using 1:2:($3-$1):($4-$2) w vectors lt 7 nohead title \"transformed\"\n";
          first +=5;
          second +=5;
        }
      }
      std::cout << "The consensus elements are written to Consensus.dat.\n"
      << "You can visualize the result using the gnuplot script \"Consensus.gp\" (Type \"gnuplot Consensus.gp -\")" << std::endl;
#endif

      delete pairwise_matcher_;
    }
  }
  ; // StarAlignment
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_STARALIGNMENT_H
