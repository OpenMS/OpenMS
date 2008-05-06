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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_STARALIGNMENT_H
#define OPENMS_ANALYSIS_MAPMATCHING_STARALIGNMENT_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapMatcherRegression.h>
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
     
		 @note If you use consensus maps, the consensus elements are used as normal elements and you will
		 loose the former consensus information.
     
		 @ingroup Analysis
	*/
	template < typename ConsensusElementT = ConsensusFeature< FeatureMap < > > >
	class StarAlignment : public BaseAlignment< ConsensusElementT >
	{
	public:
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
    typedef typename Base::ConsensusVectorType ConsensusVectorType;

    /// Pointer vector
    typedef DPeakConstReferenceArray< ElementContainerType > PeakConstReferenceMapType;
    
    /// Position
    typedef DPosition < 2 > PositionType;

    /// Quality
    typedef DoubleReal QualityType;
    
    //// Intensity
    typedef DoubleReal IntensityType;

    /// Type of element pairs
    typedef ElementPair < ConsensusElementType > ElementPairType;

    /// Container for generated consensus element pairs
    typedef std::vector < ConsensusElementType > ConsensusElementPairVectorType;

    /// Container for generated element pairs
    typedef std::vector < ElementType > ElementPairVectorType;

    using Base::param_;
    using Base::final_consensus_map_;
    using Base::transformations_;
    using Base::map_type_;

    /// Constructor
    StarAlignment()
			: Base(),
				reference_map_index_(0)
    {
    	//set the name for DefaultParamHandler error messages
    	this->setName("StarAlignment");	
      Base::subsections_.push_back("matching_algorithm");
      Base::subsections_.push_back("consensus_algorithm");
      
      Base::defaultsToParam_();
    }

    /// Destructor
    virtual ~StarAlignment()
    {}

    /// Mutable access to the index of the reference map
    void setReferenceMapIndex(UInt index) throw (Exception::InvalidValue)
    {
      if (index > final_consensus_map_.getMapVector().size())
				{
					throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The index is not contained in the vector of element containers.","") ;
				}
      else
				{
					reference_map_index_ = index;
				}
    }
    /// Non-mutable access to the index of the reference map
    UInt getReferenceMapIndex() const
    {
      return reference_map_index_;
    }

    /// Estimates the transformation for each grid cell
    virtual void run() throw (Exception::InvalidValue)
    {
      if (final_consensus_map_.getMapVector().size() < 2)
				{
					throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Take at least 2 maps for alignment.","") ;
				}
      else
				{
					assignReferenceMap_();
				}
#ifdef DEBUG_ALIGNMENT
      std::cout << "*** Reference Map is " << final_consensus_map_.getFilenames()[reference_map_index_] << " ***" <<  std::endl;
#endif

      if (map_type_ == "feature_map")
				{
					alignMultipleFeatureMaps_();
				}
      else
				{
					if (map_type_ == "consensus_map")
						{
							alignMultipleConsensusMaps_();
						} 
            else 
            {
              if (map_type_ == "peak_map")
						  {
							 alignMultiplePeakMaps_();
						  }
              else
              {
                throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Wrong map type. Choose \"feature_map\", \"peak_map\" or \"consensus_map\".","") ;
              }
            }
            
				}
    }
        
    
    /// Return the Alignment tree in Newick Tree format
    virtual String getAlignmentTree() const
    {
      String tree;
      UInt n = final_consensus_map_.getMapVector().size();
      tree = '(';
      UInt j = 0;
      for (UInt i = 0; i < n; ++i)
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
    
    /// Merge the elements of the final consensus map
    void merge(ConsensusMap < ConsensusElementType >& new_map)
    {
      final_consensus_map_.merge(new_map);
    }

  protected:
		/// Int of the reference map
    UInt reference_map_index_;


    /// Define the map with the most elements as the reference map
    void assignReferenceMap_()
    {
      UInt n = final_consensus_map_.getMapVector().size();
      UInt ref_index = 0;
      UInt max_number = final_consensus_map_.getMapVector()[ref_index]->size();

      for (UInt i = 1; i < n; ++i)
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
      std::vector < ElementContainerType* >& element_map_vector = final_consensus_map_.getMapVector();
      
      // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
      ConsensusVectorType cons_ref_map;
      const ElementContainerType& map = *(element_map_vector[reference_map_index_]);
      UInt n = map.size();
      for (UInt i=0; i < n; ++i)
      {
        ConsensusElementType c(reference_map_index_,i,map[i]);
        final_consensus_map_.push_back(c);
        cons_ref_map.push_back(c);
      }
      
      // compute the consensus map of all pairwise alignment
      /// Pairwise map matcher
			BasePairwiseMapMatcher< ConsensusVectorType >* pairwise_matcher_;
			pairwise_matcher_ = Factory<BasePairwiseMapMatcher< ConsensusVectorType > >::create("poseclustering_pairwise");
			pairwise_matcher_->setParameters(param_.copy("matching_algorithm:",true));

      pairwise_matcher_->setElementMap(MODEL,cons_ref_map);

      MapMatcherRegression<ConsensusElementType> lin_regression;
      UInt number_maps = element_map_vector.size();
      transformations_.resize(number_maps);

      for (UInt i = 0; i < number_maps; ++i)
				{
					std::cout.precision(10);
					if (i != reference_map_index_)
						{

							//build a consensus map of map i
							ConsensusVectorType map;
							buildConsensusVectorType_(i,map);
							
							// compute a transformation for each grid cell and find pairs in the reference_map_ and map_i
							pairwise_matcher_->setElementMap(SCENE, map);
							pairwise_matcher_->initGridTransformation(map);
							pairwise_matcher_->run();

							// use the linear regression only if there are more than 2 pairs
							if (pairwise_matcher_->getElementPairs().size() > 2)  
								{
									// estimate for each grid cell a better transformation using the element pairs
									lin_regression.setElementPairs(pairwise_matcher_->getElementPairs());
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

							// iterate over all Elements...
							UInt n = map.size();
							for (UInt j = 0; j < n; ++j)
							{
								IndexTuple index_tuple(i,j,(*(element_map_vector[i]))[j].getIntensity(),(*(element_map_vector[i]))[j].getPosition());
								PositionType pos = (*(element_map_vector[i]))[j].getPosition();
								// apply transform for the singleton group element
								transformations_[i].apply(pos[RawDataPoint2D::RT]);

								index_tuple.setTransformedPosition(pos);

								map[j].getPosition() = pos;
								map[j].insert(index_tuple);
							} // end for (Elements)


							// compute the consensus of the reference map and map i
							DelaunayPairFinder<ConsensusVectorType, ElementContainerType> pair_finder;
							pair_finder.setParameters(param_.copy("consensus_algorithm:",true));
							pair_finder.computeConsensusMap(map,final_consensus_map_);
						}
				}

      delete pairwise_matcher_;

    }

    /// Align all peak maps to the reference map
    void alignMultiplePeakMaps_()
		{
      std::vector < ElementContainerType* >& element_map_vector = final_consensus_map_.getMapVector();
      
      // compute the consensus map of all pairwise alignment
      // take the n-th most intensive Peaks of the reference map
			UInt n = 400;
			PeakConstReferenceMapType reference_pointer_map((element_map_vector[reference_map_index_])->begin(), (element_map_vector[reference_map_index_])->end());
			reference_pointer_map.sortByIntensity();
			UInt number = (reference_pointer_map.size() > n) ? n : reference_pointer_map.size();
			PeakConstReferenceMapType reference_most_intense(reference_pointer_map.end() - number, reference_pointer_map.end());

			BasePairwiseMapMatcher< PeakConstReferenceMapType >* pairwise_matcher_;
			pairwise_matcher_ = Factory<BasePairwiseMapMatcher< PeakConstReferenceMapType > >::create("poseclustering_pairwise");
			pairwise_matcher_->setParameters(param_.copy("matching_algorithm:",true));
			
			pairwise_matcher_->setElementMap(MODEL,reference_most_intense);

			MapMatcherRegression< ElementType > lin_regression;
			UInt number_maps = element_map_vector.size();
			transformations_.resize(number_maps);

			for (UInt i = 0; i < number_maps; ++i)
				{
					std::cout.precision(10);
					if (i != reference_map_index_)
						{
							PeakConstReferenceMapType pointer_map((element_map_vector[i])->begin(), (element_map_vector[i])->end());
							// compute a transformation for each grid cell and find pairs in the reference_map_ and map_i
							pairwise_matcher_->setElementMap(SCENE, pointer_map);
							pairwise_matcher_->initGridTransformation(pointer_map);
							pairwise_matcher_->run();

#ifdef DEBUG_ALIGNMENT
							std::cout << "*** Estimate for each grid cell a better transformation using the element pairs. number of pairs: " << pairwise_matcher_->getElementPairs().size() << " ***" << std::endl;
#endif
							// estimate for each grid cell a better transformation using the element pairs
							lin_regression.setElementPairs(pairwise_matcher_->getElementPairs());
							lin_regression.setGrid(pairwise_matcher_->getGrid());
							lin_regression.setMinQuality(-1.);
							lin_regression.estimateTransform();
							
							transformations_[i] = lin_regression.getGrid();
						}
				}

			delete pairwise_matcher_;
		}
      
		/// Align all peak maps to the reference map
    void alignMultipleConsensusMaps_()
		{
      std::vector < ElementContainerType* >& element_map_vector = final_consensus_map_.getMapVector();
      
#ifdef DEBUG_ALIGNMENT
      std::cout << "*** Build a consensus map of the elements of the reference map (contains only singleton consensus elements) ***" << std::endl;
#endif
      // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
      ConsensusVectorType cons_ref_map;
      const ElementContainerType& map = *(element_map_vector[reference_map_index_]);
      UInt m = map.size();
      for (UInt i=0; i < m; ++i)
      {
        ConsensusElementType c(reference_map_index_,i,map[i]);
        final_consensus_map_.push_back(c);
        cons_ref_map.push_back(c);
      }
   
#ifdef DEBUG_ALIGNMENT

      std::ofstream out("reference_map.dat", std::ios::out);
      for (UInt i = 0; i < cons_ref_map.size(); ++i)
				{
					out << cons_ref_map[i].getRT() << ' ' << cons_ref_map[i].getMZ() << '\n';
				}
      out.flush();

      std::cout << "*** Compute the consensus map of all pairwise alignment ***" << std::endl;
#endif
      // compute the consensus map of all pairwise alignment
      // take the n-th most intensive Peaks of the reference map
      UInt n = 50;
      PeakConstReferenceMapType reference_pointer_map((element_map_vector[reference_map_index_])->begin(), (element_map_vector[reference_map_index_])->end());
      reference_pointer_map.sortByIntensity();
      UInt number = (reference_pointer_map.size() > n) ? n : reference_pointer_map.size();
      PeakConstReferenceMapType reference_most_intense(reference_pointer_map.end() - number, reference_pointer_map.end());

			BasePairwiseMapMatcher< PeakConstReferenceMapType >* pairwise_matcher_;
			pairwise_matcher_ = Factory<BasePairwiseMapMatcher< PeakConstReferenceMapType > >::create("poseclustering_pairwise");
			pairwise_matcher_->setParameters(param_.copy("matching_algorithm:",true));
				
      pairwise_matcher_->setElementMap(MODEL,reference_most_intense);

      MapMatcherRegression< ElementType > lin_regression;
      UInt number_maps = element_map_vector.size();
      transformations_.resize(number_maps);
#ifdef DEBUG_ALIGNMENT

      UInt number_alignments = 0;
#endif

      for (UInt i = 0; i < number_maps; ++i)
				{
					std::cout.precision(10);
					if (i != reference_map_index_)
						{
#ifdef DEBUG_ALIGNMENT
							std::cout << "*** Build a consensus map of map " << i << " *** " << std::endl;
#endif
							//build a consensus map of map i
							ConsensusVectorType map;
							buildConsensusVectorType_(i,map);

							PeakConstReferenceMapType pointer_map((element_map_vector[i])->begin(), (element_map_vector[i])->end());
							// compute a transformation for each grid cell and find pairs in the reference_map_ and map_i
							pairwise_matcher_->setElementMap(SCENE, pointer_map);
							pairwise_matcher_->initGridTransformation(pointer_map);
							pairwise_matcher_->run();

#ifdef DEBUG_ALIGNMENT

							std::cout << "*** Estimate for each grid cell a better transformation using the element pairs. number of pairs: " << pairwise_matcher_->getElementPairs().size() << " ***" << std::endl;
#endif
							// estimate for each grid cell a better transformation using the element pairs
							lin_regression.setElementPairs(pairwise_matcher_->getElementPairs());
							lin_regression.setGrid(pairwise_matcher_->getGrid());
							lin_regression.setMinQuality(-1.);
							lin_regression.estimateTransform();

							transformations_[i] = lin_regression.getGrid();

#ifdef DEBUG_ALIGNMENT

							String name = "map_" + (String)number_alignments + ".dat";
							std::ofstream out(name.c_str(), std::ios::out);
#endif
							// iterate over all Elements...
							UInt n = map.size();
							for (UInt j = 0; j < n; ++j)
								{
									// apply transform for the singleton group element
									IndexTuple index_tuple(i,j,(*(element_map_vector[i]))[j].getIntensity(),(*(element_map_vector[i]))[j].getPosition());
									PositionType pos = (*(element_map_vector[i]))[j].getPosition();

									lin_regression.getGrid().apply(pos[RawDataPoint2D::RT]);
									index_tuple.setTransformedPosition(pos);

									map[j].getPosition() = pos;
									map[j].insert(index_tuple);
								} // end for

#ifdef DEBUG_ALIGNMENT
							out.flush();

							std::cout << "*** Compute the consensus of the reference map and map " << i << " ***" << std::endl;
#endif
							// compute the consensus of the reference map and map i
							DelaunayPairFinder<ConsensusVectorType, ElementContainerType> pair_finder;
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
      for (UInt i = 0; i < final_consensus_map_.size(); ++i)
				{
					ConsensusElementType* c = &(final_consensus_map_[i]);
					out_cons << c->getRT() << ' '
									 << c->getMZ() << ' '
									 << c->getIntensity() << ' ';

					for (typename ConsensusElementType::Group::const_iterator it = c->begin(); it != c->end(); ++it)
						{
							out_cons << it->getElement().getRT() << ' '
											 << it->getTransformedPosition()[RawDataPoint2D::RT] << ' '
											 << it->getElement().getMZ() << ' '
											 << it->getTransformedPosition()[RawDataPoint2D::MZ] << ' '
											 << it->getElement().getIntensity() << ' ';
						}
					out_cons << std::endl;
				}

      std::ofstream out_gp("Consensus.gp",std::ios::out);
      UInt first=5;
      UInt second=7;
      out_gp << "plot \"reference_map.dat\" using 1:2 title \"reference_map\"  w points pointtype 20 lt 1\n"
						 << "replot \"Consensus.dat\" using 1:2 title \"consensus\"  w points pointtype 20 lt 2\n"
						 << "replot \"Consensus.dat\" using " << first << ':' << second << " title \"\" w points pointtype 20 lt 1\n";
      n = element_map_vector.size();
      for (UInt i=0; i < n; ++i)
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
