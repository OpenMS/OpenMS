// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julia Thueringer$
// $Authors: Julia Thueringer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

namespace OpenMS
{
  /**
    @brief A map alignment algorithm based on peptide identifications from MS2 spectra.

    ID groups with the same sequence in different maps represent points of correspondence in RT between the maps. They are used to evaluate the distances between the maps for hierarchical clustering and form the basis for the alignment.
    Only the best PSM per spectrum is considered as the correct identification.

    For each pair of maps, the similarity is determined based on the intersection of the contained identifications using Pearson correlation. For small intersections, the Pearson value is reduced by multiplying the ratio of the intersection size to the union size: \f$\texttt{PearsonValue(map1}\cap \texttt{map2)}*\Bigl(\frac{\texttt{N(map1 }\cap\texttt{ map2})}{\texttt{N(map1 }\cup\texttt{ map2})}\Bigr)\f$
    Using hierarchical clustering together with average linkage a binary tree is produced.
    Following the tree, the maps are aligned, resulting in a transformed feature map that contains both the original and the transformed retention times.
    As long as there are at least two clusters, the alignment is done as follows:
    Of every pair of clusters, the one with the larger 10/90 percentile retention time range is selected as reference for the align() method of @ref OpenMS::MapAlignmentAlgorithmIdentification.
    align() aligns the median retention time of each ID group in the second cluster to the reference retention time of this group.
    Cubic spline smoothing is used to convert this mapping to a smooth function.
    Retention times in the second cluster are transformed to the reference scale by applying this function.
    Additionally, the original retention times are stored in the meta information of each feature.
    The reference is combined with the transformed cluster.

    The resulting map is used to extract transformation descriptions for each input map.
    For each map cubic spline smoothing is used to convert the mapping to a smooth function.
    Retention times of each map are transformed by applying the smoothed function.

    @htmlinclude OpenMS_MapAlignmentAlgorithmTreeGuided.parameters

    @ingroup MapAlignment

  */
  class OPENMS_DLLAPI MapAlignmentAlgorithmTreeGuided :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    MapAlignmentAlgorithmTreeGuided();

    /// Destructor
    ~MapAlignmentAlgorithmTreeGuided() override;

    /**
     * @brief Extract RTs given for individual features of each map, calculate distances for each pair of maps and cluster hierarchical using average linkage.
     *
     * @param feature_maps Vector of input maps (FeatureMap) whose distance is to be calculated.
     * @param tree Vector of BinaryTreeNodes that will be computed
     * @param maps_ranges Vector to store all sorted RTs of extracted identifications for each map in @p feature_maps; needed to determine the 10/90 percentiles
    */
    static void buildTree(std::vector<FeatureMap>& feature_maps, std::vector<BinaryTreeNode>& tree, std::vector<std::vector<double>>& maps_ranges);

    /**
     * @brief Align feature maps tree guided using align() of @ref OpenMS::MapAlignmentAlgorithmIdentification and use TreeNode with larger 10/90 percentile range as reference.
     *
     * @param tree Vector of BinaryTreeNodes that contains order for alignment.
     * @param feature_maps_transformed Vector with input maps for transformation process. Because the transformed maps are stored within this vector it's not const.
     * @param maps_ranges Vector that contains all sorted RTs of extracted identifications for each map; needed to determine the 10/90 percentiles.
     * @param map_transformed FeatureMap to store all features of combined maps with original and transformed RTs in order of alignment.
     * @param trafo_order Vector to store indices of maps in order of alignment.
    */
    void treeGuidedAlignment(const std::vector<BinaryTreeNode>& tree, std::vector<FeatureMap>& feature_maps_transformed,
                             std::vector<std::vector<double>>& maps_ranges, FeatureMap& map_transformed,
                             std::vector<Size>& trafo_order);

    /**
     * @brief Align feature maps tree guided using align() of @ref OpenMS::MapAlignmentAlgorithmIdentification and use TreeNode with larger 10/90 percentile range as reference.
    */
    void align(std::vector<FeatureMap>& data,
               std::vector<TransformationDescription>& transformations);

    /**
     * @brief Extract original RT ("original_RT" MetaInfo) and transformed RT for each feature to compute RT transformations.
     *
     * @param feature_maps Vector of input maps for size information.
     * @param map_transformed FeatureMap that contains all features of combined maps with original and transformed RTs in order of alignment.
     * @param transformations Vector to store transformation descriptions for each map. (output)
     * @param trafo_order Vector that contains the indices of aligned maps in order of alignment.
    */
    void computeTrafosByOriginalRT(std::vector<FeatureMap>& feature_maps, FeatureMap& map_transformed,
                                   std::vector<TransformationDescription>& transformations, const std::vector<Size>& trafo_order);

    /**
     * @brief Apply transformations on input maps.
     *
     * @param feature_maps Vector of maps to be transformed (output)
     * @param transformations Vector that contains TransformationDescriptions that are applied to input maps
    */
    static void computeTransformedFeatureMaps(std::vector<FeatureMap>& feature_maps, const std::vector<TransformationDescription>& transformations);

protected:
    /// Type to store feature retention times given for individual peptide sequence
    typedef std::map<String, DoubleList> SeqAndRTList;

    // Update defaults model_type_, model_param_ and align_algorithm_
    void updateMembers_() override;

    /// Type of transformation model
    String model_type_;

    /// Default params of transformation models linear, b_spline, lowess and interpolated
    Param model_param_;

    /// Instantiation of alignment algorithm
    MapAlignmentAlgorithmIdentification align_algorithm_;

    /**
     * @brief Similarity functor that provides similarity calculations with the ()-operator for protected type SeqAndRTList.
     * SeqAndRTList stores retention times given for individual peptide sequences of a feature map.

      Using pearson correlation, calculate the retention time similarity of two maps from their intersection of the peptide identifications.
      Small intersections are penalized by multiplication with the quotient of intersection to union.
    */
    class PeptideIdentificationsPearsonDistance_;

    /**
     * @brief For given peptide identifications extract sequences and store with associated feature RT.
     *
     * @param peptides Vector of peptide identifications to extract sequences.
     * @param peptide_rts Map to store a list of feature RTs for each peptide sequence as key.
     * @param map_range Vector in which all feature RTs are stored for given peptide identifications.
     * @param feature_rt RT value of the feature to which the peptide identifications to be analysed belong.
     */
    static void addPeptideSequences_(const std::vector<PeptideIdentification>& peptides, SeqAndRTList& peptide_rts,
            std::vector<double>& map_range, double feature_rt);

    /**
     * @brief For each input map, extract peptide identifications (sequences) of existing features with associated feature RT.
     *
     * @param feature_maps Vector of original maps containing peptide identifications.
     * @param maps_seq_and_rt Vector of maps to store feature RTs given for individual peptide sequences for each feature map.
     * @param maps_ranges Vector to store all feature RTs of extracted identifications for each map; needed to determine the 10/90 percentiles.
     */
    static void extractSeqAndRt_(const std::vector<FeatureMap>& feature_maps, std::vector<SeqAndRTList>& maps_seq_and_rt,
            std::vector<std::vector<double>>& maps_ranges);

private:
    /// Copy constructor intentionally not implemented -> private
    MapAlignmentAlgorithmTreeGuided(const MapAlignmentAlgorithmTreeGuided&);

    /// Assignment operator intentionally not implemented -> private
    MapAlignmentAlgorithmTreeGuided& operator=(const MapAlignmentAlgorithmTreeGuided&);
  };
} // namespace OpenMS
