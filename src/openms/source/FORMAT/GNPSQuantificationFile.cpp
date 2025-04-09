// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------
//-------------------------------------------------------------------------

#include <OpenMS/FORMAT/GNPSQuantificationFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <fstream>
#include <iostream>
#include <unordered_map>

namespace OpenMS
{
    /**
    @brief Generates a feature quantification file required for GNPS FBMN, as defined here: https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/#feature-quantification-table
    */
    void GNPSQuantificationFile::store(const ConsensusMap& consensus_map, const String& output_file)
    {
        // IIMN meta values will be exported, if first feature contains mv Constants::UserParam::IIMN_ROW_ID
        bool iimn = false;
        if (consensus_map[0].metaValueExists(Constants::UserParam::IIMN_ROW_ID)) iimn = true;

        // meta values for ion identity molecular networking
        std::vector<String> iimn_mvs{Constants::UserParam::IIMN_ROW_ID,
                                    Constants::UserParam::IIMN_BEST_ION,
                                    Constants::UserParam::IIMN_ADDUCT_PARTNERS,
                                    Constants::UserParam::IIMN_ANNOTATION_NETWORK_NUMBER};
        
        // initialize SVOutStream with tab separation
        std::ofstream outstr(output_file.c_str());
        SVOutStream out(outstr, "\t", "_", String::NONE);
        
        // write headers for MAP and CONSENSUS
        out << "#MAP" << "id" << "filename" << "label" << "size" << std::endl;
        out << "#CONSENSUS" << "rt_cf" << "mz_cf" << "intensity_cf" << "charge_cf" << "width_cf" << "quality_cf";
        if (iimn)
        {
        for (const auto& mv : iimn_mvs) out << mv;
        }
        for (size_t i = 0; i < consensus_map.getColumnHeaders().size(); i++)
        {
        out << "rt_" + String(i) << "mz_" + String(i) << "intensity_" + String(i) << "charge_" + String(i) << "width_" + String(i);
        }
        out << std::endl;

        // write MAP information
        for (const auto& h: consensus_map.getColumnHeaders())
        {
        out << "MAP" << h.first << h.second.filename << h.second.label << h.second.size << std::endl;
        }

        // write ConsensusFeature information
        for (const auto& cf: consensus_map)
        {
        out << "CONSENSUS" << cf.getRT() << cf.getMZ() << cf.getIntensity() << cf.getCharge() << cf.getWidth() << cf.getQuality();
        if (iimn)
        {
            for (const auto& mv : iimn_mvs) out << cf.getMetaValue(mv, "");
        }
        // map index to feature handle and write feature information on correct position, if feature is missing write empty strings
        std::unordered_map<size_t, FeatureHandle> index_to_feature;
        for (const auto& fh: cf.getFeatures()) index_to_feature[fh.getMapIndex()] = fh;
        for (size_t i = 0; i < consensus_map.getColumnHeaders().size(); i++)
        {
            if (index_to_feature.count(i))
            {
            out << index_to_feature[i].getRT() << index_to_feature[i].getMZ() << index_to_feature[i].getIntensity() << index_to_feature[i].getCharge() << index_to_feature[i].getWidth();
            }
            else
            {
            out << "" << "" << "" << "" << "";
            }
        }
        out << std::endl;
        }
        outstr.close();
    }
} // namespace