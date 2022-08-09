// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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