// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------
//-------------------------------------------------------------------------

#include <OpenMS/FORMAT/GNPSMetaValueFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <fstream>
#include <iostream>
#include <unordered_map>

namespace OpenMS
{
    /**
    @brief Generates a meta value table required for GNPS FBMN, as defined here: https://ccms-ucsd.github.io/GNPSDocumentation/metadata/
    */
    void GNPSMetaValueFile::store(const ConsensusMap& consensus_map, const String& output_file)
    {   
        StringList mzML_file_paths;
        consensus_map.getPrimaryMSRunPath(mzML_file_paths);
        std::ofstream outstr(output_file.c_str());
        SVOutStream out(outstr, "\t", "_", String::NONE);

        out << "" << "filename" << "ATTRIBUTE_MAPID" << std::endl;
        Size i = 0;
        for (const auto& path: mzML_file_paths)
        {
            out << String(i) << path.substr(path.find_last_of("/\\")+1) << "MAP"+String(i) << std::endl;
            i++;
        }
    }
} // namespace