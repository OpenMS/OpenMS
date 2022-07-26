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