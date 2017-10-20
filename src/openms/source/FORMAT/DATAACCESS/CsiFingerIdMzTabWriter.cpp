// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/CsvFile.h>

#include <OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>

using namespace OpenMS;
using namespace std;

void CsiFingerIdMzTabWriter::read(const std::vector<String> & paths, Size number, MzTab & result)
{

  CsiFingerIdMzTabWriter::CsiAdapterRun csi_result;

  for (std::vector<String>::const_iterator it = paths.begin(); it != paths.end(); ++it)
  {

    const std::string pathtocsicsv = *it + "/summary_csi_fingerid.csv";

    ifstream file(pathtocsicsv);

    if (file) 
    {
      CsvFile compounds(pathtocsicsv, '\t');
      const UInt rowcount = compounds.rowCount();
    
      if (rowcount > 1)
      {
        
        // fill identification structure containing all candidate hits for a single spectrum
        CsiFingerIdMzTabWriter::CsiAdapterIdentification csi_id;

        //Extract scan_index from path
        OpenMS::String str = File::path(pathtocsicsv);
        std::string scan_index = SiriusMzTabWriter::extract_scan_index(str);

        const UInt number_cor = (number > rowcount) ? rowcount : number;
        for (Size j = 1; j < number_cor; ++j)
        {
          
          StringList sl;
          compounds.getRow(j, sl);
          CsiFingerIdMzTabWriter::CsiAdapterHit csi_hit;
          csi_hit.inchikey2D = sl[0];
          csi_hit.inchi = sl[1];
          csi_hit.molecular_formula = sl[2];
          csi_hit.rank = sl[3].toInt();
          csi_hit.score = sl[4].toDouble();
          csi_hit.name = sl[5];
          csi_hit.smiles = sl[6];
          sl[8].split(';', csi_hit.pubchemids);
          sl[9].split(';', csi_hit.links);

          csi_id.hits.push_back(csi_hit);
        }

        csi_id.scan_index = scan_index;
        csi_result.identifications.push_back(csi_id);

        // write metadata to mzTab file
        MzTabFile mztab_out;
        MzTabMetaData md;
        MzTabMSRunMetaData md_run;
        md_run.location = MzTabString(str);
        md.ms_run[1] = md_run;
        md.description = MzTabString("CSI:FingerID-3.5");

        //needed for header generation (score)
        std::map<Size, MzTabParameter> smallmolecule_search_engine_score;
        smallmolecule_search_engine_score[1].setName("score");
        md.smallmolecule_search_engine_score = smallmolecule_search_engine_score;
        result.setMetaData(md);

        // write results to mzTab file
        MzTabSmallMoleculeSectionRows smsd;
        for (Size i = 0; i < csi_result.identifications.size(); ++i)
        {
          const CsiFingerIdMzTabWriter::CsiAdapterIdentification &id = csi_result.identifications[i];
          for (Size j = 0; j < id.hits.size(); ++j)
          {
            const CsiFingerIdMzTabWriter::CsiAdapterHit &hit = id.hits[j];
            MzTabSmallMoleculeSectionRow smsr;

            map <Size, MzTabDouble> engine_score = {{1, MzTabDouble(hit.score)}};
            smsr.best_search_engine_score = engine_score;

            smsr.chemical_formula = MzTabString(hit.molecular_formula);
            smsr.description = MzTabString(hit.name);
            std::vector <MzTabString> pubchemids;
            for (Size k = 0; k < hit.pubchemids.size(); ++k)
            {
              pubchemids.push_back(MzTabString(hit.pubchemids[k]));
            }  
            smsr.identifier.set(pubchemids);
            smsr.inchi_key = MzTabString(hit.inchikey2D);
            smsr.smiles = MzTabString(hit.smiles);
            std::vector < MzTabString > uri;
            for (Size k = 0; k < hit.links.size(); ++k)
            {
              uri.push_back(MzTabString(hit.links[k]));
            }  

            MzTabOptionalColumnEntry rank;
            rank.first = "rank";
            rank.second = MzTabString(hit.rank);

            MzTabOptionalColumnEntry compoundId;
            compoundId.first = "compoundId";
            compoundId.second = MzTabString(id.scan_index);

            smsr.opt_.push_back(rank);
            smsr.opt_.push_back(compoundId);
            smsd.push_back(smsr);
          } 
        }
        result.setSmallMoleculeSectionRows(smsd);
      }
    }
  }
}

/// @endcond
