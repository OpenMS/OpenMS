// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/SYSTEM/File.h>
#include <boost/regex.hpp>

using namespace OpenMS;
using namespace std;

void CsiFingerIdMzTabWriter::read(const std::vector<String>& sirius_output_paths,
                                  const String& original_input_mzml,
                                  const Size& top_n_hits,
                                  MzTab& result)
{

  CsiFingerIdMzTabWriter::CsiAdapterRun csi_result;

  for (const auto& it : sirius_output_paths)
  {
    // extract mz, rt of the precursor and the nativeID of the corresponding MS2 spectra in the spectrum.ms file
    SiriusMzTabWriter::SiriusSpectrumMSInfo info = SiriusMzTabWriter::extractSpectrumMSInfo(it);

    const std::string pathtocsicsv = it + "/structure_candidates.tsv";

    ifstream file(pathtocsicsv);

    if (file) 
    {
      CsvFile compounds(pathtocsicsv, '\t');
      const UInt rowcount = compounds.rowCount();

      if (rowcount > 1)
      {
        // correction if the rowcount is smaller than the number of hits used as parameter
        // rowcount-1 because the csv header will be skipped in the loop later on.
        int header = 1;
        const UInt top_n_hits_cor = (top_n_hits >= rowcount) ? rowcount-header : top_n_hits;

        // fill identification structure containing all candidate hits for a single spectrum
        CsiFingerIdMzTabWriter::CsiAdapterIdentification csi_id{};

        // extract scan_index from path
        OpenMS::String str = File::path(pathtocsicsv);
        int scan_index = SiriusMzTabWriter::extractScanIndex(str);
    
        // extract scan_number from string
        int scan_number = SiriusMzTabWriter::extractScanNumber(str);

        // extract feature_id from string
        String feature_id = SiriusMzTabWriter::extractFeatureId(str);

        // extract column name and index from header
        std::map< std::string, Size > columnname_to_columnindex = SiriusMzTabWriter::extract_columnname_to_columnindex(compounds);

        // j = 1 because of .csv file format (header)
        for (Size j = 1; j <= top_n_hits_cor; ++j)
        {
          StringList sl;
          compounds.getRow(j, sl);
          CsiFingerIdMzTabWriter::CsiAdapterHit csi_hit;
          csi_hit.inchikey2D = sl[columnname_to_columnindex.at("InChIkey2D")];
          csi_hit.inchi = sl[columnname_to_columnindex.at("InChI")];
          csi_hit.molecular_formula = sl[columnname_to_columnindex.at("molecularFormula")];
          csi_hit.rank = sl[columnname_to_columnindex.at("rank")].toInt();
          csi_hit.formula_rank = sl[columnname_to_columnindex.at("formulaRank")].toInt();
          csi_hit.adduct = sl[columnname_to_columnindex.at("adduct")];
          csi_hit.score = sl[columnname_to_columnindex.at("CSI:FingerIDScore")].toDouble();
          csi_hit.name = sl[columnname_to_columnindex.at("name")];
          csi_hit.smiles = sl[columnname_to_columnindex.at("smiles")];
          csi_hit.xlogp = sl[columnname_to_columnindex.at("xlogp")];
          csi_hit.dbflags = sl[columnname_to_columnindex.at("dbflags")];
          sl[columnname_to_columnindex.at("pubchemids")].split(';', csi_hit.pubchemids);
          sl[columnname_to_columnindex.at("links")].split(';', csi_hit.links);

          csi_id.hits.push_back(csi_hit);
        }

        csi_id.mz = info.ext_mz;
        csi_id.rt = info.ext_rt;
        csi_id.native_ids = info.ext_n_id;
        csi_id.scan_index = scan_index;
        csi_id.scan_number = scan_number;
        csi_id.feature_id = feature_id;
        csi_result.identifications.push_back(csi_id);

        // write metadata to mzTab file
        MzTabFile mztab_out;
        MzTabMetaData md;
        MzTabMSRunMetaData md_run;
        md_run.location = MzTabString(original_input_mzml);
        md.ms_run[1] = md_run;
        md.description = MzTabString("CSI:FingerID-4.6.0");

        //needed for header generation (score)
        std::map<Size, MzTabParameter> smallmolecule_search_engine_score;
        smallmolecule_search_engine_score[1].setName("CSI:FingerIDScore");
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

            map <Size, MzTabDouble> engine_score;
            engine_score[1] = MzTabDouble(hit.score);
            smsr.best_search_engine_score = engine_score;

            smsr.chemical_formula = MzTabString(hit.molecular_formula);
            smsr.description = MzTabString(hit.name);
            std::vector <MzTabString> pubchemids;
            for (Size k = 0; k < hit.pubchemids.size(); ++k)
            {
              pubchemids.emplace_back(MzTabString(hit.pubchemids[k]));
            }  
            smsr.identifier.set(pubchemids);
            smsr.inchi_key = MzTabString(hit.inchikey2D);
            smsr.smiles = MzTabString(hit.smiles);
            std::vector < MzTabString > links;
            MzTabStringList m_links;
            m_links.setSeparator('|');
            for (Size k = 0; k < hit.links.size(); ++k)
            {
              links.emplace_back(MzTabString(hit.links[k]));
            }
            m_links.set(links);

            smsr.exp_mass_to_charge = MzTabDouble(id.mz);

            vector<MzTabDouble> v_rt;
            MzTabDoubleList rt_list;
            v_rt.emplace_back(id.rt);
            rt_list.set(v_rt);
            smsr.retention_time = rt_list;
            
            MzTabOptionalColumnEntry rank = make_pair("opt_global_rank", MzTabString(hit.rank));
            MzTabOptionalColumnEntry formula_rank = make_pair("opt_global_formulaRank", MzTabString(hit.formula_rank));
            MzTabOptionalColumnEntry compoundId = make_pair("opt_global_compoundId", MzTabString(id.scan_index));
            MzTabOptionalColumnEntry compoundScanNumber = make_pair("opt_global_compoundScanNumber", MzTabString(id.scan_number));
            MzTabOptionalColumnEntry featureId = make_pair("opt_global_featureId", MzTabString(id.feature_id));
            MzTabOptionalColumnEntry adduct = make_pair("opt_global_adduct", MzTabString(hit.adduct));
            MzTabOptionalColumnEntry xlogp = make_pair("opt_global_rank", MzTabString(hit.xlogp));
            MzTabOptionalColumnEntry dblinks = make_pair("opt_global_dblinks", MzTabString(m_links.toCellString()));
            MzTabOptionalColumnEntry dbflags = make_pair("opt_global_dbflags", MzTabString(hit.dbflags));

            vector<MzTabString> m_native_ids;
            MzTabStringList ml_native_ids;
            ml_native_ids.setSeparator('|');
            for (auto& element : id.native_ids)
            {
              m_native_ids.emplace_back(MzTabString(element));
            }
            ml_native_ids.set(m_native_ids);

            MzTabOptionalColumnEntry native_ids = make_pair("opt_global_native_id", MzTabString(ml_native_ids.toCellString()));

            smsr.opt_.push_back(rank);
            smsr.opt_.push_back(compoundId);
            smsr.opt_.push_back(compoundScanNumber);
            smsr.opt_.push_back(featureId);
            smsr.opt_.push_back(native_ids);
            smsr.opt_.push_back(adduct);
            smsr.opt_.push_back(xlogp);
            smsr.opt_.push_back(dblinks);
            smsr.opt_.push_back(dbflags);
            smsd.push_back(smsr);
          } 
        }
        result.setSmallMoleculeSectionRows(smsd);
      }
    }
    file.close();
  }
}

/// @endcond
