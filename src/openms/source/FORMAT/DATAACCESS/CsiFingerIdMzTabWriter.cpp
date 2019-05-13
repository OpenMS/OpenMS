// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/METADATA/SpectrumLookup.h>

#include <boost/regex.hpp>

#include <OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>

using namespace OpenMS;
using namespace std;

void CsiFingerIdMzTabWriter::read(const std::vector<String> & sirius_output_paths,
                                  const String & original_input_mzml,
                                  const Size & top_n_hits,
                                  MzTab & result)
{

  CsiFingerIdMzTabWriter::CsiAdapterRun csi_result;

  for (const auto& it : sirius_output_paths)
  {
   
    // extract mz, rt and nativeID of the corresponding precursor spectrum in the spectrum.ms file
    String ext_nid;
    double ext_mz = 0.0;
    double ext_rt = 0.0;
    const String sirius_spectrum_ms = it + "/spectrum.ms";
    ifstream spectrum_ms_file(sirius_spectrum_ms);
    if (spectrum_ms_file)
    {
      const String nid_prefix = "##nid";
      const String rt_prefix = "#rt";
      const String pmass_prefix = ">parentmass";
      const String ms1peaks = ">ms1peaks";
      String line;
      while (getline(spectrum_ms_file, line))
      {
        if (line.hasPrefix(pmass_prefix))
        {
           ext_mz = String(line.erase(line.find(pmass_prefix), pmass_prefix.size())).toDouble();
        }
        else if (line.hasPrefix(rt_prefix))
        {
           ext_rt = String(line.erase(line.find(rt_prefix), rt_prefix.size())).toDouble();
        }
        else if (line.hasPrefix(nid_prefix))
        {
           ext_nid = line.erase(line.find(nid_prefix), nid_prefix.size());
        }
        else if (line.hasPrefix(ms1peaks))
        {
           break; // only run till >ms1peaks
        }
      }
      spectrum_ms_file.close();
    } 
   
    const std::string pathtocsicsv = it + "/summary_csi_fingerid.csv";

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
        CsiFingerIdMzTabWriter::CsiAdapterIdentification csi_id;

        // extract scan_index from path
        OpenMS::String str = File::path(pathtocsicsv);
        int scan_index = SiriusMzTabWriter::extract_scan_index(str);
    
        // extract scan_number from string
        boost::regex regexp("-(?<SCAN>\\d+)-");
        int scan_number = SpectrumLookup::extractScanNumber(str, regexp, false);

        // extract feature_id from string
        boost::smatch match;
        String feature_id;
        boost::regex regexp_feature("_(?<SCAN>\\d+)-");
        bool found = boost::regex_search(str, match, regexp_feature);
        if (found && match["SCAN"].matched) {feature_id = "id_" + match["SCAN"].str();}
        String unassigned = "null";

        // j = 1 because of .csv file format (header)
        for (Size j = 1; j <= top_n_hits_cor; ++j)
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

        csi_id.mz = ext_mz;
        csi_id.rt = ext_rt;
        csi_id.native_id = ext_nid;
        csi_id.scan_index = scan_index;
        csi_id.scan_number = scan_number;
        // check if results were assigned to a feature
        if (feature_id != "id_0")
        {
          csi_id.feature_id = feature_id;
        }
        else
        {
          csi_id.feature_id = unassigned;
        }
        csi_result.identifications.push_back(csi_id);

        // write metadata to mzTab file
        MzTabFile mztab_out;
        MzTabMetaData md;
        MzTabMSRunMetaData md_run;
        md_run.location = MzTabString(original_input_mzml);
        md.ms_run[1] = md_run;
        md.description = MzTabString("CSI:FingerID-4.0.1");

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

            map <Size, MzTabDouble> engine_score;
            engine_score[1] = MzTabDouble(hit.score);
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

            smsr.exp_mass_to_charge = MzTabDouble(id.mz);

            vector<MzTabDouble> v_rt;
            MzTabDoubleList rt_list;
            v_rt.emplace_back(id.rt);
            rt_list.set(v_rt);
            smsr.retention_time = rt_list;
            
            MzTabOptionalColumnEntry rank = make_pair("opt_global_rank", MzTabString(hit.rank));
            MzTabOptionalColumnEntry compoundId = make_pair("opt_global_compoundId", MzTabString(id.scan_index));
            MzTabOptionalColumnEntry compoundScanNumber = make_pair("opt_global_compoundScanNumber", MzTabString(id.scan_number));
            MzTabOptionalColumnEntry featureId = make_pair("opt_global_featureId", MzTabString(id.feature_id));
            MzTabOptionalColumnEntry native_id = make_pair("opt_global_native_id", MzTabString(id.native_id));

            smsr.opt_.push_back(rank);
            smsr.opt_.push_back(compoundId);
            smsr.opt_.push_back(compoundScanNumber);
            smsr.opt_.push_back(featureId);
            smsr.opt_.push_back(native_id);
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
