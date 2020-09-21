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
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/SYSTEM/File.h>
#include <boost/regex.hpp>

using namespace OpenMS;
using namespace std;

int SiriusMzTabWriter::extract_scan_index(const String &path)
{
  boost::regex regexp_ind("--(?<SCAN>\\d+)--");
  return SpectrumLookup::extractScanNumber(path, regexp_ind, false);;
}

std::map< String, Size > SiriusMzTabWriter::extract_columnname_to_columnindex(CsvFile& csvfile)
{
  StringList header_row;
  std::map< String, Size > columnname_to_columnindex;
  csvfile.getRow(0, header_row);

  for (size_t i = 0; i < header_row.size(); i++)
  {
    columnname_to_columnindex.insert(make_pair(header_row[i], i));
  }

  return columnname_to_columnindex;
};

void SiriusMzTabWriter::read(const std::vector<String>& sirius_output_paths,
                             const String& original_input_mzml,
                             const Size& top_n_hits,
                             MzTab& result)
{
  SiriusMzTabWriter::SiriusAdapterRun sirius_result;

  for (const auto& it : sirius_output_paths)
  {
    // extract mz, rt of the precursor and the nativeID of the corresponding MS2 spectra in the spectrum.ms file
    StringList ext_nid; // multiple possible MS2 spectra
    double ext_mz = 0.0;
    double ext_rt = 0.0;
    const String sirius_spectrum_ms = it + "/spectrum.ms";
    ifstream spectrum_ms_file(sirius_spectrum_ms);
    if (spectrum_ms_file)
    {
      const String nid_prefix = "##nid";
      const String rt_prefix = ">rt";
      const String pmass_prefix = ">parentmass";
      String line;
      while (getline(spectrum_ms_file, line))
      {
        if (line.hasPrefix(pmass_prefix))
        {
           ext_mz = String(line.erase(line.find(pmass_prefix), pmass_prefix.size())).toDouble();
        }
        else if (line.hasPrefix(rt_prefix))
        {
           line = line.erase(line.find("s"), 1); // >rt 418.39399999998s - remove unit "s"
           ext_rt = String(line.erase(line.find(rt_prefix), rt_prefix.size())).toDouble();
        }
        else if (line.hasPrefix(nid_prefix))
        {
           ext_nid.emplace_back(line.erase(line.find(nid_prefix), nid_prefix.size()));
        }
      }
      spectrum_ms_file.close();
    }

    // extract data from formula_candidates.csv
    const std::string pathtosiriuscsv = it + "/formula_candidates.tsv";

    ifstream file(pathtosiriuscsv);
    if (file) 
    {
      CsvFile compounds(pathtosiriuscsv, '\t');
      const UInt rowcount = compounds.rowCount();

      if (rowcount > 1)
      {
        // correction if the rowcount is smaller than the number of hits used as parameter
        // rowcount-1 because the csv header will be skipped in the loop later on.
        int header = 1;
        const UInt top_n_hits_cor = (top_n_hits >= rowcount) ? rowcount-header : top_n_hits;
        
        // fill identification structure containing all candidate hits for a single spectrum
        SiriusMzTabWriter::SiriusAdapterIdentification sirius_id;

        // extract scan_number from path
        OpenMS::String str = File::path(pathtosiriuscsv);
        int scan_index = SiriusMzTabWriter::extract_scan_index(str);

        // extract scan_number from string
        boost::regex regexp("-(?<SCAN>\\d+)--");
        int scan_number = SpectrumLookup::extractScanNumber(str, regexp, false);

        // extract feature_id from string
        boost::smatch match;
        String feature_id;
        boost::regex regexp_feature("_(?<SCAN>\\d+)-");

        bool found = boost::regex_search(str, match, regexp_feature);
        if (found && match["SCAN"].matched) {feature_id = "id_" + match["SCAN"].str();}
        String unassigned = "null";

        std::map< String, Size > columnname_to_columnindex = SiriusMzTabWriter::extract_columnname_to_columnindex(compounds);

        // formula	adduct	precursorFormula	rank	rankingScore	IsotopeScore	TreeScore	siriusScore	explainedPeaks	explainedIntensity
        // j = 1 because of .csv file format (header)
        for (Size j = 1; j <= top_n_hits_cor; ++j)
        {
          StringList sl;
          compounds.getRow(j, sl);
          SiriusMzTabWriter::SiriusAdapterHit sirius_hit;

          // maybe should check columnname instead?
          // rank	molecularFormula	adduct	precursorFormula	rankingScore	TreeIsotope_Score	Tree_Score	Isotope_Score	explainedPeaks	explainedIntensity
          // parse single candidate hit
          sirius_hit.formula = sl[columnname_to_columnindex.at("molecularFormula")];
          sirius_hit.adduct = sl[columnname_to_columnindex.at("adduct")];
          sirius_hit.precursor_formula = sl[columnname_to_columnindex.at("precursorFormula")];
          sirius_hit.rank = sl[columnname_to_columnindex.at("rank")].toInt();
          // sirius_hit.rank_score = sl[columnname_to_columnindex.at("rankingScore")].toDouble(); // removed SIRIUS 4.4.7
          sirius_hit.iso_score = sl[columnname_to_columnindex.at("Isotope_Score")].toDouble();
          sirius_hit.tree_score = sl[columnname_to_columnindex.at("Tree_Score")].toDouble();
          sirius_hit.sirius_score = sl[columnname_to_columnindex.at("TreeIsotope_Score")].toDouble();
          sirius_hit.explainedpeaks = sl[columnname_to_columnindex.at("explainedPeaks")].toInt();
          sirius_hit.explainedintensity = sl[columnname_to_columnindex.at("explainedIntensity")].toDouble();

          sirius_id.hits.push_back(sirius_hit);
        }

        sirius_id.mz = ext_mz;
        sirius_id.rt = ext_rt;
        sirius_id.native_ids = ext_nid;
        sirius_id.scan_index = scan_index;
        sirius_id.scan_number = scan_number;
        // check if results were assigned to a feature
        if (feature_id != "id_0")
        {
          sirius_id.feature_id = feature_id;
        }
        else
        {
          sirius_id.feature_id = unassigned;
        }
        sirius_result.identifications.push_back(sirius_id);

        // write metadata to mzTab file
        MzTabMetaData md;
        MzTabMSRunMetaData md_run;
        md_run.location = MzTabString(original_input_mzml);
        md.ms_run[1] = md_run;
        md.description = MzTabString("Sirius-4.4.29");

        //needed for header generation (score)
        std::map<Size, MzTabParameter> smallmolecule_search_engine_score;
        smallmolecule_search_engine_score[1].setName("sirius_score");
        smallmolecule_search_engine_score[2].setName("tree_score");
        smallmolecule_search_engine_score[3].setName("iso_score");
        //smallmolecule_search_engine_score[4].setName("rank_score");
        md.smallmolecule_search_engine_score = smallmolecule_search_engine_score;
        result.setMetaData(md);

        // write results to mzTab file
        MzTabSmallMoleculeSectionRows smsd;
        for (Size i = 0; i < sirius_result.identifications.size(); ++i)
        {
          const SiriusMzTabWriter::SiriusAdapterIdentification &id = sirius_result.identifications[i];
          for (Size j = 0; j < id.hits.size(); ++j)
          {
            const SiriusMzTabWriter::SiriusAdapterHit &hit = id.hits[j];
            MzTabSmallMoleculeSectionRow smsr;

            map<Size, MzTabDouble> engine_score;
            engine_score[1] = MzTabDouble(hit.sirius_score);
            engine_score[2] = MzTabDouble(hit.tree_score);
            engine_score[3] = MzTabDouble(hit.iso_score);
            //engine_score[4] = MzTabDouble(hit.rank_score);
            smsr.best_search_engine_score = engine_score;

            smsr.chemical_formula = MzTabString(hit.formula);
            smsr.exp_mass_to_charge = MzTabDouble(id.mz);
           
            vector<MzTabDouble> v_rt;
            MzTabDoubleList rt_list;
            v_rt.emplace_back(id.rt); 
            rt_list.set(v_rt);
            smsr.retention_time = rt_list;
            
            MzTabOptionalColumnEntry adduct = make_pair("opt_global_adduct", MzTabString(hit.adduct));
            MzTabOptionalColumnEntry precursor_formula = make_pair("opt_gobal_precursorFormula", MzTabString(hit.precursor_formula));
            MzTabOptionalColumnEntry rank = make_pair("opt_global_rank", MzTabString(hit.rank));
            MzTabOptionalColumnEntry explainedPeaks = make_pair("opt_global_explainedPeaks", MzTabString(hit.explainedpeaks));
            MzTabOptionalColumnEntry explainedIntensity = make_pair("opt_global_explainedIntensity", MzTabString(hit.explainedintensity));
            MzTabOptionalColumnEntry compoundId = make_pair("opt_global_compoundId", MzTabString(id.scan_index));
            MzTabOptionalColumnEntry compoundScanNumber = make_pair("opt_global_compoundScanNumber", MzTabString(id.scan_number));
            MzTabOptionalColumnEntry featureId = make_pair("opt_global_featureId", MzTabString(id.feature_id));

            vector<MzTabString> m_native_ids;
            MzTabStringList ml_native_ids;
            ml_native_ids.setSeparator('|');
            for (auto& element : id.native_ids)
            {
              m_native_ids.emplace_back(MzTabString(element));
            }
            ml_native_ids.set(m_native_ids);

            MzTabOptionalColumnEntry native_ids = make_pair("opt_global_native_id", MzTabString(ml_native_ids.toCellString()));

            smsr.opt_.push_back(adduct);
            smsr.opt_.push_back(precursor_formula);
            smsr.opt_.push_back(rank);
            smsr.opt_.push_back(explainedPeaks);
            smsr.opt_.push_back(explainedIntensity);
            smsr.opt_.push_back(compoundId);
            smsr.opt_.push_back(compoundScanNumber);
            smsr.opt_.push_back(featureId);
            smsr.opt_.push_back(native_ids);
            smsd.push_back(smsr);
          }
        }  
        result.setSmallMoleculeSectionRows(smsd);
      }
      file.close();
    }
  }
} // namespace OpenMS

/// @endcond
