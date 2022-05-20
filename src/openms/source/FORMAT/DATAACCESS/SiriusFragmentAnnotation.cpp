// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <fstream>
#include <QtCore/QDir>
#include <QtCore/QString>

using namespace std;

namespace OpenMS
{

  std::vector<MSSpectrum> SiriusFragmentAnnotation::extractSiriusAnnotationsTgtOnly
    (const std::vector<String>& sirius_workspace_subdirs, double score_threshold, bool use_exact_mass, bool resolve = true)
  {
    Size max_rank = 10000; // this should be enough for any search
    if (resolve) max_rank = 1;
    std::unordered_map<String, MSSpectrum> native_ids_annotated_spectra;
    std::vector<MSSpectrum> annotated_spectra;
    double score = 0.0;
    // There is one subdir for every candidate that we pass to Sirius
    // Currently it is not possible to run Sirius with multiple candidates, see https://github.com/OpenMS/OpenMS/issues/5882
    for (const auto& subdir : sirius_workspace_subdirs)
    {
      std::vector<MSSpectrum> best_annotated_spectra = extractAnnotationsFromSiriusFile(subdir, max_rank, false, use_exact_mass);

      if (!resolve) annotated_spectra.reserve(sirius_workspace_subdirs.size());

      for (auto& spectrum : best_annotated_spectra)
      {
        score = double(spectrum.getMetaValue(Constants::UserParam::SIRIUS_SCORE));
        // only use spectra over a certain score threshold (0-1)
        if (score >= score_threshold)
        {
          if (resolve)
          {
            // resolve multiple use of the same concatenated native ids based on the sirius score (used for multiple features/identifications)
            unordered_map<String, MSSpectrum>::iterator it;
            it = native_ids_annotated_spectra.find(spectrum.getNativeID());
            if (it != native_ids_annotated_spectra.end())
            {
              if (score >= double(it->second.getMetaValue(Constants::UserParam::SIRIUS_SCORE)))
              {
                it->second = spectrum;
              }
            }
            else
            {
              native_ids_annotated_spectra.insert(make_pair(spectrum.getNativeID(), spectrum));
            }
          }
          else
          {
            annotated_spectra.push_back(std::move(spectrum));
          }
        }
      }
    }

    if (resolve)
    {
      // convert temporary map to vector
      annotated_spectra.reserve(native_ids_annotated_spectra.size());
      for (auto& id_spec : native_ids_annotated_spectra)
      {
        annotated_spectra.emplace_back(std::move(id_spec.second));
      }
    }
    else
    {
      annotated_spectra.shrink_to_fit();
    }

    return annotated_spectra;
  }

  std::vector<SiriusFragmentAnnotation::SiriusTargetDecoySpectra> SiriusFragmentAnnotation::extractAndResolveSiriusAnnotations(
    const std::vector<String>& sirius_workspace_subdirs, double score_threshold, bool use_exact_mass)
  {
    std::map<String, SiriusFragmentAnnotation::SiriusTargetDecoySpectra> native_ids_annotated_spectra;
    std::vector<SiriusFragmentAnnotation::SiriusTargetDecoySpectra> annotated_spectra;
    MSSpectrum best_annotated_spectrum;
    double score = 0.0;
    for (const auto& subdir : sirius_workspace_subdirs)
    {
      std::vector<MSSpectrum> ann_spec_tmp = extractAnnotationsFromSiriusFile(subdir, 1, false, use_exact_mass);
      if (ann_spec_tmp.empty())
      {
        continue;
      }
      else
      {
        // max_rank 1 will get the best.
        best_annotated_spectrum = extractAnnotationsFromSiriusFile(subdir, 1, false, use_exact_mass)[0];
      }

      ann_spec_tmp = extractAnnotationsFromSiriusFile(subdir, 1, true, use_exact_mass);
      // if no spectrum can be extracted we add an empty spectrum
      //  to the TD pair for backwards compatibility with AssayGeneratorMetabo
      MSSpectrum annotated_decoy_for_best_tgt = MSSpectrum();
      if (!ann_spec_tmp.empty())
      {
        //I ASSUME that decoys are in the same ranking order as their corresponding targets.
        annotated_decoy_for_best_tgt = ann_spec_tmp[0];
      }
      else // fill with basics for backwards-compatibility. IMHO this should be solved differently.
      {
        OpenMS::String concat_native_ids = SiriusFragmentAnnotation::extractConcatNativeIDsFromSiriusMS_(subdir);
        OpenMS::String concat_m_ids = SiriusFragmentAnnotation::extractConcatMIDsFromSiriusMS_(subdir);
        annotated_decoy_for_best_tgt.setNativeID(concat_native_ids);
        annotated_decoy_for_best_tgt.setName(concat_m_ids);
      }

      score = double(best_annotated_spectrum.getMetaValue(Constants::UserParam::SIRIUS_SCORE));
      // only use spectra over a certain score threshold (0-1)
      if (score >= score_threshold)
      {
        // resolve multiple use of the same concatenated nativeids based on the sirius score (used for multiple features/identifications)
        map<String, SiriusFragmentAnnotation::SiriusTargetDecoySpectra>::iterator it;
        it = native_ids_annotated_spectra.find(best_annotated_spectrum.getNativeID());
        if (it != native_ids_annotated_spectra.end())
        {
          if (score >= double(it->second.target.getMetaValue(Constants::UserParam::SIRIUS_SCORE)))
          {
            SiriusFragmentAnnotation::SiriusTargetDecoySpectra target_decoy(best_annotated_spectrum, annotated_decoy_for_best_tgt);
            it->second = target_decoy;
          }
        }
        else
        {
          SiriusFragmentAnnotation::SiriusTargetDecoySpectra target_decoy(best_annotated_spectrum, annotated_decoy_for_best_tgt);
          native_ids_annotated_spectra.insert(make_pair(best_annotated_spectrum.getNativeID(), target_decoy));
        }
      }

    }

    // convert to vector
    annotated_spectra.reserve(native_ids_annotated_spectra.size());
    for (auto& id_spec : native_ids_annotated_spectra)
    {
      annotated_spectra.emplace_back(std::move(id_spec.second));
    }

    return annotated_spectra;
  }
  
  // extract native id from SIRIUS spectrum.ms output file (workspace - compound specific)
  // first native id in the spectrum.ms
  OpenMS::String SiriusFragmentAnnotation::extractConcatNativeIDsFromSiriusMS_(const String& path_to_sirius_workspace)
  {
    vector< String > ext_n_ids;
    String ext_n_id;
    const String sirius_spectrum_ms = path_to_sirius_workspace + "/spectrum.ms";
    ifstream spectrum_ms_file(sirius_spectrum_ms);
    if (spectrum_ms_file)
    {
      const OpenMS::String n_id_prefix = "##n_id ";
      String line;
      while (getline(spectrum_ms_file, line))
      {
        if (line.hasPrefix(n_id_prefix))
        {
          String n_id = line.erase(line.find(n_id_prefix), n_id_prefix.size());
          ext_n_ids.emplace_back(n_id);
        }
        else if (spectrum_ms_file.eof())
        {
          OPENMS_LOG_WARN << "No native id was found - please check your input mzML. " << std::endl;
          break;
        }
      }
      spectrum_ms_file.close();
    }
    ext_n_id = ListUtils::concatenate(ext_n_ids, "|");
    return ext_n_id;
  }

  // extract native id from SIRIUS spectrum.ms output file (workspace - compound specific)
  // first native id in the spectrum.ms
  OpenMS::String SiriusFragmentAnnotation::extractFeatureIDFromSiriusMS_(const String& path_to_sirius_workspace)
  {
    String fid = "";
    const String sirius_spectrum_ms = path_to_sirius_workspace + "/spectrum.ms";
    ifstream spectrum_ms_file(sirius_spectrum_ms);
    if (spectrum_ms_file)
    {
      const OpenMS::String fid_prefix = "##fid ";
      String line;
      while (getline(spectrum_ms_file, line))
      {
        if (line.hasPrefix(fid_prefix))
        {
          fid = line.erase(line.find(fid_prefix), fid_prefix.size());
          break; // one file can only have one fid
        }
        else if (spectrum_ms_file.eof())
        {
          return ""; // fid is optional and only there if the origin of the sirius run was an OpenMS feature (not a single MS2).
        }
      }
      spectrum_ms_file.close();
    }
    return fid;
  }

  // extract m_id from SIRIUS spectrum.ms output file (workspace - compound specific)
  // and concatenates them if multiple spectra have been used in case of the compound.
  OpenMS::String SiriusFragmentAnnotation::extractConcatMIDsFromSiriusMS_(const String& path_to_sirius_workspace)
  {
    vector< String > ext_m_ids;
    String ext_m_id;
    const String sirius_spectrum_ms = path_to_sirius_workspace + "/spectrum.ms";
    ifstream spectrum_ms_file(sirius_spectrum_ms);
    if (spectrum_ms_file)
    {
      const OpenMS::String m_id_prefix = "##m_id ";
      String line;
      while (getline(spectrum_ms_file, line))
      {
        if (line.hasPrefix(m_id_prefix))
        {
          String m_id = line.erase(line.find(m_id_prefix), m_id_prefix.size());
          ext_m_ids.emplace_back(m_id);
        }
        else if (spectrum_ms_file.eof())
        {
          OPENMS_LOG_WARN << "No SiriusAdapter m_id was found - please check your input mzML. " << std::endl;
          break;
        }
      }
      spectrum_ms_file.close();
    }
    ext_m_id = ListUtils::concatenate(ext_m_ids, "|");
    return ext_m_id;
  }

  // provides a mapping of rank and the file it belongs to since this is not encoded in the directory structure/filename
  std::map< Size, String > SiriusFragmentAnnotation::extractCompoundRankingAndFilename_(const String& path_to_sirius_workspace)
  {
    map< Size, String > rank_filename;
    String line;

    const String sirius_formula_candidates = path_to_sirius_workspace + "/formula_candidates.tsv"; // based on SIRIUS annotation
    ifstream fcandidates(sirius_formula_candidates);
    if (fcandidates)
    {
      CsvFile candidates(sirius_formula_candidates, '\t');
      const UInt rowcount = candidates.rowCount();

      std::map< std::string, Size > columnname_to_columnindex = SiriusMzTabWriter::extract_columnname_to_columnindex(candidates);

      // i starts at 1, due to header
      for (size_t i = 1; i < rowcount; i++)
      {
        StringList sl;
        candidates.getRow(i, sl);
        String adduct = sl[columnname_to_columnindex.at("adduct")];
        adduct.erase(std::remove_if(adduct.begin(), adduct.end(), ::isspace), adduct.end());
        rank_filename.emplace(std::make_pair(sl[columnname_to_columnindex.at("rank")].toInt(),
                              String(sl[columnname_to_columnindex.at("molecularFormula")] + "_" + adduct + ".tsv")));
      }
    }
    fcandidates.close();
    return rank_filename;
  }

  // provides a mapping of rank and the TreeIsotope_Score which can be used to resolve ambiguities in mapping
  // 43_Cyazofamid_[M+H]+_3 sample=1 period=1 cycle=683 experiment=3|sample=1 period=1 cycle=684 experiment=3|sample=1 period=1 cycle=685 experiment=5
  // 44_Ethofumesate_[M+K]+_3 sample=1 period=1 cycle=683 experiment=3|sample=1 period=1 cycle=684 experiment=3|sample=1 period=1 cycle=685 experiment=5
  // which have the same mass within 25 ppm due to their adduct in AccurateMassSearch
  std::map< Size, double > SiriusFragmentAnnotation::extractCompoundRankingAndScore_(const String& path_to_sirius_workspace)
  {
    map< Size, double > rank_score;
    String line;

    const String sirius_formula_candidates = path_to_sirius_workspace + "/formula_candidates.tsv"; // based on SIRIUS annotation
    ifstream fcandidates(sirius_formula_candidates);
    if (fcandidates)
    {
      CsvFile candidates(sirius_formula_candidates, '\t');
      const UInt rowcount = candidates.rowCount();

      std::map< std::string, Size > columnname_to_columnindex = SiriusMzTabWriter::extract_columnname_to_columnindex(candidates);

      // i starts at 1, due to header
      for (size_t i = 1; i < rowcount; i++)
      {
        StringList sl;
        candidates.getRow(i, sl);
        rank_score.emplace(std::make_pair(sl[columnname_to_columnindex.at("rank")].toInt(),
                                             sl[columnname_to_columnindex.at("explainedIntensity")].toDouble()));
      }
    }
    fcandidates.close();
    return rank_score;
  }

  std::vector<MSSpectrum> SiriusFragmentAnnotation::extractAnnotationsFromSiriusFile(const String& path_to_sirius_workspace, Size max_rank, bool decoy, bool use_exact_mass)
  {
    if (decoy) use_exact_mass = false;
    std::vector<MSSpectrum> result;
    std::string subfolder = decoy ? "/decoys/" : "/spectra/";
    const std::string sirius_spectra_dir = path_to_sirius_workspace + subfolder;
    QDir dir(QString::fromStdString(sirius_spectra_dir));

    if (dir.exists())
    {
      // TODO this probably can and should be extracted in one go.
      OpenMS::String concat_native_ids = SiriusFragmentAnnotation::extractConcatNativeIDsFromSiriusMS_(path_to_sirius_workspace);
      OpenMS::String concat_m_ids = SiriusFragmentAnnotation::extractConcatMIDsFromSiriusMS_(path_to_sirius_workspace);
      OpenMS::String fid = SiriusFragmentAnnotation::extractFeatureIDFromSiriusMS_(path_to_sirius_workspace);
      std::map< Size, String > rank_filename = SiriusFragmentAnnotation::extractCompoundRankingAndFilename_(path_to_sirius_workspace);
      std::map< Size, double > rank_score = SiriusFragmentAnnotation::extractCompoundRankingAndScore_(path_to_sirius_workspace);
      Size max_found_rank = rank_filename.rbegin()->first;
      if (rank_filename.empty() || rank_score.empty())
      {
        OPENMS_LOG_WARN << "Extraction of the compound ranking, filename and score failed for, please check if the SIRIUS project space is correct for." << sirius_spectra_dir << std::endl;
      }
      else
      {
        result.resize(std::min(max_rank, max_found_rank));
      }

      String suffix = "";
      for (unsigned int i = 1; i <= result.size(); ++i)
      {
        MSSpectrum& msspectrum_to_fill = result[i-1];
        // to be backwards compatible, the best rank does not have a suffix
        // TODO should be changed. But I need to figure out where.
        if (i > 1)
        {
          suffix = "_" + String(i);
        }
        msspectrum_to_fill.setNativeID(concat_native_ids + suffix);
        msspectrum_to_fill.setName(concat_m_ids + suffix);
        String filename = rank_filename.at(i); // rank 1
        double score = rank_score.at(i); // rank 1
        QFileInfo sirius_result_file(dir,filename.toQString());

        if (use_exact_mass)
        {
          msspectrum_to_fill.setMetaValue(Constants::UserParam::SIRIUS_PEAKMZ, DataValue(Constants::UserParam::SIRIUS_EXACTMASS));
        }
        else
        {
          msspectrum_to_fill.setMetaValue(Constants::UserParam::SIRIUS_PEAKMZ, DataValue(Constants::UserParam::SIRIUS_MZ));
        }

        // filename: sumformula_adduct.csv - save sumformula and adduct as metavalue
        String current_sumformula = filename.substr(0, filename.find_last_of("_"));
        String current_adduct = filename.substr(filename.find_last_of("_") + 1, filename.find_last_of(".") - filename.find_last_of("_") - 1);
        msspectrum_to_fill.setMetaValue(Constants::UserParam::SIRIUS_ANNOTATED_SUMFORMULA, DataValue(current_sumformula));
        msspectrum_to_fill.setMetaValue(Constants::UserParam::SIRIUS_ANNOTATED_ADDUCT, DataValue(current_adduct));
        msspectrum_to_fill.setMetaValue(Constants::UserParam::SIRIUS_DECOY, decoy);
        if (!fid.empty())
        {
          msspectrum_to_fill.setMetaValue(Constants::UserParam::SIRIUS_FEATURE_ID, fid);
        }

        // read file and save in MSSpectrum
        ifstream fragment_annotation_file(sirius_result_file.absoluteFilePath().toStdString());
        if (fragment_annotation_file)
        {
          // Target schema
          //mz  intensity   rel.intensity   exactmass   formula ionization
          //51.023137   713.15  9.36    51.022927   C4H2    [M + H]+

          //Decoy schema
          //mz	      rel.intensity	formula	ionization
          //46.994998	0.71	        CH2S	  [M + H]+

          std::vector<Peak1D> fragments_mzs_ints;
          MSSpectrum::FloatDataArray fragments_extra_masses;
          MSSpectrum::StringDataArray fragments_explanations;
          MSSpectrum::StringDataArray fragments_ionization;
          Size main_mz_col = 0;
          Size extra_mz_col = 3;

          // option to use the exact mass as peak MZ (e.g. for library preparation).
          if (use_exact_mass)
          {
            main_mz_col = 3;
            extra_mz_col = 0;
            fragments_extra_masses.setName(Constants::UserParam::SIRIUS_MZ);
          }
          else
          {
            fragments_extra_masses.setName(Constants::UserParam::SIRIUS_EXACTMASS);
          }

          fragments_explanations.setName(Constants::UserParam::SIRIUS_EXPLANATION);
          String line;
          std::getline(fragment_annotation_file, line); // skip header
          while (std::getline(fragment_annotation_file, line))
          {
            StringList splitted_line;
            line.split("\t",splitted_line);

            if (!decoy)
            {
              fragments_extra_masses.push_back(splitted_line[extra_mz_col].toFloat());
            }
            fragments_mzs_ints.emplace_back(splitted_line[main_mz_col].toDouble(), splitted_line[1].toFloat());
            fragments_explanations.push_back(splitted_line[splitted_line.size() - 2]);
            fragments_ionization.push_back(splitted_line[splitted_line.size() - 1]);
          }
          msspectrum_to_fill.setMetaValue(Constants::UserParam::SIRIUS_SCORE, DataValue(score));

          msspectrum_to_fill.setMSLevel(2);
          msspectrum_to_fill.swap(fragments_mzs_ints);
          msspectrum_to_fill.getFloatDataArrays().push_back(fragments_extra_masses);
          msspectrum_to_fill.getStringDataArrays().push_back(fragments_explanations);
          msspectrum_to_fill.getStringDataArrays().push_back(fragments_ionization);
        }
      }
    }
    else
    {
      OPENMS_LOG_DEBUG << "Directory '" + subfolder + "' was not found for: " << sirius_spectra_dir << std::endl;
    }
    return result;
  }

} // namespace OpenMS

/// @endcond
