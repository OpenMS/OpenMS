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
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <fstream>
#include <QtCore/QDir>
#include <QtCore/QString>

using namespace std;

namespace OpenMS
{

  std::vector<SiriusFragmentAnnotation::SiriusTargetDecoySpectra> SiriusFragmentAnnotation::extractAndResolveSiriusAnnotations(const std::vector<String>& sirius_workspace_subdirs, const double& score_threshold, bool use_exact_mass)
  {
    std::map<String, SiriusFragmentAnnotation::SiriusTargetDecoySpectra> native_ids_annotated_spectra;
    std::vector<SiriusFragmentAnnotation::SiriusTargetDecoySpectra> annotated_spectra;
    for (const auto& subdir : sirius_workspace_subdirs)
    {
      MSSpectrum annotated_spectrum;
      MSSpectrum annotated_decoy;

      SiriusFragmentAnnotation::extractSiriusFragmentAnnotationMapping(subdir,
                                                                       annotated_spectrum,
                                                                       use_exact_mass);

      SiriusFragmentAnnotation::extractSiriusDecoyAnnotationMapping(subdir,
                                                                    annotated_decoy);

      // if no target was generated - no decoy will be available (due to SIRIUS)
      if (annotated_spectrum.empty())
      {
        continue;
      }
      else
      {
        // only use spectra over a certain score threshold (0-1)
        if (double(annotated_spectrum.getMetaValue("score")) >= score_threshold)
        {
          // resolve multiple use of the same concatenated nativeids based on the sirius score (used for multiple features/identifications)
          map<String, SiriusFragmentAnnotation::SiriusTargetDecoySpectra>::iterator it;
          it = native_ids_annotated_spectra.find(annotated_spectrum.getNativeID());
          if (it != native_ids_annotated_spectra.end())
          {
            if (double(annotated_spectrum.getMetaValue("score")) >= double(it->second.target.getMetaValue("score")))
            {
              SiriusFragmentAnnotation::SiriusTargetDecoySpectra target_decoy(annotated_spectrum, annotated_decoy);
              native_ids_annotated_spectra.insert(make_pair(annotated_spectrum.getNativeID(), target_decoy));
            }
          }
          else
          {
            SiriusFragmentAnnotation::SiriusTargetDecoySpectra target_decoy(annotated_spectrum, annotated_decoy);
            native_ids_annotated_spectra.insert(make_pair(annotated_spectrum.getNativeID(), target_decoy));
          }
        }
      }
    }

    // convert to vector
    for (const auto& it : native_ids_annotated_spectra)
    {
      annotated_spectra.emplace_back(std::move(it.second));
    }

    return annotated_spectra;
  }

  void SiriusFragmentAnnotation::extractSiriusFragmentAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass)
  {
    OpenMS::String concat_native_ids = SiriusFragmentAnnotation::extractConcatNativeIDsFromSiriusMS_(path_to_sirius_workspace);
    OpenMS::String concat_m_ids = SiriusFragmentAnnotation::extractConcatMIDsFromSiriusMS_(path_to_sirius_workspace);
    SiriusFragmentAnnotation::extractAnnotationFromSiriusFile_(path_to_sirius_workspace, msspectrum_to_fill, use_exact_mass);
    msspectrum_to_fill.setNativeID(concat_native_ids);
    msspectrum_to_fill.setName(concat_m_ids);
  }

  void SiriusFragmentAnnotation::extractSiriusDecoyAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill)
  {
    OpenMS::String concat_native_ids = SiriusFragmentAnnotation::extractConcatNativeIDsFromSiriusMS_(path_to_sirius_workspace);
    OpenMS::String concat_m_ids = SiriusFragmentAnnotation::extractConcatMIDsFromSiriusMS_(path_to_sirius_workspace);
    SiriusFragmentAnnotation::extractAnnotationFromDecoyFile_(path_to_sirius_workspace, msspectrum_to_fill);
    msspectrum_to_fill.setNativeID(concat_native_ids);
    msspectrum_to_fill.setName(concat_m_ids);
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
  
  // use the first ranked sumformula (works for known and known_unkowns)
  // currently only supports sumformula from rank 1
  void SiriusFragmentAnnotation::extractAnnotationFromSiriusFile_(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass)
  { 
    if (!msspectrum_to_fill.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Non empty MSSpectrum was provided");
    }

    const std::string sirius_spectra_dir = path_to_sirius_workspace + "/spectra/";
    QDir dir(QString::fromStdString(sirius_spectra_dir));
    if (dir.exists())
    {
      std::map< Size, String > rank_filename = SiriusFragmentAnnotation::extractCompoundRankingAndFilename_(path_to_sirius_workspace);
      std::map< Size, double > rank_score = SiriusFragmentAnnotation::extractCompoundRankingAndScore_(path_to_sirius_workspace);

      if (rank_filename.empty() || rank_score.empty())
      {
        OPENMS_LOG_WARN << "Extraction of the compound ranking, filename and score failed for, please check if the SIRIUS project space is correct for." << sirius_spectra_dir << std::endl;
      }

      // use first file in folder (rank 1)
      String filename = rank_filename.at(1); // rank 1
      double score = rank_score.at(1); // rank 1
      QFileInfo firstfile(dir,filename.toQString());

      if (use_exact_mass)
      {
        msspectrum_to_fill.setMetaValue("peak_mz", DataValue("exact_mass"));
      }
      else
      {
        msspectrum_to_fill.setMetaValue("peak_mz", DataValue("mz"));
      }

      // filename: sumformula_adduct.csv - save sumformula and adduct as metavalue
      String current_sumformula = filename.substr(0, filename.find_last_of("_"));
      String current_adduct = filename.substr(filename.find_last_of("_") + 1, filename.find_last_of(".") - filename.find_last_of("_") - 1);
      msspectrum_to_fill.setMetaValue("annotated_sumformula", DataValue(current_sumformula));
      msspectrum_to_fill.setMetaValue("annotated_adduct", DataValue(current_adduct));
      msspectrum_to_fill.setMetaValue("decoy",0);

      // read file and save in MSSpectrum
      ifstream fragment_annotation_file(firstfile.absoluteFilePath().toStdString());
      if (fragment_annotation_file)
      {
        //mz  intensity   rel.intensity   exactmass   formula ionization
        //51.023137   713.15  9.36    51.022927   C4H2    [M + H]+

        std::vector<Peak1D> fragments_mzs_ints;
        MSSpectrum::FloatDataArray fragments_exactmasses;
        MSSpectrum::StringDataArray fragments_explanations;
        MSSpectrum::StringDataArray fragments_ionization;
        if (use_exact_mass)
        {
          fragments_exactmasses.setName("mz");
        }
        else
        {
          fragments_exactmasses.setName("exact_mass");
        }
        fragments_explanations.setName("explanation");
        String line;
        std::getline(fragment_annotation_file, line); // skip header
        while (std::getline(fragment_annotation_file, line))
        {
          Peak1D fragment_mz_int;
          StringList splitted_line;
          line.split("\t",splitted_line);     
          // option to use the exact mass as peakMZ (e.g. for library preparation).
          if (use_exact_mass)
          {
            fragment_mz_int.setMZ(splitted_line[3].toDouble());
            fragments_exactmasses.push_back(splitted_line[0].toDouble());
          }
          else
          {
            fragment_mz_int.setMZ(splitted_line[0].toDouble());
            fragments_exactmasses.push_back(splitted_line[3].toDouble());
          }
          fragment_mz_int.setIntensity(splitted_line[1].toDouble());
          fragments_mzs_ints.push_back(fragment_mz_int);
          fragments_explanations.push_back(splitted_line[4]);
          fragments_ionization.push_back(splitted_line[5]);
        }
        msspectrum_to_fill.setMSLevel(2);
        msspectrum_to_fill.setMetaValue("score", DataValue(score));
        msspectrum_to_fill.insert(msspectrum_to_fill.begin(), fragments_mzs_ints.begin(), fragments_mzs_ints.end());
        msspectrum_to_fill.getFloatDataArrays().push_back(fragments_exactmasses);
        msspectrum_to_fill.getStringDataArrays().push_back(fragments_explanations);
        msspectrum_to_fill.getStringDataArrays().push_back(fragments_ionization);
      } 
    }
    else
    {
      OPENMS_LOG_DEBUG << "Directory 'spectra' was not found for: " << sirius_spectra_dir << std::endl;
    }
  }

  // use the first ranked sumformula (works for known and known_unkowns)
  // currently only supports sumformula from rank 1
  void SiriusFragmentAnnotation::extractAnnotationFromDecoyFile_(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill)
  {
    if (!msspectrum_to_fill.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Non empty MSSpectrum was provided");
    }

    const std::string sirius_spectra_dir = path_to_sirius_workspace + "/decoys/";
    QDir dir(QString::fromStdString(sirius_spectra_dir));
    if (dir.exists())
    {
      std::map< Size, String > rank_filename = SiriusFragmentAnnotation::extractCompoundRankingAndFilename_(path_to_sirius_workspace);
      std::map< Size, double > rank_score = SiriusFragmentAnnotation::extractCompoundRankingAndScore_(path_to_sirius_workspace);

      if (rank_filename.empty() || rank_score.empty())
      {
        OPENMS_LOG_WARN << "Extraction of the compound ranking, filename and score failed for, please check if the SIRIUS project space is correct for." << sirius_spectra_dir << std::endl;
      }

      // use first file in folder (rank 1)
      String filename = rank_filename.at(1); // rank 1
      double score = rank_score.at(1); // rank 1
      QFileInfo firstfile(dir,filename.toQString());

      msspectrum_to_fill.setMetaValue("peak_mz", DataValue("mz"));

      // filename: sumformula_adduct.csv - save sumformula and adduct as metavalue
      String current_sumformula = filename.substr(0, filename.find_last_of("_"));
      String current_adduct = filename.substr(filename.find_last_of("_") + 1, filename.find_last_of(".") - filename.find_last_of("_") - 1);
      msspectrum_to_fill.setMetaValue("annotated_sumformula", DataValue(current_sumformula));
      msspectrum_to_fill.setMetaValue("annotated_adduct", DataValue(current_adduct));
      msspectrum_to_fill.setMetaValue("decoy",1);

      // read file and save in MSSpectrum
      ifstream fragment_annotation_file(firstfile.absoluteFilePath().toStdString());
      if (fragment_annotation_file)
      {
        //mz	      rel.intensity	formula	ionization
        //46.994998	0.71	        CH2S	  [M + H]+

        std::vector<Peak1D> fragments_mzs_ints;
        MSSpectrum::StringDataArray fragments_explanations;
        MSSpectrum::StringDataArray fragments_ionization;

        fragments_explanations.setName("explanation");
        String line;
        std::getline(fragment_annotation_file, line); // skip header
        while (std::getline(fragment_annotation_file, line))
        {
          Peak1D fragment_mz_int;
          StringList splitted_line;
          line.split("\t",splitted_line);
          fragment_mz_int.setMZ(splitted_line[0].toDouble());
          fragment_mz_int.setIntensity(splitted_line[1].toDouble());
          fragments_mzs_ints.push_back(fragment_mz_int);
          fragments_explanations.push_back(splitted_line[2]);
          fragments_ionization.push_back(splitted_line[3]);
        }
        msspectrum_to_fill.setMSLevel(2);
        msspectrum_to_fill.setMetaValue("score", DataValue(score));
        msspectrum_to_fill.insert(msspectrum_to_fill.begin(), fragments_mzs_ints.begin(), fragments_mzs_ints.end());
        msspectrum_to_fill.getStringDataArrays().push_back(fragments_explanations);
        msspectrum_to_fill.getStringDataArrays().push_back(fragments_ionization);
      }
    }
    else
    {
      OPENMS_LOG_DEBUG << "Directory 'decoys' was not found for: " << sirius_spectra_dir << std::endl;
    }
  }

} // namespace OpenMS

/// @endcond
