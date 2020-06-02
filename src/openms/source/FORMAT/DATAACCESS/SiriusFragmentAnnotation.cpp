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

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

  void SiriusFragmentAnnotation::extractSiriusFragmentAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass)
  {
    OpenMS::String native_id = SiriusFragmentAnnotation::extractNativeIDFromSiriusMS_(path_to_sirius_workspace);
    OpenMS::String mid = SiriusFragmentAnnotation::extractMIDFromSiriusMS_(path_to_sirius_workspace);
    SiriusFragmentAnnotation::extractAnnotationFromSiriusFile_(path_to_sirius_workspace, msspectrum_to_fill, use_exact_mass);
    msspectrum_to_fill.setNativeID(native_id);
    msspectrum_to_fill.setName(mid);
  }
  
  // extract native id from SIRIUS spectrum.ms output file (workspace - compound specific)
  // first native id in the spectrum.ms
  OpenMS::String SiriusFragmentAnnotation::extractNativeIDFromSiriusMS_(const String& path_to_sirius_workspace)
  {
    vector< String > ext_nids;
    String ext_nid;
    const String sirius_spectrum_ms = path_to_sirius_workspace + "/spectrum.ms";
    ifstream spectrum_ms_file(sirius_spectrum_ms);
    if (spectrum_ms_file)
    {
      const OpenMS::String nid_prefix = "##nid ";
      String line;
      while (getline(spectrum_ms_file, line))
      {
        if (line.hasPrefix(nid_prefix))
        {
           String nid = line.erase(line.find(nid_prefix), nid_prefix.size());
           ext_nids.emplace_back(nid);
        }
        else if (spectrum_ms_file.eof())
        {
           OPENMS_LOG_WARN << "No native id was found - please check your input mzML. " << std::endl;
           break;
        }
      }
      spectrum_ms_file.close();
    }
    ext_nid = ListUtils::concatenate(ext_nids, "|");
    return ext_nid;
  }

  // extract mid from SIRIUS spectrum.ms output file (workspace - compound specific)
  // first mid id in the spectrum.ms
  OpenMS::String SiriusFragmentAnnotation::extractMIDFromSiriusMS_(const String& path_to_sirius_workspace)
  {
    vector< String > ext_mids;
    String ext_mid;
    const String sirius_spectrum_ms = path_to_sirius_workspace + "/spectrum.ms";
    ifstream spectrum_ms_file(sirius_spectrum_ms);
    if (spectrum_ms_file)
    {
      const OpenMS::String mid_prefix = "##mid ";
      String line;
      while (getline(spectrum_ms_file, line))
      {
        if (line.hasPrefix(mid_prefix))
        {
          String mid = line.erase(line.find(mid_prefix), mid_prefix.size());
          ext_mids.emplace_back(mid);
        }
        else if (spectrum_ms_file.eof())
        {
          OPENMS_LOG_WARN << "No SiriusAdapter_MID id was found - please check your input mzML. " << std::endl;
          break;
        }
      }
      spectrum_ms_file.close();
    }
    ext_mid = ListUtils::concatenate(ext_mids, "|");
    return ext_mid;
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

      std::map< String, Size > columnname_to_columnindex = SiriusMzTabWriter::extract_columnname_to_columnindex(candidates);

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
  
  // use the first ranked sumformula (works for known and known_unkowns)
  // currently only supports sumformula from rank 1
  void SiriusFragmentAnnotation::extractAnnotationFromSiriusFile_(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass)
  { 
    if (!msspectrum_to_fill.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Non empty MSSpectrum was provided");
    }

    std::map< Size, String > rank_filename = SiriusFragmentAnnotation::extractCompoundRankingAndFilename_(path_to_sirius_workspace);

    const std::string sirius_spectra_dir = path_to_sirius_workspace + "/spectra/";
    QDir dir(QString::fromStdString(sirius_spectra_dir));
    if (dir.exists())
    {
      if (use_exact_mass)
      {
        msspectrum_to_fill.setMetaValue("peak_mz", DataValue("exact_mass"));
      }
      else
      {
        msspectrum_to_fill.setMetaValue("peak_mz", DataValue("mz"));
      }

      // use first file in folder (rank 1)
      String filename = rank_filename.at(1); // rank 1
      QFileInfo firstfile(dir,filename.toQString());

      // filename: sumformula_adduct.csv - save sumformula and adduct as metavalue
      String current_sumformula = filename.substr(0, filename.find_last_of("_"));
      String current_adduct = filename.substr(filename.find_last_of("_") + 1, filename.find_last_of(".") - filename.find_last_of("_") - 1);
      msspectrum_to_fill.setMetaValue("annotated_sumformula", DataValue(current_sumformula));
      msspectrum_to_fill.setMetaValue("annotated_adduct", DataValue(current_adduct));

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
        msspectrum_to_fill.insert(msspectrum_to_fill.begin(), fragments_mzs_ints.begin(), fragments_mzs_ints.end());
        msspectrum_to_fill.getFloatDataArrays().push_back(fragments_exactmasses);
        msspectrum_to_fill.getStringDataArrays().push_back(fragments_explanations);
        msspectrum_to_fill.getStringDataArrays().push_back(fragments_ionization);
      } 
    }
    else
    {
      OPENMS_LOG_WARN << "Directory 'spectra' was not found for: " << sirius_spectra_dir << std::endl;
    }
  }
} // namespace OpenMS

/// @endcond
