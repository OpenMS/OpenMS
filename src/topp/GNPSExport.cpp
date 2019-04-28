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
// $Maintainer: UCSD Dorrestein Lab $
// $Authors: Abinesh Sarvepalli $
// --------------------------------------------------------------------------

//----------------------------------------------------------
// Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_GNPSExport GNPSExport
  @brief Export MS/MS data in .MGF format for GNPS (http://gnps.ucsd.edu).
GNPS (Global Natural Products Social Molecular Networking, http://gnps.ucsd.edu) is an open-access knowledge base for community-wide organisation and sharing of raw, processed or identified tandem mass (MS/MS) spectrometry data. The GNPS web-platform makes possible to perform spectral library search against public MS/MS spectral libraries, as well as to perform various data analysis such as MS/MS molecular networking, network annotation propagation (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006089), and the Dereplicator-based annotation (https://www.nature.com/articles/nchembio.2219). The GNPS manuscript is available here: https://www.nature.com/articles/nbt.3597

This tool was developed for the OpenMS-GNPS workflow. It can be accessed on GNPS (https://gnps.ucsd.edu/ProteoSAFe/static/gnps-experimental.jsp). The steps used by that workflow are as following:
  1. Input mzML files
  2. Run the FeatureFinderMetabo tool
  3. Run the IDMapper tool
  4. Run the MapAlignerPoseClustering tool
  5. Run the FeautureLinkerUnlabeledKD tool
  6. Run the FileConverter tool (output FeatureXML format)
  7. Run the MetaboliteAdductDecharger
  8. Run the GNPSExport to export an .MGF
  9. Upload the .MGF file on http://gnps.ucsd.edu and follow the instructions here:
https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <iostream>
#include <fstream>

using namespace OpenMS;
using namespace std;

class TOPPGNPSExport : public TOPPBase
{
public:
  TOPPGNPSExport() :
  TOPPBase("GNPSExport", "Tool to export consensus features into MGF format", false) {}

private:
  double DEF_COSINE_SIMILARITY = 0.95;
  double DEF_PRECURSOR_MZ_TOLERANCE = 0.0001;
  double DEF_PRECURSOR_RT_TOLERANCE = 5;

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in_cm", "<file>", "", "input file containing consensus elements with \'peptide\' annotations");
    setValidFormats_("in_cm", ListUtils::create<String>("consensusXML"));

    registerInputFileList_("in_mzml", "<files>", ListUtils::create<String>(""), "original mzml files containing ms/ms spectrum information");
    setValidFormats_("in_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Output MGF file");
    setValidFormats_("out", ListUtils::create<String>("mgf"));

    registerStringOption_("output_type", "<choice>", "full_spectra", "specificity of mgf output information", false);
    setValidStrings_("output_type", ListUtils::create<String>("full_spectra,merged_spectra"));

    registerDoubleOption_("precursor_mz_tolerance", "<num>", DEF_PRECURSOR_MZ_TOLERANCE, "Tolerance mz window for precursor selection", false);
    registerDoubleOption_("precursor_rt_tolerance", "<num>", DEF_PRECURSOR_RT_TOLERANCE, "Tolerance rt window for precursor selection", false);

    registerTOPPSubsection_("merged_spectra", "Options for exporting mgf file with merged spectra per feature");
    registerDoubleOption_("merged_spectra:cos_similarity", "<num>", DEF_COSINE_SIMILARITY, "Cosine similarity threshold for merged_spectra output", false);
  }

  // the main function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
    ProgressLogger progress_logger;
    progress_logger.setLogType(log_type_);

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String consensus_file_path(getStringOption_("in_cm"));
    StringList mzml_file_paths(getStringList_("in_mzml"));

    String out(getStringOption_("out"));

    String output_type(getStringOption_("output_type"));

    double prec_mz_tol(getDoubleOption_("precursor_mz_tolerance"));
    double prec_rt_tol(getDoubleOption_("precursor_rt_tolerance"));

    double cos_sim(getDoubleOption_("merged_spectra:cos_similarity"));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    // ConsensusMap
    ConsensusXMLFile consensus_file;
    consensus_file.setLogType(log_type_);
    ConsensusMap consensus_map;
    consensus_file.load(consensus_file_path, consensus_map);

    // MSExperiment
    vector<MSExperiment> ms_maps;
    for (const auto& mzml_file_path : mzml_file_paths)
    {
      MzMLFile mzml_file;
      MSExperiment map;
      mzml_file.setLogType(log_type_);
      mzml_file.load(mzml_file_path, map);

      ms_maps.push_back(map);
    }


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    progress_logger.startProgress(0, consensus_map.size(), "parsing features and ms2 identifications...");
    std::stringstream output_stream;
    Size feature_count = 1;
    for (Size i = 0; i != consensus_map.size(); ++i)
    {
      progress_logger.setProgress(i);
      // current feature
      const ConsensusFeature& feature = consensus_map[i];

      // store "mz rt" information from each scan
      stringstream scans_output;
      scans_output << setprecision(2) << setfill('0') << fixed;

      // determining charge and most intense feature for header
      BaseFeature::ChargeType charge = feature.getCharge();
      for (ConsensusFeature::HandleSetType::const_iterator feature_iter = feature.begin();
      feature_iter != feature.end(); ++feature_iter)
      {
        if (feature_iter->getCharge() > charge)
        {
          charge = feature_iter->getCharge();
        }
      }

      // print spectra information (PeptideIdentification tags)
      vector<PeptideIdentification> peptide_identifications = feature.getPeptideIdentifications();


      // clean peptide identifications outside mz rt tol ranges

      // vector of <<map index, spectrum index>, most intense ms2 scan>
      vector<pair<pair<double,PeptideIdentification>,pair<int,int>>> peptides;

      bool should_skip_feature;
      if (!(should_skip_feature = peptide_identifications.empty()))
      {
        for (Size peptide_index = 0; peptide_index < peptide_identifications.size(); peptide_index++)
        {
          auto peptide_identification = peptide_identifications[peptide_index];

          // append spectra information to scans_output
          int map_index = -1, spectrum_index = -1;
          if (peptide_identification.metaValueExists("spectrum_index"))
          {
            spectrum_index = peptide_identification.getMetaValue("spectrum_index");
          }
          if (peptide_identification.metaValueExists("map_index"))
          {
            map_index = peptide_identification.getMetaValue("map_index");
          }

          if (map_index != -1 && spectrum_index != -1)
          {
            // TEMP: log debug map index and spectrum index values once they are found
            LOG_DEBUG << "map index\t" << map_index << "\tspectrum index\t" << spectrum_index << endl;

            // retrieve spectrum for current peptide annotation
            auto ms2_scan = ms_maps[map_index][spectrum_index];
            ms2_scan.sortByIntensity(true);

            if (ms2_scan.getMSLevel() == 2 && !ms2_scan.empty())
            {
              should_skip_feature = false;

              // DEBUG determine if within user rt and mz tol range
              if (abs(feature.getMZ() - peptide_identification.getMZ()) > prec_mz_tol
              && abs(feature.getRT() - peptide_identification.getRT()) > prec_rt_tol)
              {
                continue;
              }

              double similarity_index = 5 * abs(feature.getMZ() - peptide_identification.getMZ()) +
              abs(feature.getRT() - peptide_identification.getRT());

              pair<double,PeptideIdentification> first_pair = pair<double,PeptideIdentification>(similarity_index, peptide_identification);
              pair<int,int> second_pair = pair<int,int>(map_index, spectrum_index);

              peptides.push_back(pair<pair<double,PeptideIdentification>,pair<int,int>>(first_pair,second_pair));
            }
          }
          else
          {
            should_skip_feature = true;
          }
        }
      }

      // peptides list of < <similarity_index, PeptideIdentification>, <map_index, feature_index> >

      // with the remaining peptides left within mz/rt tol of most intense
      if (!should_skip_feature && !peptides.empty())
      {
        // prepare peptides for output with highest mz value at top
        sort (peptides.begin(), peptides.end(), [](const pair<pair<double,PeptideIdentification>,pair<int,int>> &a, const pair<pair<double,PeptideIdentification>,pair<int,int>> &b)
        {
          return a.first.first < b.first.first;
        });

        // tmp stream for current feature
        stringstream feature_stream;
        feature_stream << setprecision(4) << fixed;

        // full spectra
        if (output_type == "full_spectra")
        {
          for (const auto& peptide : peptides)
          {
            feature_stream << "BEGIN IONS" << endl;

            feature_stream << "FEATURE_ID=" << to_string(feature_count) << endl;

            string filename = mzml_file_paths[peptide.second.first];
            Size parse_index = filename.rfind("/") + 1;
            filename = filename.substr(parse_index);
            feature_stream << "CONSENSUSID=e_" << feature.getUniqueId() << endl;

            feature_stream << "MSLEVEL=2" << endl;
            feature_stream << "CHARGE=" << to_string(charge == 0 ? 1 : charge) << "+" << endl;
            feature_stream << "PEPMASS=" << peptide.first.second.getMZ() << endl;
            feature_stream << "FILE_INDEX=" << peptide.second.second << endl;
            feature_stream << "RTINSECOND=" << peptide.first.second.getRT() << endl;

            auto ms2_scan = ms_maps[peptide.second.first][peptide.second.second];
            // sort spectra
            sort (ms2_scan.begin(), ms2_scan.end(), [](const Peak1D& a, const Peak1D& b)
            {
              return a.getMZ() > b.getMZ();
            });
            // ms2_scan.sortByIntensity(true);

            for (Size l = 0; l < ms2_scan.size(); l++)
            {
              feature_stream << ms2_scan[l].getMZ() << "\t" << (int) ms2_scan[l].getIntensity() << endl;
            }

            feature_stream << "END IONS" << endl << endl;
          }
          feature_count++;
        }
        else // merged spectra
        {
          // map mz to intensity
          map<double,int> ms2_block;

          // MSExperiment exp;

          const BinnedSpectrum binned_highest_int(ms_maps[peptides[0].second.first][peptides[0].second.second], BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

          for (const auto& peptide : peptides)
          {
            int map_index = peptide.second.first;
            int spectra_index = peptide.second.second;
            // int highest_binned_intensity = peptide.first.first;
            // auto highest_peptide_identification = peptide.first.second;

            auto spectrum = ms_maps[map_index][spectra_index];

            const BinnedSpectrum binned_spectrum(spectrum, BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

            BinnedSpectralContrastAngle bsca;
            double cosine_sim = bsca(binned_highest_int, binned_spectrum);
            // LOG_DEBUG << cosine_sim << " >= " << cos_sim << endl;

            // compare calculated cosine sim to binned highest int
            if (cosine_sim >= cos_sim)
            {
              for (Size spectrum_index = 0; spectrum_index < spectrum.size(); ++spectrum_index)
              {
                // exp.addSpectrum(spectrum);
                // exp.addSpectrum(spectrum);
                auto curr_spectrum = spectrum[spectrum_index];
                if (ms2_block[curr_spectrum.getMZ()] < curr_spectrum.getIntensity())
                {
                  ms2_block[curr_spectrum.getMZ()] = curr_spectrum.getIntensity();
                }
              }
            }
          }


          feature_stream << "BEGIN IONS" << endl;

          feature_stream << "FEATURE_ID=" << feature_count++ << endl;
          feature_stream << "CONSENSUSID=e_" << feature.getUniqueId() << endl;

          feature_stream << "MSLEVEL=2" << endl;
          feature_stream << "CHARGE=" << std::to_string(charge == 0 ? 1 : charge) << "+" << endl;
          feature_stream << "PEPMASS=" << peptides[0].first.second.getMZ() << endl;
          feature_stream << "FILE_INDEX=" << peptides[0].second.second << endl;
          feature_stream << "RTINSECOND=" << peptides[0].first.second.getRT() << endl;

          for (auto ms2_iter = ms2_block.rbegin(); ms2_iter != ms2_block.rend(); ++ms2_iter)
          {
            feature_stream << ms2_iter->first << "\t" << (int) ms2_iter->second << endl;
          }
          feature_stream << "END IONS" << endl << endl;
        }

        // output feature information to general outputStream
        output_stream << feature_stream.str() << endl;
      }
    }
    progress_logger.endProgress();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    ofstream output_file(out);
    progress_logger.startProgress(0, 1, "writing mgf file");
    output_file << output_stream.str();
    progress_logger.endProgress();
    output_file.close();

    return EXECUTION_OK;
  }
};

// the actual main functioned needed to create an executable
int main (int argc, const char** argv)
{
  TOPPGNPSExport tool;
  return tool.main(argc, argv);
}
/// @endcond
