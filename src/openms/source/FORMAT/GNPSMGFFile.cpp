// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Dorrestein Lab - University of California San Diego - https://dorresteinlab.ucsd.edu/$
// $Authors: Abinesh Sarvepalli and Louis Felix Nothias$
// $Contributors: Fabian Aicheler and Oliver Alka from Oliver Kohlbacher's group at Tubingen University$
// --------------------------------------------------------------------------
//

#include <OpenMS/FORMAT/GNPSMGFFile.h>
#include <OpenMS/COMPARISON/BinnedSpectralContrastAngle.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS
{
  GNPSMGFFile::GNPSMGFFile() :
    DefaultParamHandler("GNPSMGFFile"),
    ProgressLogger()
  {
    defaults_.setValue("output_type", "most_intense", "specificity of mgf output information");
    defaults_.setValidStrings("output_type", {"merged_spectra","most_intense"});

    defaults_.setValue("peptide_cutoff", DEF_PEPT_CUTOFF, "Number of most intense peptides to consider per consensus element; '-1' to consider all identifications.");
    defaults_.setMinInt("peptide_cutoff", -1);

    defaults_.setValue("ms2_bin_size", DEF_MERGE_BIN_SIZE, "Bin size (Da) for fragment ions when merging ms2 scans.");
    defaults_.setMinFloat("ms2_bin_size", 0);

    defaults_.setValue("merged_spectra:cos_similarity", DEF_COSINE_SIMILARITY, "Cosine similarity threshold for merged_spectra output.");
    defaults_.setMinFloat("merged_spectra:cos_similarity", 0);

    defaults_.setSectionDescription("merged_spectra", "Options for exporting mgf file with merged spectra per consensusElement");

    defaultsToParam_(); // copy defaults into param_
  }

  /**
   * @brief Bin peaks by similar m/z position and averaged intensities
   * @param peaks Vector of Peak1D peaks sorted by m/z position
   * @param bin_width Size of bin
   * @param binned_peaks Result vector with binned peaks passed in by reference
   */
  void binPeaks_(
    const vector<Peak1D> &peaks,
    const double bin_width,
    vector<Peak1D> &binned_peaks
  )
  {
    double last_mz = peaks.at(0).getMZ();
    double sum_mz{};
    double sum_intensity{};
    int count{};
    for (auto& spec : peaks)
    {
      if (count > 0 && spec.getMZ() - last_mz > bin_width)
      {
        if (sum_intensity > 0)
        {
          Peak1D curr(sum_mz/count, sum_intensity/count);
          binned_peaks.push_back(curr);
        }
        last_mz = spec.getMZ();
        sum_mz = 0;
        sum_intensity = 0;
        count = 0;
      }

      sum_mz += spec.getMZ();
      sum_intensity += spec.getIntensity();
      count += 1;
    }
    if (count > 0)
    {
      Peak1D curr(sum_mz/count, sum_intensity/count);
      binned_peaks.push_back(curr);
    }
  }

  /**
   * @brief Flatten spectra from MSExperiment into a single vector of Peak1D peaks
   * @param exp MSExperiment containing at least 1 spectrum
   * @param bin_width Size of binned scan (m/z)
   * @param merged_peaks Result vector of peaks passed in by reference
   */
  void flattenAndBinSpectra_(
    MSExperiment &exp,
    const double bin_width,
    vector<Peak1D> &merged_peaks
  )
  {
    // flatten spectra
    vector<Peak1D> flat_spectra;
    for (auto& spec : exp.getSpectra())
    {
      for (auto& peak : spec)
      {
        flat_spectra.push_back(peak);
      }
    }

    sort(flat_spectra.begin(), flat_spectra.end(), [](const Peak1D &a, const Peak1D &b)
      {
        return a.getMZ() < b.getMZ();
      }
    );

    // bin peaks
    binPeaks_(flat_spectra, bin_width, merged_peaks);

    // return value is modified merged_peaks passed by reference
  }

  /**
   * @brief Private function that outputs MS/MS Block Header
   * @param output_file Stream that will write to file
   * @param scan_index Current scan index in GNPSExport formatted output
   * @param feature_id ConsensusFeature Id found in input mzXML file
   * @param feature_charge ConsensusFeature's highest charge as mentioned in the input mzXML file
   * @param feature_mz m/z position of PeptideIdentification with highest intensity
   * @param spec_index Spectrum index of PeptideIdentification with highest intensity
   * @param feature_rt ConsensusFeature's retention time specified in input mzXML file
   */
  void writeMSMSBlockHeader_(
    ofstream &output_file,
    const String &output_type,
    const int &scan_index,
    const String &feature_id,
    const int &feature_charge,
    const String &feature_mz,
    const String &spec_index,
    const String &feature_rt
  )
  {
    if (output_file.is_open())
    {
      output_file << "BEGIN IONS" << "\n"
                  << "OUTPUT=" << output_type << "\n"
                  << "SCANS=" << scan_index << "\n"
                  << "FEATURE_ID=e_" << feature_id << "\n"
                  << "MSLEVEL=2" << "\n"
                  << "CHARGE=" << to_string(feature_charge == 0 ? 1 : abs(feature_charge))+(feature_charge >= 0 ? "+" : "-") << "\n"
                  << "PEPMASS=" << feature_mz << "\n"
                  << "FILE_INDEX=" << spec_index << "\n"
                  << "RTINSECONDS=" << feature_rt << "\n";
    }
  }

  /**
   * @brief Private function to write peak mass and intensity to output file
   * @param output_file Stream that will write to file
   * @param peaks Vector of peaks that will be outputted
   */
  void writeMSMSBlock_(
    ofstream &output_file,
    const vector<Peak1D> &peaks
  )
  {
    if (output_file.is_open())
    {
      output_file << setprecision(4) << fixed;
      for (auto& peak : peaks)
      {
        output_file << peak.getMZ() << "\t" << peak.getIntensity() << "\n";
      }

      output_file << "END IONS" << "\n" << endl;
    }
  }

  /**
   * @brief Private method used to sort PeptideIdentification map indices in order of annotation's intensity
   * @param feature ConsensusFeature annotated with PeptideIdentifications
   * @param featureMaps_sortedByInt Result vector of map indices in order of PeptideIdentification intensity
   */
  void sortElementMapsByIntensity_(const ConsensusFeature& feature, vector<pair<int,double>>& element_maps)
  {
    // convert element maps to vector of pair<int,double>(map, intensity)
    for (const auto& feature_handle : feature)
    {
      element_maps.emplace_back(feature_handle.getMapIndex(), feature_handle.getIntensity());
    }

    // sort elements by intensity
    sort(element_maps.begin(), element_maps.end(), [](const pair<int,double> &a, const pair<int,double> &b)
    {
      return a.second > b.second;
    });

    // return value will be reformatted vector<int> element_maps passed in by value
  }

  /**
   * @brief Retrieve list of PeptideIdentification parameters from ConsensusFeature metadata, sorted by map intensity
   * @param feature ConsensusFeature feature containing PeptideIdentification annotations
   * @param sorted_element_maps Sorted list of element_maps
   * @param pepts Result vector of <map_index,spectrum_index> of PeptideIdentification annotations sorted by map intensity in feature
   */
  void getElementPeptideIdentificationsByElementIntensity_(
    const ConsensusFeature& feature,
    vector<pair<int,double>>& sorted_element_maps,
    vector<pair<int,int>>& pepts
  )
  {
    for (pair<int,double>& element_pair : sorted_element_maps)
    {
      int element_map = element_pair.first;
      vector<PeptideIdentification> feature_pepts = feature.getPeptideIdentifications();
      for (PeptideIdentification& pept_id : feature_pepts)
      {
        if (pept_id.metaValueExists("spectrum_index") && pept_id.metaValueExists("map_index")
            && (int)pept_id.getMetaValue("map_index") == element_map)
        {
          int map_index = pept_id.getMetaValue("map_index");
          int spec_index = pept_id.getMetaValue("spectrum_index");
          pepts.emplace_back(map_index,spec_index);
          break;
        }
      }
    }
    // return will be reformatted vector<PeptideIdentification> pepts passed in by value
  }

  void GNPSMGFFile::store(const String& consensus_file_path, const StringList& mzml_file_paths, const String& out) const
  {
    std::string output_type = getParameters().getValue("output_type");

    double bin_width = getParameters().getValue("ms2_bin_size");

    int pept_cutoff((output_type == "merged_spectra") ? (int)getParameters().getValue("peptide_cutoff") : 1);

    double cos_sim_threshold = getParameters().getValue("merged_spectra:cos_similarity");

    ofstream output_file(out);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    // ConsensusMap
    ConsensusMap consensus_map;
    FileHandler().loadConsensusFeatures(consensus_file_path, consensus_map, {FileTypes::CONSENSUSXML});

    //-------------------------------------------------------------
    // open on-disc data (=spectra are only loaded on demand to safe memory)
    //-------------------------------------------------------------
    vector<OnDiscMSExperiment> specs_list(mzml_file_paths.size(), OnDiscMSExperiment());

    map<size_t, size_t> map_index2file_index; // <K, V> = <map_index, file_index>
    Size num_msmaps_cached = 0;

    //-------------------------------------------------------------
    // write output (+ merge computations)
    //-------------------------------------------------------------
    startProgress(0, consensus_map.size(), "parsing features and ms2 identifications...");

    for (Size cons_i = 0; cons_i < consensus_map.size(); ++cons_i)
    {
      setProgress(cons_i);

      const ConsensusFeature& feature = consensus_map[cons_i];

      // determine feature's charge as maximum feature handle charge
      int charge = feature.getCharge();
      for (auto& fh : feature)
      { 
        if (fh.getCharge() > charge)
        {
          charge = fh.getCharge();
        }
      }

      // compute most intense peptide identifications (based on precursor intensity)
      vector<pair<int,double>> element_maps;
      sortElementMapsByIntensity_(feature, element_maps);
      vector<pair<int, int>> pepts;
      getElementPeptideIdentificationsByElementIntensity_(feature, element_maps, pepts);

      // discard poorer precursor spectra for 'merged_spectra' and 'full_spectra' output
      if (pept_cutoff != -1 && pepts.size() > (unsigned long) pept_cutoff)
      {
        pepts.erase(pepts.begin()+pept_cutoff, pepts.end());
      }

      // validate all peptide annotation maps have been loaded
      for (const auto& pep : pepts)
      {
        int map_index = pep.first;

        // open on-disc experiments
        if (map_index2file_index.find(map_index) == map_index2file_index.end())
        {
          specs_list[num_msmaps_cached].openFile(mzml_file_paths[map_index], false); // open on-disc experiment and load meta-data
          map_index2file_index[map_index] = num_msmaps_cached;
          ++num_msmaps_cached;
        }
      }

      // identify most intense spectrum
      const int best_mapi = pepts[0].first;
      const int best_speci = pepts[0].second;
      auto best_spec = specs_list[map_index2file_index[best_mapi]][best_speci];

      if (best_spec.empty()) continue; // some Bruker files have MS2 spectra without peaks. skip those during exprot

      // write block output header
      writeMSMSBlockHeader_(
        output_file,
        output_type,
        (cons_i + 1),
        feature.getUniqueId(),
        charge,
        feature.getMZ(),
        best_speci,
        best_spec.getRT()
      );

      // OPENMS_LOG_DEBUG << "Best spectrum (index/RT): " << best_speci << "\t" << best_spec.getRT() << std::endl;

      // store outputted spectra in MSExperiment
      MSExperiment exp;

      // add most intense spectrum to MSExperiment
      exp.addSpectrum(best_spec);

      if (output_type == "merged_spectra")
      {
        // merge spectra that meet cosine similarity threshold to most intense spectrum
        BinnedSpectrum binned_highest_int(best_spec, bin_width, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

        // Retain peptide annotations that do not meet user-specified cosine similarity threshold
        for (pair<int,int> &pept : pepts)
        {
          int map_index = pept.first;
          int spec_index = pept.second;
          auto test_spec = specs_list[map_index2file_index[map_index]][spec_index];

          BinnedSpectrum binned_spectrum(test_spec, bin_width, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

          BinnedSpectralContrastAngle bsca;
          double cos_sim = bsca(binned_highest_int, binned_spectrum);

          if (cos_sim >= cos_sim_threshold)
          {
            exp.addSpectrum(test_spec);
          }
        }
      }

      // store outputted peaks in vector<Peak1D>
      vector<Peak1D> peaks;
      flattenAndBinSpectra_(
        exp,
        bin_width,
        peaks
      );

      // write peaks to output block
      writeMSMSBlock_(
        output_file,
        peaks
      );
    }

    output_file.close();
  }
}
