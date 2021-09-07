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
// $Maintainer: Dorrestein Lab - University of California San Diego - https://dorresteinlab.ucsd.edu/$
// $Authors: Abinesh Sarvepalli and Louis Felix Nothias$
// $Contributors: Fabian Aicheler and Oliver Alka from Oliver Kohlbacher's group at Tubingen University$
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
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

//----------------------------------------------------------
// Doxygen docu
//----------------------------------------------------------
/**
  @page TOPP_GNPSExport GNPSExport

  @brief Export MS/MS data in .MGF format for GNPS (http://gnps.ucsd.edu).

GNPS (Global Natural Products Social Molecular Networking, http://gnps.ucsd.edu) is an open-access knowledge base for community-wide organization and sharing of raw, processed or identified tandem mass (MS/MS) spectrometry data. The GNPS web-platform makes possible to perform spectral library search against public MS/MS spectral libraries, as well as to perform various data analysis such as MS/MS molecular networking, network annotation propagation (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006089), and the Dereplicator-based annotation (https://www.nature.com/articles/nchembio.2219). The GNPS manuscript is available here: https://www.nature.com/articles/nbt.3597

This tool was developed for the Feature Based Molecular Networking (FBMN) workflow on GNPS (https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash2.jsp)

Please cite our preprint:
Nothias, LF., Petras, D., Schmid, R. et al. Feature-based molecular networking in the GNPS analysis environment.
Nat Methods 17, 905–908 (2020). https://doi.org/10.1038/s41592-020-0933-6

See the FBMN workflow documentation here (https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)

In brief, after running an OpenMS "metabolomics" pipeline, the GNPSExport TOPP tool can be used
on the consensusXML file and corresponding mzML files to generate the files needed for FBMN on GNPS.
These two files are:

	- The MS/MS spectral data file (.MGF format) which is generated  with the GNPSExport util.
	- The feature quantification table (.CSV format) which is generated with the TextExport util.

For each consensusElement in the consensusXML file, the GNPSExport produces one representative consensus
MS/MS spectrum (named peptide annotation in OpenMS jargon) outputted in the MS/MS spectral file (.MGF file).
Several modes for the generation of the consensus MS/MS spectrum are available and described below.
Note that these parameters are defined in the GNPSExport INI parameters file.

Representative command:
@code
GNPSExport -ini iniFile-GNPSExport.ini -in_cm filefilter.consensusXML -in_mzml inputFile0.mzML inputFile1.mzML -out GNPSExport_output.mgf
@endcode

The GNPSExport TOPP tool can be run on a consensusXML file and the corresponding mzML files to generate a MS/MS spectral file (MGF format)
and corresponding feature quantification table (.TXT format) that contains the LC-MS peak area intensity.

Requirements:
	- The IDMapper has to be run on the featureXML files, in order to associate MS2 scan(s) (peptide annotation) with each
	features. These peptide annotations are used by the GNPSExport.
	- The FileFilter has to be run on the consensusXML file, prior to the GNPSExport, in order to remove consensusElements
	without MS2 scans (peptide annotation).

Parameters:
	- Binning (ms2_bin_size): Defines the binning width of fragment ions during the merging of eligible MS/MS spectra.
	- Cosine Score Threshold (merged_spectra:cos_similarity): Defines the necessary pairwise cosine similarity with the highest precursor intensity MS/MS scan.

  - Output Type (output_type):
Options for outputting GNPSExport spectral processing are:
    -# [RECOMMENDED] merged_spectra
      For each consensusElement, the GNPSExport will merge all the eligible MS/MS scans into one representative consensus MS/MS spectrum.
      Eligible MS/MS scans have a pairwise cosine similarity with the MS/MS scan of highest precursor intensity above the Cosine Similarity Threshold.
	    The fragment ions of merged MS/MS scans are binned in m/z (or Da) range defined by the Binning width parameter.
      .
	  -# Most intense: most_intense - For each consensusElement, the GNPSExport will output the most intense MS/MS scan (with the highest precursor ion intensity) as consensus MS/MS spectrum.
      .

Note that mass accuracy and the retention time window for the pairing between MS/MS scans and a LC-MS feature
or consensusElement is defined at the IDMapper tool step.

A representative OpenMS-GNPS workflow would sequentially use these OpenMS TOPP tools:
  1. Input mzML files
  2. Run the @ref TOPP_FeatureFinderMetabo tool on the mzML files.
  3. Run the @ref TOPP_IDMapper tool on the featureXML and mzML files.
  4. Run the @ref TOPP_MapAlignerPoseClustering tool on the featureXML files.
  5. Run the @ref TOPP_MetaboliteAdductDecharger on the featureXML files.
  6. Run the @ref TOPP_FeatureLinkerUnlabeledKD tool or FeatureLinkerUnlabeledQT, on the featureXML files and output a consensusXML file.
  8. Run the @ref TOPP_FileFilter on the consensusXML file to keep only consensusElements with at least MS/MS scan (peptide identification).
  9. Run the @ref TOPP_GNPSExport on the "filtered consensusXML file" to export an .MGF file.
  10. Run the @ref TOPP_TextExporter on the "filtered consensusXML file" to export an .TXT file.
  11. Upload your files to GNPS and run the Feature-Based Molecular Networking workflow. Instructions are here:
https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/

The GitHub for that ProteoSAFe workflow and an OpenMS python wrappers is available here:
https://github.com/Bioinformatic-squad-DorresteinLab/openms-gnps-workflow

An online version of the OpenMS-GNPS pipeline for FBMN running on CCMS server (http://proteomics.ucsd.edu/) is available on GNPS:
https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-openms/

GNPS (Global Natural Products Social Molecular Networking, https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash2.jsp)
is an open-access knowledge base for community-wide organization and sharing of raw, processed
or identified tandem mass (MS/MS) spectrometry data.
The GNPS web-platform makes possible to perform spectral library search against public MS/MS spectral libraries,
as well as to perform various data analysis such as MS/MS molecular networking, Network Annotation Propagation
Network Annotation Propagation (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006089)
and the DEREPLICATOR (https://www.nature.com/articles/nchembio.2219)
The GNPS paper is available here (https://www.nature.com/articles/nbt.3597)

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_GNPSExport.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_GNPSExport.html
 */

class TOPPGNPSExport : public TOPPBase
{
public:
  TOPPGNPSExport() : TOPPBase(
    "GNPSExport",
    "Tool to export representative consensus MS/MS scan per consensusElement into a .MGF file format.\nSee the documentation on https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking_with_openms",
    true,
    {
      {
        "Nothias L.F. et al.", // authors
        "Feature-based Molecular Networking in the GNPS Analysis Environment", // title
        "bioRxiv 812404 (2019)", // when_where
        "10.1101/812404" // doi
      }
    }
  ) {}

private:
  static constexpr double DEF_COSINE_SIMILARITY = 0.9;
  static constexpr double DEF_MERGE_BIN_SIZE = static_cast<double>(BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES);

  static constexpr double DEF_PREC_MASS_TOL = 0.5;
  static constexpr bool DEF_PREC_MASS_TOL_ISPPM = false;

  static constexpr double DEF_PEPT_CUTOFF = 5;
  static constexpr double DEF_MSMAP_CACHE = 50;

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
    double sum_mz = 0, sum_intensity = 0;
    int count = 0;
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
    for (unsigned long i = 0; i < exp.getSpectra().size(); ++i)
    {
      MSExperiment::SpectrumType &spec = exp.getSpectrum(i);
      for (auto& spec : spec)
      {
        Peak1D curr(spec.getMZ(), spec.getIntensity());
        flat_spectra.push_back(curr);
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
                  << "CHARGE=" << to_string(feature_charge == 0 ? 1 : feature_charge) << "+" << "\n"
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
    for (ConsensusFeature::HandleSetType::const_iterator feature_iter = feature.begin();\
          feature_iter != feature.end(); ++feature_iter)
    {
      element_maps.push_back(pair<int,double>(feature_iter->getMapIndex(), feature_iter->getIntensity()));
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
          pepts.push_back(pair<int,int>(map_index,spec_index));
          break;
        }
      }
    }

    // return will be reformatted vector<PeptideIdentification> pepts passed in by value
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_cm", "<file>", "", "Input consensusXML file containing only consensusElements with \"peptide\" annotations.");
    setValidFormats_("in_cm", ListUtils::create<String>("consensusXML"));

    registerInputFileList_("in_mzml", "<files>", ListUtils::create<String>(""), "Original mzml files containing the ms2 spectra (aka peptide annotation). \nMust be in order that the consensusXML file maps the original mzML files.");
    setValidFormats_("in_mzml", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Output MGF file");
    setValidFormats_("out", ListUtils::create<String>("mgf"));

    registerStringOption_("output_type", "<choice>", "most_intense", "specificity of mgf output information", false);
    setValidStrings_("output_type", ListUtils::create<String>("merged_spectra,most_intense"));

    addEmptyLine_();

    // registerIntOption_("msmap_cache", "<num>", DEF_MSMAP_CACHE, "Number of msmaps that can be cached during export for optimized performance", false, true);
    // setMinInt_("msmap_cache", 0);
    registerIntOption_("peptide_cutoff", "<num>", DEF_PEPT_CUTOFF, "Number of most intense peptides to consider per consensus element; '-1' to consider all identifications", false, true);
    setMinInt_("peptide_cutoff", -1);
    registerDoubleOption_("ms2_bin_size", "<num>", DEF_MERGE_BIN_SIZE, "Bin size (Da) for fragment ions when merging ms2 scans", false, false);
    setMinFloat_("ms2_bin_size", 0);

    registerTOPPSubsection_("merged_spectra", "Options for exporting mgf file with merged spectra per consensusElement");
    registerDoubleOption_("merged_spectra:precursor_mass_tolerance", "<num>", DEF_PREC_MASS_TOL, "Precursor mass tolerance (Da) for ms annotations", false);
    setMinFloat_("merged_spectra:precursor_mass_tolerance", 0);
    registerDoubleOption_("merged_spectra:cos_similarity", "<num>", DEF_COSINE_SIMILARITY, "Cosine similarity threshold for merged_spectra output", false);
    setMinFloat_("merged_spectra:cos_similarity", 0);
  }

  // the main function is called after all parameters are read
  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    double cos_sim_threshold(getDoubleOption_("merged_spectra:cos_similarity"));
    double bin_width(getDoubleOption_("ms2_bin_size"));

    String consensus_file_path(getStringOption_("in_cm"));
    StringList mzml_file_paths = getStringList_("in_mzml");
    String out(getStringOption_("out"));
    String output_type(getStringOption_("output_type"));

    int pept_cutoff((output_type == "merged_spectra") ? getIntOption_("peptide_cutoff") : 1);

    ofstream output_file(out);

    ProgressLogger progress_logger;
    progress_logger.setLogType(log_type_);


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    // ConsensusMap
    ConsensusXMLFile consensus_file;
    ConsensusMap consensus_map;
    consensus_file.load(consensus_file_path, consensus_map);


    //-------------------------------------------------------------
    // preprocessing: allocate memory
    //-------------------------------------------------------------
    // max_msmap_cache = std::min(max_msmap_cache, static_cast<int>(mzml_file_paths.size()));
    int max_msmap_cache = static_cast<int>(mzml_file_paths.size());
    MzMLFile *mzml_files = new MzMLFile[max_msmap_cache];
    MSExperiment *specs_list = new MSExperiment[max_msmap_cache];

    map<int,int> msmaps_cached; // <K, V> = <map_index, mzml_files index>
    Size num_msmaps_cached = 0;


    //-------------------------------------------------------------
    // write output (+ merge computations)
    //-------------------------------------------------------------
    progress_logger.startProgress(0, consensus_map.size(), "parsing features and ms2 identifications...");
    for (Size cons_i = 0; cons_i < consensus_map.size(); ++cons_i)
    {
      progress_logger.setProgress(cons_i);

      const ConsensusFeature &feature = consensus_map[cons_i];

      // determine feature's charge
      BaseFeature::ChargeType charge = feature.getCharge();
      for (ConsensusFeature::HandleSetType::const_iterator feature_iter = feature.begin();\
            feature_iter != feature.end(); ++feature_iter)
      {
        if (feature_iter->getCharge() > charge)
        {
          charge = feature_iter->getCharge();
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
      for (vector<pair<int,int>>::iterator pepts_iter=pepts.begin(); pepts_iter!=pepts.end(); pepts_iter++)
      {
        int map_index = pepts_iter->first;

        // load missing MzMLFile files + MSExperiment specs
        if (msmaps_cached.find(map_index) == msmaps_cached.end())
        {
          MzMLFile mzmlfile;
          mzml_files[num_msmaps_cached] = mzmlfile;

          MSExperiment exp;
          specs_list[num_msmaps_cached] = exp;
          mzml_files[num_msmaps_cached].load(mzml_file_paths[map_index], specs_list[num_msmaps_cached]);

          msmaps_cached[map_index] = num_msmaps_cached;
          num_msmaps_cached += 1;
        }
      }

      // identify most intense spectrum
      const int best_mapi = (pepts[0]).first;
      const int best_speci = (pepts[0]).second;
      auto best_spec = specs_list[msmaps_cached[best_mapi]][best_speci];

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
          auto test_spec = specs_list[msmaps_cached[map_index]][spec_index];

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

    delete [] mzml_files;
    delete [] specs_list;

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
